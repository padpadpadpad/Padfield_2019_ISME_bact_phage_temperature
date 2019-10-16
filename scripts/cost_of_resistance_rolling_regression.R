#--------------------------------------------------------------------------#
# analysis of resistant and susceptible clones using rolling regression ####
#--------------------------------------------------------------------------#

# clear workspace
rm(list = ls())

# load in packages
library(nls.multstart)
library(tidyverse)
library(broom)
library(modelr)
library(nlsMicrobio)
library(brms)
library(tidybayes)
library(patchwork)
library(zoo)

# read in data
d <- readRDS('data/cost_of_resistance.rds')

# load in functions
source('scripts/growth_curve_models.R')

#-----------------#
# wrangle data ####
#-----------------#

# if od_cor is negative or less than the blank, or at the minimum possible difference between control and sample (0.001), replace with NA
# try different od_cors. Here use the min, control3
d <- mutate(d, od_cor = od - control,
            od_cor2 = ifelse(od_cor <= 0.001, NA, od_cor))

# create log abundance
d <- mutate(d,
            log_od_cor = log(od_cor2),
            # parse timepoint number,
            time_point = readr::parse_number(time_point),
            time = ymd_hms(as.character(time)))

# add time in hours
d <- group_by(d, temp, phage, well) %>%
  arrange(., time_point) %>%
  mutate(., lag = lag(time, order_by = time_point), 
         t = as.numeric(time - lag(time, order_by = time_point), units = 'hours'),
         t = ifelse(is.na(t), 0, t),
         t = cumsum(t)) %>%
  ungroup() %>%
  select(-lag) %>%
  mutate(., id = group_indices(., temp, phage),
         log_od_cor = ifelse(is.nan(log_od_cor), NA, log_od_cor))

# remove NAs and infinites
d2 <- filter(d, !is.na(log_od_cor) & !is.infinite(log_od_cor))

# delete the first time point - prone to bubbles and problems of being fixed
d2 <- filter(d2, time_point >= 1)

#-----------------------------------------------------------------#
# run a rolling regression for each sample every 4 time points ####
#-----------------------------------------------------------------#

# define function for rolling regression
rolling_coefs <- . %>% data.frame %>% lm() %>% broom::tidy()

# run rolling regression on ln od_cor ~ time
# every 4 time points
models <- d2 %>%
  group_by(temp, clone, rep, phage) %>%
  do(cbind(model = select(., log_od_cor, t) %>% 
           zoo::rollapplyr(width = 4, rolling_coefs, by.column = FALSE, fill = NA),
           t = select(., t)))

# for each curve, calculate maximum slope of rolling regression - this represents the exponential growth rate
d_growth_rate <- filter(models, model.term == 't') %>%
  mutate_at(., vars(starts_with('model')), as.numeric) %>%
  group_by(temp, clone, rep, phage) %>%
  filter(., model.estimate == max(model.estimate, na.rm = TRUE)) %>%
  ungroup()
  
# clean up data set
d_growth_rate <- janitor::clean_names(d_growth_rate) %>%
  mutate(., rep2 = group_indices(., rep, phage),
         K = temp + 273.15)

# calculate mean growth rate, average over clones within replicate
d_means <- group_by(d_growth_rate, rep2, phage, temp, K) %>%
  summarise(gr_rate = mean(model_estimate),
            sd = sd(model_estimate)) %>%
  ungroup()

#-----------------------------------------------#
# model thermal performance curve using brms ####
#-----------------------------------------------#

# define the model - uses the Sharpe-Schoolfield model
nlform <- bf(gr_rate ~  c*exp(E/(8.62 * 10^-5) * (1/293.15 - 1/K)) * (1/(1 + exp(Eh/(8.62 * 10^-5) * (1/Th - 1/K)))), 
             E ~ 0 + phage + (1|rep2), 
             c ~ 0 + phage + (1|rep2), 
             Eh ~ 0 + phage + (1|rep2), 
             Th ~ 0 + phage + (1|rep2), 
             nl = TRUE)

# set up priors
nlprior <- c(prior(normal(1, 5), nlpar = "c"), 
             prior(normal(0.6, 5), nlpar = "E", lb = 0), 
             prior(normal(4, 5), nlpar = "Eh",lb = 0, ub = 10), 
             prior(normal(300, 5), nlpar = "Th", lb = 273.15))

# run model
fit_bayes <- brm(formula = nlform, data = d_means, 
                  family = gaussian(), 
                  prior = nlprior, 
                  control = list(adapt_delta = 0.9), 
                  chains = 3, 
                  iter = 5000)

# look at model
fit_bayes

# extract and clean up TPC "traits" 
d_params_bact <- fit_bayes %>%
  spread_draws(`b_.*`, regex = TRUE) %>%
  mutate(., b_Topt_phageB = get_topt(b_Eh_phageB, b_Th_phageB, b_E_phageB),
         b_Topt_phageBP = get_topt(b_Eh_phageBP, b_Th_phageBP, b_E_phageBP),
         b_rmax_phageB = sharpeschoolhigh_1981(b_Topt_phageB, b_c_phageB, b_E_phageB, b_Eh_phageB, b_Th_phageB, 20),
         b_rmax_phageBP = sharpeschoolhigh_1981(b_Topt_phageBP, b_c_phageBP, b_E_phageBP, b_Eh_phageBP, b_Th_phageBP, 20),
         b_relrmax_rmax = 1-b_rmax_phageBP/b_rmax_phageB) %>%
  gather(., 'param', 'estimate', 4:ncol(.)) %>%
  separate(., param, c('blah', 'term', 'treat'), sep = '_') %>%
  select(., -starts_with('blah')) %>%
  mutate(., treat = substr(treat, 6, nchar(treat)),
         estimate = ifelse(term %in% c('Th', 'Topt'), estimate - 273.15, estimate),
         treat = ifelse(treat == 'B', 'bact', 'bact_phage')) %>%
  group_by(., term, treat) %>%
  mean_qi()

# plot TPC "traits"
fig_params <- filter(d_params_bact, term %in% c('E', 'Eh', 'Topt', 'c', 'rmax')) %>%
  ggplot(., aes(treat, estimate, col = treat)) + 
  geom_point(size = 3) + 
  geom_linerange(aes(ymin = .lower, ymax = .upper)) +
  facet_wrap(~term, scale = "free", ncol = 1) + 
  scale_color_manual(values = c("black", "blue")) +
  theme_bw(base_size = 14) +
  theme(legend.position = 'none') +
  NULL

# get predictions over the treatment level, ignoring the random effects
preds_bact <- modelr::data_grid(d_means, K = modelr::seq_range(K, n = 100), phage) %>%
  add_fitted_draws(fit_bayes, re_formula = NA) %>%
  data.frame() %>%
  group_by(K, phage) %>%
  mean_qi(estimate = .value)

# get difference in predictions between resistant and susceptible clones
preds_diff <- modelr::data_grid(d_growth_rate, K = modelr::seq_range(K, n = 100), phage) %>%
  add_fitted_draws(fit_bayes, re_formula = NA) %>%
  data.frame() %>%
  select(., K, phage, .value, .draw) %>%
  spread(., phage, .value) %>%
  mutate(., .value = BP / B) %>%
  group_by(K) %>%
  mean_qi(estimate = .value)

# plot TPCs
plot_bact <- ggplot() +
  geom_point(aes(K - 273.15, gr_rate, col = phage), d_means, alpha = 0.5, size = 3) +
  geom_line(aes(K - 273.15, estimate, col = phage, group = interaction(phage)), data = preds_bact) +
  geom_ribbon(aes(K - 273.15, ymin = .lower, ymax = .upper, fill = phage), data = preds_bact, alpha = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c('blue', 'black')) +
  scale_color_manual(values = c('blue', 'black')) +
  xlab('Temperature (ºC)') +
  ylab(expression('Bacteria growth rate ( hr'^-1~')')) +
  ggtitle('(a)') +
  theme(legend.position = 'none')

# plot difference between TPcs
plot_diff <- ggplot(preds_diff) +
  geom_line(aes(K - 273.15, estimate), data = preds_diff) +
  geom_ribbon(aes(K - 273.15, ymin = .lower, ymax = .upper), data = preds_diff, alpha = 0.2) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  theme_bw(base_size = 14) +
  xlab('Temperature (ºC)') +
  ylab('Relative fitness (resistant / susceptible bacteria)') +
  ggtitle('(b)')

p1 <- plot_bact + plot_diff

ggsave('plots/Figure_4.png', p1, width = 13, height = 6)
