#--------------------------------------------------#
# Bayesian model fits of TPCs of phage replication #
# and bacteria in the presence/absence of phage    #
#--------------------------------------------------#

# load in packages ####
library(tidyverse)
library(janitor)
library(patchwork)
library(brms)
library(tidybayes)

source('scripts/growth_curve_models.R')

#--------------------------#
# Phage replication TPC ####
#--------------------------#

# load in data ####
d <- read.csv('data/phage_replication.csv', stringsAsFactors = FALSE) %>%
  mutate(., abund_ml = MicrobioUoE::calc_PFU_ml(count, dil_fac, vol_phage = vol_plated),
         logabund = log10(abund_ml))

# wrangle data ####
d_inoc <- filter(d, temp == 'inoculate') %>%
  # convert number of phage into amount put into flasks. 20 µl into 6 ml
  mutate(., abund_ml = abund_ml / 50 / 6) %>%
  summarise(., inoc_density = mean(abund_ml)) %>%
  mutate(., logabund = log(inoc_density))

# calculate abundance for each temperature
d <- filter(d, temp != 'inoculate') %>%
  mutate(., temp = as.numeric(temp)) %>%
  group_by(temp, rep, dil_fac, vol_plated) %>%
  summarise(., count = mean(count)) %>%
  ungroup() %>%
  mutate(abund_ml = MicrobioUoE::calc_PFU_ml(count, dil_fac, vol_phage = vol_plated),
         logabund = log(abund_ml),
         K = temp + 273.15,
         logabund = ifelse(is.infinite(logabund), NA, logabund),
         inoc_density = d_inoc$inoc_density,
         log_inoc = d_inoc$logabund,
         growth_rate = (logabund - log_inoc) / 4,
         # add buffer for modelling
         gr_cor = growth_rate + 3)

# fit phage model using brms
nlform <- bf(growth_rate ~ -3 + c*exp(E/(8.62 * 10^-5) * (1/293.15 - 1/K)) * (1/(1 + exp(Eh/(8.62 * 10^-5) * (1/Th - 1/K)))), 
             E ~ 1, 
             c ~ 1, 
             Eh ~ 1, 
             Th ~ 1, 
             nl = TRUE)

# set up priors
nlprior <- c(prior(normal(3, 5), nlpar = "c"), 
             prior(normal(0.6, 5), nlpar = "E", lb = 0), 
             prior(normal(4, 5), nlpar = "Eh",lb = 0, ub = 25), 
             prior(normal(300, 5), nlpar = "Th", lb = 273.15))

# run model
# only run on first time there are negative growth rates - 30 ºC
d2 <- filter(d, temp <= 30)

fit_bayes_phage <- brm(formula = nlform, data = d2, 
                  family = gaussian(), 
                  prior = nlprior, 
                  control = list(adapt_delta = 0.9), 
                  chains = 3, 
                  iter = 5000)

summary(fit_bayes_phage)

#------------------------------------#
# Bacteria with and without phage ####
#------------------------------------#

d_bact <- readRDS('data/bacteria_with_without_phage_growth_rates.rds') %>%
  mutate(., K = temp + 273.15) %>%
  rename(., r = estimate)
d_bact_raw <- readRDS('data/bacteria_with_without_phage_growth_rates_raw.rds') %>%
  mutate(., K = temp + 273.15) %>%
  rename(., r = estimate)

# fit bacteria model using brms
nlform <- bf(r ~  c*exp(E/(8.62 * 10^-5) * (1/293.15 - 1/K)) * (1/(1 + exp(Eh/(8.62 * 10^-5) * (1/Th - 1/K)))), 
             E ~ 0 + phage, 
             c ~ 0 + phage, 
             Eh ~ 0 + phage, 
             Th ~ 0 + phage, 
             nl = TRUE)

# set up priors
nlprior <- c(prior(normal(1, 5), nlpar = "c"), 
             prior(normal(0.6, 5), nlpar = "E", lb = 0), 
             prior(normal(4, 5), nlpar = "Eh",lb = 0, ub = 25), 
             prior(normal(300, 5), nlpar = "Th", lb = 273.15))

# run model on cleaned and raw datasets
fit_bayes <- brm(formula = nlform, data = d_bact, 
                  family = gaussian(), 
                  prior = nlprior, 
                  control = list(adapt_delta = 0.9), 
                  chains = 3, 
                  iter = 5000)

fit_bayes_raw <- brm(formula = nlform, data = d_bact_raw, 
                 family = gaussian(), 
                 prior = nlprior, 
                 control = list(adapt_delta = 0.9), 
                 chains = 3, 
                 iter = 5000)

#---------------------------------------------------#
# extract parameters and predictions from models ####
#---------------------------------------------------#

# parameters from phage model
d_params_phage <- fit_bayes_phage %>%
  spread_draws(`b_.*`, regex = TRUE) %>%
  mutate(., b_Topt_Intercept = get_topt(b_Eh_Intercept, b_Th_Intercept, b_E_Intercept)) %>%
  gather(., 'param', 'estimate', 4:ncol(.)) %>%
  separate(., param, c('blah', 'term', 'blah2'), sep = '_') %>%
  select(., -starts_with('blah')) %>%
  mutate(., estimate = ifelse(term %in% c('Th', 'Topt'), estimate - 273.15, estimate),
         treat = 'phage') %>%
  filter(., !is.nan(estimate)) %>%
  group_by(., term, treat) %>%
  mean_qi()

# parameters from bacteria model
d_params_bact <- fit_bayes %>%
  spread_draws(`b_.*`, regex = TRUE) %>%
  mutate(., b_Topt_phagenophage = get_topt(b_Eh_phagenophage, b_Th_phagenophage, b_E_phagenophage),
         b_Topt_phagephage = get_topt(b_Eh_phagephage, b_Th_phagephage, b_E_phagephage),
         b_rmax_phagenophage = sharpeschoolhigh_1981(b_Topt_phagenophage, b_c_phagenophage, b_E_phagenophage, b_Eh_phagenophage, b_Th_phagenophage, 20),
         b_rmax_phagephage = sharpeschoolhigh_1981(b_Topt_phagephage, b_c_phagephage, b_E_phagephage, b_Eh_phagephage, b_Th_phagephage, 20),
         b_relrmax_rmax = 1 - b_rmax_phagephage/b_rmax_phagenophage) %>%
  gather(., 'param', 'estimate', 4:ncol(.)) %>%
  separate(., param, c('blah', 'term', 'treat'), sep = '_') %>%
  select(., -starts_with('blah')) %>%
  mutate(., treat = substr(treat, 6, nchar(treat)),
         estimate = ifelse(term %in% c('Th', 'Topt'), estimate - 273.15, estimate),
         treat = ifelse(treat == 'nophage', 'bact', 'bact_phage')) %>%
  group_by(., term, treat) %>%
  mean_qi()

# predictions from phage model
preds_phage <- data.frame(K = seq(15 + 273.15, 30 + 273.15, length.out = 200)) %>%
  add_fitted_draws(fit_bayes_phage, re_formula = NA) %>%
  data.frame() %>%
  group_by(K) %>%
  mean_qi(estimate = .value)

# predictions from bacteria model
preds_bact <-modelr::data_grid(d_bact, K = modelr::seq_range(K, n = 100), phage) %>%
  add_fitted_draws(fit_bayes, re_formula = NA) %>%
  data.frame() %>%
  group_by(K, phage) %>%
  mean_qi(estimate = .value)

# calculate the closest place where abund is closest to 5.12 (the inoculate)
CT_max <- filter(preds_phage, K > 273.15 + 26) %>%
  group_by() %>%
  do(., tibble(CTmax_low = .[which.min(abs(.$.lower - 0)),]$K - 273.15,
               CTmax = .[which.min(abs(.$estimate - 0)),]$K - 273.15,
               CTmax_high = .[which.min(abs(.$.upper - 0)),]$K - 273.15))

# calculate relative fitness between bacteria in the presence / absence of phage
preds_diff <- modelr::data_grid(d_bact, K = modelr::seq_range(K, n = 100), phage) %>%
  add_fitted_draws(fit_bayes, re_formula = NA) %>%
  data.frame() %>%
  select(., K, phage, .value, .draw) %>%
  spread(., phage, .value) %>%
  mutate(., .value = phage / nophage) %>%
  group_by(K) %>%
  mean_qi(estimate = .value)

#-----------------------------#
# plot TPCs and parameters ####
#-----------------------------#

# plot predictions
plot_phage <- ggplot() +
  geom_point(aes(K - 273.15, growth_rate), filter(d, K <= 303.15), size = 3, alpha = 0.5) +
  geom_line(aes(K - 273.15, estimate), data = preds_phage) +
  geom_ribbon(aes(K - 273.15, ymin = .lower, ymax = .upper), data = preds_phage, alpha = 0.2) +
  theme_bw(base_size = 14) +
  xlab('Temperature (ºC)') +
  ylab(expression('Phage replication rate ( hr'^-1~')')) +
  ggtitle('(a) Thermal performance of phage replication rate') +
  theme(legend.position = 'none') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  NULL

plot_bact <- ggplot() +
  geom_point(aes(K - 273.15, r, col = phage), d_bact, alpha = 0.5, size = 3) +
  geom_line(aes(K - 273.15, estimate, col = phage, group = interaction(phage)), data = preds_bact) +
  geom_ribbon(aes(K - 273.15, ymin = .lower, ymax = .upper, fill = phage), data = preds_bact, alpha = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c('blue', 'black')) +
  scale_color_manual(values = c('blue', 'black')) +
  xlab('Temperature (ºC)') +
  ylab(expression('Bacteria growth rate ( hr'^-1~')')) +
  ggtitle('(b) Thermal performance of bacterial growth rate') +
  theme(legend.position = 'none') +
  geom_vline(aes(xintercept = 29.2), linetype = 2, col = 'black') 

p1 <- plot_phage + plot_bact

ggsave('plots/Figure_1.png', p1, width = 12, height = 5)

# make Figure 2
plot_diff <- ggplot(preds_diff) +
  geom_line(aes(K - 273.15, estimate)) +
  geom_ribbon(aes(K - 273.15, ymin = .lower, ymax = .upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  theme_bw(base_size = 14) +
  xlab('Temperature (ºC)') +
  ylab('Relative fitness (phage presence / phage absence)') +
  ggtitle('(a) Change in bacterial growth rate in the presence of phage')

x_labs <- c('bacteria\nalone', 'bacteria\nwith phage')

plot_e <- filter(d_params_bact, term %in% c('E')) %>%
  ggplot(., aes(treat, estimate, col = treat)) + 
  geom_point(size = 3) + 
  geom_linerange(aes(ymin = .lower, ymax = .upper, col = treat)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  ylab('Activation energy (eV)') +
  scale_color_manual(values = c('blue', 'black')) +
  ggtitle('(b)') +
  xlab('') +
  scale_x_discrete(labels = x_labs)
  
plot_Topt <- filter(d_params_bact, term %in% c('Topt')) %>%
  ggplot(., aes(treat, estimate, col = treat)) + 
  geom_point(size = 3) + 
  geom_linerange(aes(ymin = .lower, ymax = .upper, col = treat)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  ylab('Optimum temperature (ºC)') +
  scale_color_manual(values = c('blue', 'black')) +
  ggtitle('(c)') +
  xlab('') +
  scale_x_discrete(labels = x_labs)

plot_rmax <- filter(d_params_bact, term %in% c('rmax')) %>%
  ggplot(., aes(treat, estimate, col = treat)) + 
  geom_point(size = 3) + 
  geom_linerange(aes(ymin = .lower, ymax = .upper, col = treat)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  ylab(expression('Optimal growth rate ('~hr^-1~')')) +
  scale_color_manual(values = c('blue', 'black')) +
  ggtitle('(d)') +
  xlab('') +
  scale_x_discrete(labels = x_labs)

plot_eh <- filter(d_params_bact, term %in% c('Eh')) %>%
  ggplot(., aes(treat, estimate, col = treat)) + 
  geom_point(size = 3) + 
  geom_linerange(aes(ymin = .lower, ymax = .upper, col = treat)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  ylab('Deactivation energy (eV)') +
  scale_color_manual(values = c('blue', 'black')) +
  ggtitle('(e)') +
  xlab('') +
  scale_x_discrete(labels = x_labs)
  
p2 <- plot_diff + {
  plot_e + plot_Topt + plot_rmax + plot_eh + plot_layout(ncol = 4)
} +
  plot_layout(nrow = 2, heights = c(0.7, 0.3))

ggsave('plots/Figure_2.png', p2, width = 9, height = 8)

#------------------------------------------------------#
# plot raw data for bacteria with and without phage ####
#------------------------------------------------------#

# extract parameters ####
d_params_bact <- fit_bayes_raw %>%
  spread_draws(`b_.*`, regex = TRUE) %>%
  mutate(., b_Topt_phagenophage = get_topt(b_Eh_phagenophage, b_Th_phagenophage, b_E_phagenophage),
         b_Topt_phagephage = get_topt(b_Eh_phagephage, b_Th_phagephage, b_E_phagephage),
         b_rmax_phagenophage = sharpeschoolhigh_1981(b_Topt_phagenophage, b_c_phagenophage, b_E_phagenophage, b_Eh_phagenophage, b_Th_phagenophage, 20),
         b_rmax_phagephage = sharpeschoolhigh_1981(b_Topt_phagephage, b_c_phagephage, b_E_phagephage, b_Eh_phagephage, b_Th_phagephage, 20),
         b_relrmax_rmax = 1 - b_rmax_phagephage/b_rmax_phagenophage) %>%
  gather(., 'param', 'estimate', 4:ncol(.)) %>%
  separate(., param, c('blah', 'term', 'treat'), sep = '_') %>%
  select(., -starts_with('blah')) %>%
  mutate(., treat = substr(treat, 6, nchar(treat)),
         estimate = ifelse(term %in% c('Th', 'Topt'), estimate - 273.15, estimate),
         treat = ifelse(treat == 'nophage', 'bact', 'bact_phage')) %>%
  group_by(., term, treat) %>%
  mean_qi()

# extract predictions ####
preds_bact <- modelr::data_grid(d_mumax, K = modelr::seq_range(K, n = 100), phage) %>%
  add_fitted_draws(fit_bayes_raw, re_formula = NA) %>%
  data.frame() %>%
  group_by(K, phage) %>%
  mean_qi(estimate = .value)

# get difference in predictions
preds_diff <- modelr::data_grid(d_mumax, K = modelr::seq_range(K, n = 100), phage) %>%
  add_fitted_draws(fit_bayes_raw, re_formula = NA) %>%
  data.frame() %>%
  select(., K, phage, .value, .draw) %>%
  spread(., phage, .value) %>%
  mutate(., .value = phage / nophage) %>%
  group_by(K) %>%
  mean_qi(estimate = .value)

# plot preds
plot_bact <- ggplot() +
  geom_point(aes(K - 273.15, r, col = phage), d_bact_raw, alpha = 0.5, size = 3) +
  geom_line(aes(K - 273.15, estimate, col = phage, group = interaction(phage)), data = preds_bact) +
  geom_ribbon(aes(K - 273.15, ymin = .lower, ymax = .upper, fill = phage), data = preds_bact, alpha = 0.2) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = c('blue', 'black')) +
  scale_color_manual(values = c('blue', 'black')) +
  xlab('Temperature (ºC)') +
  ylab(expression('Bacteria growth rate ( hr'^-1~')')) +
  ggtitle('(a) Thermal performance of bacterial growth rate') +
  theme(legend.position = 'none')

# plot diff
plot_diff <- ggplot(preds_diff) +
  geom_line(aes(K - 273.15, estimate)) +
  geom_ribbon(aes(K - 273.15, ymin = .lower, ymax = .upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  theme_bw(base_size = 12) +
  xlab('Temperature (ºC)') +
  ylab('Relative fitness (phage presence / phage absence)') +
  ggtitle('(b) Change in bacterial growth rate in the presence of phage')

# plot params

x_labs <- c('bacteria\nalone', 'bacteria\nwith phage')

plot_e <- filter(d_params_bact, term %in% c('E')) %>%
  ggplot(., aes(treat, estimate, col = treat)) + 
  geom_point(size = 3) + 
  geom_linerange(aes(ymin = .lower, ymax = .upper, col = treat)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  ylab('Activation energy (eV)') +
  scale_color_manual(values = c('blue', 'black')) +
  ggtitle('(c)') +
  xlab('') +
  scale_x_discrete(labels = x_labs)

plot_Topt <- filter(d_params_bact, term %in% c('Topt')) %>%
  ggplot(., aes(treat, estimate, col = treat)) + 
  geom_point(size = 3) + 
  geom_linerange(aes(ymin = .lower, ymax = .upper, col = treat)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  ylab('Optimum temperature (ºC)') +
  scale_color_manual(values = c('blue', 'black')) +
  ggtitle('(d)') +
  xlab('') +
  scale_x_discrete(labels = x_labs)

plot_rmax <- filter(d_params_bact, term %in% c('rmax')) %>%
  ggplot(., aes(treat, estimate, col = treat)) + 
  geom_point(size = 3) + 
  geom_linerange(aes(ymin = .lower, ymax = .upper, col = treat)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  ylab(expression('Optimal growth rate ('~hr^-1~')')) +
  scale_color_manual(values = c('blue', 'black')) +
  ggtitle('(e)') +
  xlab('') +
  scale_x_discrete(labels = x_labs)

plot_eh <- filter(d_params_bact, term %in% c('Eh')) %>%
  ggplot(., aes(treat, estimate, col = treat)) + 
  geom_point(size = 3) + 
  geom_linerange(aes(ymin = .lower, ymax = .upper, col = treat)) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  ylab('Deactivation energy (eV)') +
  scale_color_manual(values = c('blue', 'black')) +
  ggtitle('(f)') +
  xlab('') +
  scale_x_discrete(labels = x_labs)

p2 <- (plot_bact | plot_diff) /
  (plot_e | plot_Topt | plot_rmax | plot_eh) +
  plot_layout(heights = c(0.7, 0.3))


ggsave('plots/Figure_S9.png', p2, height = 8, width = 13)
