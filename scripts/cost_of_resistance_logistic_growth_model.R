#--------------------------------------------------------------------------------#
# analysis of resistant and susceptible clones using logistical growth curves ####
#--------------------------------------------------------------------------------#

# clear workspace
rm(list = ls())

# load in packages
library(nls.multstart)
library(broom)
library(modelr)
library(nlsMicrobio)
library(patchwork)
library(tidyverse)

# load in data
d <- readRDS('data/cost_of_resistance.rds')

# load in growth curves models
source('scripts/growth_curve_models.R')

#-----------------#
# wrangle data ####
#-----------------#

# if od_cor is negative or less than the blank, replace with the minimum value possible for difference between control and sample (0.001)
d <- mutate(d, od_cor = od - control,
            od_cor2 = ifelse(od_cor <= 0.001, NA, od_cor))

# create log10 abundance
d <- mutate(d, log10_od_cor = log10(od_cor2),
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
  # remove NAs and infinites
  mutate(., id = group_indices(., temp, phage),
         log10_od_cor = ifelse(is.nan(log10_od_cor), NA, log10_od_cor))

# remove NAs and infinites
d2 <- filter(d, !is.na(log10_od_cor) & !is.infinite(log10_od_cor))

# data cleaning
# remove the t0 time point (prone to bubbles)
# here we set t0 per curve to the first value where there was a measurement
d2 <- filter(d2, t >= 1) %>%
  group_by(., temp, phage, clone, rep) %>%
  mutate(., t2 = t - min(t)) %>%
  ungroup()

#-----------------------------#
# run models on each curve ####
#-----------------------------#

# set upper and lower limits for parameter estimates and delete those that max out at those values

models = NULL

if(is.null(models)){
  # run nls.multstart on all data
  models <- group_by(d2, temp, phage, rep, clone) %>%
    nest() %>%
    dplyr::mutate(., 
                  # first do a fit of all baranyi curves
                  fits_baranyi = purrr::map(data, ~nls_multstart(log10_od_cor ~ baranyi(log10_nmax, log10_n0, mumax, t = t2, lag),
                                                                 data = .x,
                                                                 iter = 500,
                                                                 start_lower = c(log10_nmax = -2, log10_n0 = -4, mumax = 0, lag = 0),
                                                                 start_upper = c(log10_nmax = 0, log10_n0 = -2, mumax = 5, lag = 15),
                                                                 supp_errors = 'Y',
                                                                 na.action = na.omit,
                                                                 lower = c(log10_nmax = -5, log10_n0 = -10, mumax = 0, lag = 0),
                                                                 upper = c(log10_nmax = 1, log10_n0 = -1.5, mumax = 10, lag = 40))),
                  # then do a fit of baranyi curves without lag
                  fits_baranyi_without_lag = purrr::map(data, ~nls_multstart(log10_od_cor ~ baranyi_without_lag(log10_nmax, log10_n0, mumax, t = t2),
                                                                             data = .x,
                                                                             iter = 500,
                                                                             start_lower = c(log10_nmax = -2, log10_n0 = -4, mumax = 0),
                                                                             start_upper = c(log10_nmax = 0, log10_n0 = -2, mumax = 5),
                                                                             supp_errors = 'Y',
                                                                             na.action = na.omit,
                                                                             lower = c(log10_nmax = -5, log10_n0 = -10, mumax = 0),
                                                                             upper = c(log10_nmax = 1, log10_n0 = -1.5, mumax = 10))),
                  # gompertz model fits
                  fits_gompertz = purrr::map(data, ~nls_multstart(log10_od_cor ~ gompertz(log10_nmax, log10_n0, mumax, t = t2, lag),
                                                                  data = .x,
                                                                  iter = 500,
                                                                  start_lower = c(log10_nmax = -2, log10_n0 = -4, mumax = 0, lag = 0),
                                                                  start_upper = c(log10_nmax = 0, log10_n0 = -2, mumax = 5, lag = 15),
                                                                  supp_errors = 'Y',
                                                                  na.action = na.omit,
                                                                  lower = c(log10_nmax = -5, log10_n0 = -10, mumax = 0, lag = 0),
                                                                  upper = c(log10_nmax = 1, log10_n0 = -1.5, mumax = 10, lag = 40))),
                  # buchanan fits
                  fits_buchanan = purrr::map(data, ~nls_multstart(log10_od_cor ~ buchanan(log10_nmax, log10_n0, mumax, t = t2, lag),
                                                                  data = .x,
                                                                  iter = 500,
                                                                  start_lower = c(log10_nmax = -2, log10_n0 = -4, mumax = 0, lag = 0),
                                                                  start_upper = c(log10_nmax = 0, log10_n0 = -2, mumax = 5, lag = 15),
                                                                  supp_errors = 'Y',
                                                                  na.action = na.omit,
                                                                  lower = c(log10_nmax = -5, log10_n0 = -10, mumax = 0, lag = 0),
                                                                  upper = c(log10_nmax = 1, log10_n0 = -1.5, mumax = 10, lag = 40))),
                  # buchanan without lag fits
                  fits_buchanan_without_lag = purrr::map(data, ~nls_multstart(log10_od_cor ~ buchanan_without_lag(log10_nmax, log10_n0, mumax, t = t2),
                                                                              data = .x,
                                                                              iter = 500,
                                                                              start_lower = c(log10_nmax = -2, log10_n0 = -4, mumax = 0),
                                                                              start_upper = c(log10_nmax = -0,log10_n0 = -2, mumax = 5),
                                                                              supp_errors = 'Y',
                                                                              na.action = na.omit,
                                                                              lower = c(log10_nmax = -2.5, log10_n0 = -10, mumax = 0),
                                                                              upper = c(log10_nmax = 1, log10_n0 = -1.5, mumax = 10))))
}


# this takes a long time. Go and get a cup of tea.
# save model output
saveRDS(models, 'data/models/cost_of_resistance_models.rds')

# stack models and calculate AIC score for each model
model_stack <- gather(models, 'model', 'best_model', starts_with('fits')) %>%
  filter(., !is.null(best_model)) %>%
  mutate(., aic = map_dbl(best_model, possibly(MuMIn::AICc, otherwise = NA)))

# add column for facets
model_stack <- mutate(model_stack, model_name = case_when(model == "fits_baranyi" ~ '(a) baranyi',
                                                          model == "fits_baranyi_without_lag" ~ '(b) baranyi without lag',
                                                          model == "fits_buchanan" ~ '(c) buchanan without lag',
                                                          model == "fits_buchanan_without_lag" ~ '(d) buchanan without lag',
                                                          model == "fits_gompertz" ~ '(e) gompertz'))

aic_means <- filter(model_stack, !is.na(aic)) %>%
  group_by(., model_name) %>%
  summarise(., mean_aic = mean(aic),
            median_aic = median(aic),
            count = n()) %>%
  ungroup()

# visualise each set of fits
ggplot(model_stack, aes(aic)) +
  geom_histogram(fill = 'lightgrey', col = 'black') +
  geom_vline(aes(xintercept = mean_aic), col = 'red', aic_means) +
  geom_vline(aes(xintercept = median_aic), col = 'blue', aic_means) +
  facet_wrap(~ model_name, ncol = 2) +
  theme_bw(base_size = 16) +
  ylab('Count') +
  xlab('AICc score') +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0))

# save plot
ggsave('plots/Figure_S11.png', last_plot(), width = 7, height = 10)

#---------------------------------------------------------------#
# model selection and calculating predictions and parameters ####
#---------------------------------------------------------------#

# estimate best overall model
# calculated as sum(lowest AIC for each curve)
best_model <- group_by(model_stack, temp, phage, rep, clone) %>%
  filter(aic == min(aic, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(., model) %>%
  tally() %>%
  mutate(., perc_best = round(n/sum(n)*100))
# the baranyi-without-lag model was the best 63% of the time

# filter for the best model
model_best <- filter(model_stack, model == 'fits_baranyi_without_lag')

# get parameters of best model and extract mumax
params <- model_best %>%
  mutate(., params = map(best_model, tidy)) %>%
  unnest(params)  %>%
  select(., -data, -model, -best_model, -model_name)

# get preds of best mod
new_preds <- d2 %>%
  do(., data.frame(t2 = seq(min(.$t2), max(.$t2), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(d2, temp, phage, rep, clone) %>%
  summarise(., min_t = min(t2), max_t = max(t2)) %>%
  ungroup()

# smoothed preds
preds <- model_best %>%
  mutate(., preds = map(best_model, augment, newdata = new_preds)) %>%
  unnest(preds)  %>%
  select(., -data, -model, -best_model, -model_name) %>%
  merge(., max_min, by = c('temp', 'phage', 'rep', 'clone')) %>%
  group_by(., temp, phage, rep, clone) %>%
  filter(., t2 > unique(min_t) & t2 < unique(max_t)) %>%
  ungroup() %>%
  mutate(., id = group_indices(., temp, phage))

# preds without new data
preds2 <- model_best %>%
  mutate(., preds = map(best_model, augment)) %>%
  unnest(preds)  %>%
  select(., -data, -model, -best_model, -model_name) %>%
  merge(., max_min, by = c('temp', 'phage', 'rep', 'clone')) %>%
  group_by(., temp, phage, rep, clone) %>%
  filter(., t2 > unique(min_t) & t2 < unique(max_t)) %>%
  ungroup() %>%
  mutate(., id = group_indices(., temp, phage))

#-------------------------------------------------------#
# analysis and plotting of predictions and residuals ####
#-------------------------------------------------------#

mark_facets_degree <- function(string){
  len <- length(string)
  string <- stringr::str_extract(string, '[[:digit:]]+')
  string = paste('(', letters[1:len], ') ', string, ' ºC', sep = '')
  return(string)
}

# plot predictions
p_preds <- ggplot(d2, aes(t2, log10_od_cor, col = phage)) +
  geom_point(position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  geom_line(aes(t2, .fitted, group = interaction(clone, rep, phage)), data = filter(preds), alpha = 0.5) +
  facet_wrap(~ interaction(phage, temp), labeller = labeller(`interaction(phage, temp)` = mark_facets_degree, .multi_line = FALSE), ncol = 2) +
  theme_bw(base_size = 16) +
  xlab('Time (hours)') +
  ylab(expression(log[10]~OD[600])) +
  scale_color_manual(values = c('blue', 'black')) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0))

ggsave('plots/Figure_S12.png', p_preds, width = 7, height = 10)

# plot residuals
p_resid <- ggplot(preds2, aes(t2, .resid, col = phage)) +
  geom_point(alpha = 0.5) +
  geom_vline(aes(xintercept = 7), linetype = 2) +
  facet_wrap(~ interaction(phage, temp), labeller = labeller(`interaction(phage, temp)` = mark_facets_degree, .multi_line = FALSE), ncol = 2) +
  geom_hline(aes(yintercept = 0)) +
  theme_bw(base_size = 14) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  scale_color_manual(values = c('blue', 'black')) +
  xlab('Time (hours)') +
  ylab('Residuals')

ggsave('plots/Figure_S13.png', p_resid, width = 7, height = 10)

# check slope of regressions 
slope <- filter(preds2, t2 < 7) %>%
  group_by(., temp, rep, phage, clone) %>%
  nest() %>%
  mutate(., model = purrr::map(data, ~lm(.resid ~ t2, data = .x))) %>%
  mutate(params = map(model, tidy)) %>%
  unnest(params) %>%
  select(., -data, -model) %>%
  filter(term == 't2') %>%
  mutate(temp2 = paste0('T_', temp))

# plot these slopes
ggplot(slope, aes(temp2, estimate, col = phage, fill = phage)) +
  geom_hline(aes(yintercept = 0)) +
  MicrobioUoE::geom_pretty_boxplot(width = 0.6) +
  geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.1), shape = 21, fill = 'white') +
  scale_color_manual(values = c('blue', 'black')) +
  scale_fill_manual(values = c('blue', 'black')) +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  xlab('Temperature (ªC)') +
  ylab(expression(slope~of~residuals~through~time~(hr^-1))) +
  scale_x_discrete(labels = unique(d2$temp))

ggsave('plots/Figure_S14.png', last_plot(), width = 7, height = 5)
