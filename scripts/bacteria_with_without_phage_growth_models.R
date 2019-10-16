#----------------------------------------------------------------------------------#
# growth curve analysis of bacteria and bacteria + phage logistic growth curves ####
#----------------------------------------------------------------------------------#

# clear workspace
rm(list = ls())

# load in packages 
library(nls.multstart)
library(tidyverse)
library(broom)
library(modelr)
library(nlsMicrobio)
library(dataViewer) # remotes::install_github('padpadpadpad/dataViewer')

# load in data
d <- readRDS('data/bacteria_with_without_phage.rds')

# load in growth curves models
source('scripts/growth_curve_models.R')

#-----------------#
# wrangle data ####
#-----------------#

# if od_cor is negative or less than the blank, replace with the minimum value possible for difference between control and sample (0.001)
d <- mutate(d, od_cor = od - control,
            od_cor2 = ifelse(od_cor <= 0, NA, od_cor)) %>%
  group_by(., temp, phage) %>%
  mutate(., od_cor2 = ifelse(is.na(od_cor2)|od_cor2 < 0.001, 0.001, od_cor2)) %>%
  ungroup()

# create log10 abundance
d <- mutate(d, log10_od_cor = log10(od_cor2),
            log10_od = log10(od),
            # parse timepoint number,
            time_point = readr::parse_number(time_point),
            time = ymd_hms(as.character(time)))

# add time in hours and other columns
d <- group_by(d, temp, phage, well) %>%
  arrange(., time_point) %>%
  mutate(., lag = lag(time, order_by = time_point), 
         t = as.numeric(time - lag(time, order_by = time_point), units = 'hours'),
         t = ifelse(is.na(t), 0, t),
         t = cumsum(t)) %>%
  ungroup() %>%
  select(-lag) %>%
  mutate(., id = group_indices(., temp, phage),
         log10_od = ifelse(is.nan(log10_od), NA, log10_od))

# remove curves that are just plain bad
d <- unite(d, 'id2', c(temp, well, phage), remove = FALSE) %>%
  filter(., id2 != '15_D_4_nophage')

# points to remove from the dataset

# can look at each curve using dataViewer https://github.com/padpadpadpad/dataViewer/ 
# anomalies <- dataViewer::dataViewer(d, 't', 'log10_od_cor', id_cols = Hmisc::Cs(temp, phage, nutrient, well))
anomalies <- readRDS('data/bacteria_with_without_phage_anomalies.rds')

# remove points from dataset
d2 <- dataViewer::get_unclicked(d, anomalies) 

#-----------------------------#
# run models on each curve ####
#-----------------------------#

# set upper and lower limits for parameter estimates and delete those that max out at those values

models = NULL

# run nls.multstart on all data ####
models <- group_by(d2, temp, phage, well, id) %>%
  nest() %>%
  dplyr::mutate(., 
                # first do a fit of all baranyi curves
                fits_baranyi = purrr::map(data, ~nls_multstart(log10_od_cor ~ baranyi(log10_nmax, log10_n0, mumax, t = t, lag),
                                                               data = .x,
                                                               iter = 500,
                                                               start_lower = c(log10_nmax = -2, log10_n0 = -4, mumax = 0, lag = 0),
                                                               start_upper = c(log10_nmax = 0, log10_n0 = -2, mumax = 5, lag = 15),
                                                               supp_errors = 'Y',
                                                               na.action = na.omit,
                                                               lower = c(log10_nmax = -5, log10_n0 = -10, mumax = 0, lag = 0),
                                                               upper = c(log10_nmax = 1, log10_n0 = -2.85, mumax = 10, lag = 40))),
                # then do a fit of baranyi curves without lag
                fits_baranyi_without_lag = purrr::map(data, ~nls_multstart(log10_od_cor ~ baranyi_without_lag(log10_nmax, log10_n0, mumax, t = t),
                                                                           data = .x,
                                                                           iter = 500,
                                                                           start_lower = c(log10_nmax = -2, log10_n0 = -4, mumax = 0),
                                                                           start_upper = c(log10_nmax = 0, log10_n0 = -2, mumax = 5),
                                                                           supp_errors = 'Y',
                                                                           na.action = na.omit,
                                                                           lower = c(log10_nmax = -5, log10_n0 = -10, mumax = 0),
                                                                           upper = c(log10_nmax = 1, log10_n0 = -2.85, mumax = 10))),
                # gompertz model fits
                fits_gompertz = purrr::map(data, ~nls_multstart(log10_od_cor ~ gompertz(log10_nmax, log10_n0, mumax, t = t, lag),
                                                                data = .x,
                                                                iter = 500,
                                                                start_lower = c(log10_nmax = -2, log10_n0 = -4, mumax = 0, lag = 0),
                                                                start_upper = c(log10_nmax = 0, log10_n0 = -2, mumax = 5, lag = 15),
                                                                supp_errors = 'Y',
                                                                na.action = na.omit,
                                                                lower = c(log10_nmax = -5, log10_n0 = -10, mumax = 0, lag = 0),
                                                                upper = c(log10_nmax = 1, log10_n0 = -2.85, mumax = 10, lag = 40))),
                # buchanan fits
                fits_buchanan = purrr::map(data, ~nls_multstart(log10_od_cor ~ buchanan(log10_nmax, log10_n0, mumax, t = t, lag),
                                                                data = .x,
                                                                iter = 500,
                                                                start_lower = c(log10_nmax = -2, log10_n0 = -4, mumax = 0, lag = 0),
                                                                start_upper = c(log10_nmax = 0, log10_n0 = -2, mumax = 5, lag = 15),
                                                                supp_errors = 'Y',
                                                                na.action = na.omit,
                                                                lower = c(log10_nmax = -5, log10_n0 = -10, mumax = 0, lag = 0),
                                                                upper = c(log10_nmax = 1, log10_n0 = -2.85, mumax = 10, lag = 40))),
                # buchanan without lag fits
                fits_buchanan_without_lag = purrr::map(data, ~nls_multstart(log10_od_cor ~ buchanan_without_lag(log10_nmax, log10_n0, mumax, t = t),
                                                                            data = .x,
                                                                            iter = 500,
                                                                            start_lower = c(log10_nmax = -2, log10_n0 = -4, mumax = 0),
                                                                            start_upper = c(log10_nmax = -0,log10_n0 = -2, mumax = 5),
                                                                            supp_errors = 'Y',
                                                                            na.action = na.omit,
                                                                            lower = c(log10_nmax = -5, log10_n0 = -10, mumax = 0),
                                                                            upper = c(log10_nmax = 1, log10_n0 = -2.85, mumax = 10))))

# stack models
model_stack <- gather(models, 'model', 'best_model', starts_with('fits')) %>%
  filter(., !is.null(best_model)) %>%
  mutate(., aic = map_dbl(best_model, possibly(MuMIn::AICc, otherwise = NA)))

# make model comparison plot

# add column for facets
model_stack <- mutate(model_stack, model_name = case_when(model == "fits_baranyi" ~ '(a) baranyi',
                                                          model == "fits_baranyi_without_lag" ~ '(b) buchanan without lag',
                                                          model == "fits_buchanan" ~ '(c) buchanan without lag',
                                                          model == "fits_buchanan_without_lag" ~ '(d) buchanan without lag',
                                                          model == "fits_gompertz" ~ '(e) gompertz'))

aic_means <- filter(model_stack, !is.na(aic)) %>%
  group_by(., model, model_name) %>%
  summarise(., mean_aic = mean(aic),
            median_aic = median(aic),
            count = n()) %>%
  ungroup()

# visualise aic of all curves
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
ggsave('plots/Figure_S10.png', last_plot(), width = 7, height = 10)

#--------------------#
# model selection ####
#--------------------#

# look at which model fits best overall
group_by(model_stack, temp, phage, well) %>%
  filter(aic == min(aic, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(., model) %>%
  tally() %>%
  mutate(., perc_best = round(n/sum(n)*100))
# gompertz is overwhelmingly the best

# refit with just the gompertz model
models <- group_by(d2, temp, phage, well, id) %>%
  nest() %>%
  dplyr::mutate(., # gompertz model fits
                fits_gompertz = purrr::map(data, ~nls_multstart(log10_od_cor ~ gompertz(log10_nmax, log10_n0, mumax, t = t, lag),
                                                                data = .x,
                                                                iter = 500,
                                                                start_lower = c(log10_nmax = -2, log10_n0 = -4, mumax = 0, lag = 0),
                                                                start_upper = c(log10_nmax = 0, log10_n0 = -2, mumax = 5, lag = 15),
                                                                supp_errors = 'Y',
                                                                na.action = na.omit,
                                                                lower = c(log10_nmax = -5, log10_n0 = -10, mumax = 0, lag = 0),
                                                                upper = c(log10_nmax = 1, log10_n0 = -2.85, mumax = 10, lag = 40))))

model_stack <- gather(models, 'model', 'best_model', starts_with('fits')) %>%
  filter(., !is.null(best_model)) %>%
  mutate(., aic = map_dbl(best_model, possibly(MuMIn::AICc, otherwise = NA)))

#-----------------------#
# look for poor fits ####
#-----------------------#

# do not want any parameter estimates at the limits of parameter space
start_lims <- tibble(
  log10_nmax = c(-5, 1),
  log10_n0 = c(-10, -2.75),
  mumax = c(0, 10),
  lag = c(0, 40)
)

params <- model_stack %>%
  mutate(., params = map(best_model, tidy)) %>%
  unnest(params)  %>%
  select(., -data, -model, -best_model, -model)

# drop anything with log10_nmax <= -3 - represents an OD of 0.001
# this is the smallest possible number of OD - control - cannot be right

# drop when log10_n0 > log10_nmax
curves_id_n0_higher_nmax <- select(params, temp, phage, well, term, estimate, id) %>%
  spread(., term, estimate) %>%
  filter(log10_n0 > log10_nmax)

# curves to drop based on start parameters
curve_id_drop <- select(params, temp, phage, well, term, estimate, id) %>%
  spread(., term, estimate) %>%
  filter(., lag == 40 |
           log10_nmax %in% start_lims$log10_nmax |
           log10_nmax <= -3 |
           log10_n0 %in% start_lims$log10_n0 |
           mumax %in% start_lims$mumax) %>%
  pull(., id)

# no curves to drop
curve_id_drop
curves_id_n0_higher_nmax

#-----------------------------------------#
# calculate predictions and parameters ####
#-----------------------------------------#

# get parameters of best model and extract mumax
params <- model_stack %>%
  mutate(., params = map(best_model, tidy)) %>%
  unnest(params)  %>%
  select(., -data, -model, -best_model, -model)

# get preds of gompertz model
new_preds <- d2 %>%
  do(., data.frame(t = seq(min(.$t), max(.$t), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(d2, temp, phage, well) %>%
  summarise(., min_t = min(t), max_t = max(t)) %>%
  ungroup()

# smoothed preds
preds <- model_stack %>%
  mutate(., preds = map(best_model, augment, newdata = new_preds)) %>%
  unnest(preds)  %>%
  select(., -data, -model, -best_model, -model) %>%
  merge(., max_min, by = c('temp', 'phage', 'well')) %>%
  group_by(., temp, phage, well) %>%
  filter(., t > unique(min_t) & t < unique(max_t)) %>%
  ungroup() %>%
  mutate(., id = group_indices(., temp, phage))

# plot logistical growth curves ####
mark_facets_degree <- function(string){
  len <- length(string)
  string = paste('(', letters[1:len], ') ', string, ' ºC', sep = '')
  return(string)
}

plot_bact_phage <- ggplot(filter(d), aes(t, log10_od_cor, col = phage)) +
  geom_point(position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  geom_line(aes(t, .fitted, group = interaction(well, phage)), data = filter(preds), alpha = 0.5) +
  facet_wrap(~ temp, labeller = labeller(temp = mark_facets_degree), ncol = 2) +
  theme_bw(base_size = 16) +
  xlab('Time (hours)') +
  ylab(expression(log[10]~OD[600])) +
  scale_color_manual(values = c('blue', 'black')) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0))

#--------------------------------#
# look at lag phase extension ####
#--------------------------------#

d_lag <- filter(params, term == 'lag')

mod_lag <- lm(estimate ~ temp*phage, d_lag)
mod_lag2 <- lm(estimate ~ temp+phage, d_lag)
anova(mod_lag, mod_lag2)

d_lag <- bind_cols(d_lag, predict(mod_lag, interval = 'confidence') %>% data.frame())

ggplot(d_lag, aes(temp, estimate)) +
  geom_point(aes(col = phage),size = 3) +
  geom_line(aes(temp, fit, col = phage)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = phage), alpha = 0.25) +
  theme_bw(base_size = 14) +
  theme(legend.position = 'none') +
  scale_fill_manual(values = c('blue', 'black')) +
  scale_color_manual(values = c('blue', 'black')) +
  labs(title = 'Lag time in the presence and absence of phage',
       subtitle = 'As defined by the Gompertz model',
       y = 'Lag time (hours)',
       x = 'Temperature (ºC)')

ggsave('plots/change_in_lag_time.png', last_plot(), width = 7, height = 5)

#------------------------------------------------------------#
# calculate time to stationary phase for each phage curve ####
#------------------------------------------------------------#

# calculate at 90% of log10 nmax (as it is negative)

d_ts <- filter(preds, phage == 'phage') %>%
  select(id, temp, phage, well, t, .fitted) %>%
  merge(., filter(params, phage == 'phage' & term == 'log10_nmax') %>% select(id, temp, phage, well, estimate), by = c('phage', 'well', 'id', 'temp')) %>%
  group_by(temp, id, phage, well) %>%
  filter(., .fitted > (estimate+(estimate*0.1))) %>%
  filter(., t == min(t))

ggplot(d_ts, aes(temp, t)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_hline(aes(yintercept = 24), linetype = 2) +
  geom_hline(aes(yintercept = 48), linetype = 2) +
  geom_hline(aes(yintercept = 12), linetype = 2) +
  theme_bw(base_size = 16) +
  xlab('Temperature (ºC)') +
  ylab('Time until stationary phase (hours)') +
  ylim(c(0, 70))

ggsave('plots/Figure_S15.png', last_plot(), width = 7, height = 5)

#------------------------------------------#
# save out exponential growth rate data ####
#------------------------------------------#

d_mumax <- filter(params, term == 'mumax') 

# wrangle data
d_mumax <- mutate(d_mumax, K = temp + 273.15) %>%
  ungroup()

# save out dataframe
select(d_mumax, temp, phage, estimate) %>%
  saveRDS('data/bacteria_with_without_phage_growth_rates.rds')

#----------------------------------------------------------------#
# check importance of removed points in growth rate estimates ####
#----------------------------------------------------------------#

# fit models with the raw data
models_raw <- mutate(d, source = 'raw') %>%
  group_by(., temp, phage, well, id, source) %>%
  nest() %>%
  dplyr::mutate(., # gompertz model fits
                fits_gompertz = purrr::map(data, ~nls_multstart(log10_od_cor ~ gompertz(log10_nmax, log10_n0, mumax, t = t, lag),
                                                                data = .x,
                                                                iter = 500,
                                                                start_lower = c(log10_nmax = -2, log10_n0 = -4, mumax = 0, lag = 0),
                                                                start_upper = c(log10_nmax = 0, log10_n0 = -2, mumax = 5, lag = 15),
                                                                supp_errors = 'Y',
                                                                na.action = na.omit,
                                                                lower = c(log10_nmax = -5, log10_n0 = -10, mumax = 0, lag = 0),
                                                                upper = c(log10_nmax = 1, log10_n0 = -2.85, mumax = 10, lag = 40))))

models <- mutate(models, source = 'cleaned')

# bind these together
models <- bind_rows(models, models_raw)

# get preds of gompertz model
new_preds <- d2 %>%
  do(., data.frame(t = seq(min(.$t), max(.$t), length.out = 150), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(d2, temp, phage, well) %>%
  summarise(., min_t = min(t), max_t = max(t)) %>%
  ungroup()

# smoothed preds
preds <- models %>%
  mutate(., preds = map(fits_gompertz, augment, newdata = new_preds)) %>%
  unnest(preds)  %>%
  select(., -data, -fits_gompertz) %>%
  merge(., max_min, by = c('temp', 'phage', 'well')) %>%
  group_by(., temp, phage, well) %>%
  filter(., t > unique(min_t) & t < unique(max_t)) %>%
  ungroup() %>%
  mutate(., id = group_indices(., temp, phage))

# plot anomalies vs non-anomalies
temps <- unique(d$temp)

mark_facets_custom <- function(string){
  len <- length(string)
  string = paste('(', letters[1:len], ') ', gsub('_.*', '', string), sep = '')
  return(string)
}

# add column for labelling
d2 <- unite(d2, id_plot, c(temp, phage, well), sep = '_', remove = FALSE)
preds <- unite(preds, id_plot, c(temp, phage, well), sep = '_', remove = FALSE)
anomalies <- unite(anomalies, id_plot, c(temp, phage, well), sep = '_', remove = FALSE)
d <- unite(d, id_plot, c(temp, phage, well), sep = '_', remove = FALSE)

d2 <- mutate(d2, label = paste(temp, 'ºC', phage, '_', well, sep = ' '))
anomalies <- mutate(anomalies, label = paste(temp, 'ºC', phage, '_', well, sep = ' '))
preds <- mutate(preds, label = paste(temp, 'ºC', phage,'_', well, sep = ' '))
d <- mutate(d, label = paste(temp, 'ºC', phage, '_',well, sep = ' '))

# loop through plots
for(i in 1:length(temps)){
  p_temp <- ggplot(filter(d2, temp == temps[i]), aes(t, log10_od_cor, col = phage)) +
    geom_line(data = filter(d, temp == temps[i]), alpha = 0.5) +
    geom_point(position = position_jitter(width = 0.1, height = 0), alpha = 1) +
    geom_point(col = 'red', data = filter(anomalies, temp == temps[i]), size = 3) +
    geom_line(aes(t, .fitted), col = 'red', data = filter(preds, temp == temps[i] & source == 'raw'), alpha = 1) +
    geom_line(aes(t, .fitted), data = filter(preds, temp == temps[i] & source == 'cleaned'), alpha = 1) +
    facet_wrap(~ label, labeller = labeller(label = mark_facets_custom)) +
    theme_bw(base_size = 16) +
    xlab('Time (hours)') +
    ylab(expression(log[10]~OD[600])) +
    scale_color_manual(values = c('blue', 'black')) +
    theme(legend.position = 'none',
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0))
  
  ggsave(paste('plots/Figure_S', i, '.png', sep = ''), p_temp, width = 10, height = 6)
}

# saveout raw mumax values
d_mumax <- models_raw %>%
  mutate(., params = map(fits_gompertz, tidy)) %>%
  unnest(params)  %>%
  select(., -data, -fits_gompertz) %>%
  filter(term == 'mumax')

# wrangle data
d_mumax <- mutate(d_mumax, K = temp + 273.15) %>%
  ungroup()

# save out dataframe
select(d_mumax, temp, phage, estimate) %>%
  saveRDS('data/bacteria_with_without_phage_growth_rates_raw.rds')
