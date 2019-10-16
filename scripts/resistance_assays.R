# analysis of phage streaks through time

# load packages
library(tidyverse)

#-----------------------------#
# read in data and wrangle ####
#-----------------------------#

d <- read.csv('data/phage_streaks.csv', stringsAsFactors = FALSE) %>%
  filter(., temp != 'inoc') %>%
  mutate(., infect = ifelse(infect == '-', NA, as.numeric(infect)))

# create a data frame for the inoculate
d_inoc <- data.frame(expand.grid(temp = unique(d$temp), rep = unique(d$rep), clone = 1:12, stringsAsFactors = FALSE)) %>%
  mutate(., time_point = 0,
         infect = 1) 

# add hours column
d <- bind_rows(d, d_inoc) %>%
  mutate(., temp = as.numeric(temp),
         hours = case_when(time_point == 0 ~ 0,
                           time_point == 1 ~ 12,
                           time_point == 2 ~ 24,
                           time_point == 3 ~ 48)) %>%
  filter(., !is.na(infect))

# delete 37 ºC because they do not maintain viable populations
d <- filter(d, temp < 37)

# create data ready for analysis
d_sum <- group_by(d, temp, hours, rep) %>%
  # calculate total number of clones that worked
  summarise(., n = n(),
            # calculate number of resistant clones
            infect = sum(infect)) %>%
  # calculate number of susceptible clines
  mutate(., not_infect = n - infect) %>%
  ungroup() %>%
  mutate(., obs = 1:n(),
         # make temperature and time factors
         temp_fac = as.character(temp),
         hours_fac = as.character(hours),
         # add 1 to all resistant numbers
         infect2 = 1 + infect,
         # add 1 to all susceptible numbers
         not_infect2 = 1 + not_infect,
         # calculate total number of clones used in the analysis for each replicate
         n2 = infect2 + not_infect2)

# get rid of time 0
d_sum_filt <- filter(d_sum, hours > 0)

#----------------#
# fit a model ####
#----------------#

# time and temp as factors
# weighted by the number of observations
# quasibinomial error structure
mod1 <- glm(cbind(not_infect2, infect2) ~ temp_fac * hours_fac, d_sum_filt, family = 'quasibinomial', weights = n2)
mod2 <- glm(cbind(not_infect2, infect2) ~ temp_fac + hours_fac, d_sum_filt, family = 'quasibinomial', weights = n2)

anova(mod1, test = 'Chisq')
anova(mod1, mod2, test ='F')

# check for overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(mod1) # ratio would indicate overdispersion
overdisp_fun(mod2)

# plot diagnostic plots
plot(mod1)
# not great, but seem ok

# do multiple comparisons on resistance through time at each temperature - Table S3
emmeans::emmeans(mod1, pairwise ~ hours_fac|temp_fac, type = 'response')

#-----------------------------#
# get predictions and plot ####
#-----------------------------#

# get predictions from model as a factor
preds_fac <- select(d_sum_filt, hours_fac, temp_fac) %>%
  distinct() %>%
  mutate(., pred_mean = predict(mod1, newdata = ., type = 'link'),
         pred_se = predict(mod1, newdata = ., se = TRUE)$se,
         CI_low = pred_mean - (1.96*pred_se),
         CI_high = pred_mean + (1.96*pred_se)) %>%
  mutate_at(., c('pred_mean', 'CI_low', 'CI_high'), boot::inv.logit) %>%
  mutate(., hours = as.numeric(hours_fac),
         temp = as.numeric(temp_fac))

# data for phage TPC
phage <- tibble(xmin = c(26.5, 29.1),
                xmax = c(27.5, 29.4),
                ymin = c(0, 0),
                ymax = c(1,1),
                group = c(1,2))

# plot - makes figure 3
ggplot() +
  geom_point(aes(temp, 1 - infect/n), d_sum_filt, position = position_jitter(height = 0.01, width = 0.1), alpha = 0.1) +
  geom_point(aes(temp, pred_mean), preds_fac, size = 3) +
  # geom_line(aes(temp, 1 - pred_mean), preds_fac) +
  geom_linerange(aes(temp, ymin = CI_low, ymax = CI_high), preds_fac) +
  facet_wrap(~ hours, ncol = 3, labeller = labeller(hours = c('12' = '(a) 12 hours', '24' = '(b) 24 hours', '48' = '(c) 48 hours'))) +
  xlab('Temperature (º C)') +
  ylab('Resistance') +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  theme(strip.text.x = element_text(hjust = 0, size = 14),
        strip.background.x = element_blank()) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = group), fill = 'black', phage, alpha = 0.3)

ggsave('plots/Figure_3.png', last_plot(), width = 10, height = 4)



