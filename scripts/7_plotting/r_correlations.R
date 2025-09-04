pacman::p_load(brms, tidyverse, bayesplot, ggridges)
source("scripts/theme_ggplot.R")

## load data
load(file = "output/loads.RData")

chick_gerp <- subset(loads, loadtype == "gerp" & age == "chick")
chick_high <- subset(loads, loadtype == "high" & age == "chick")
ad_gerp <- subset(loads, loadtype == "gerp" & age == "adult")
ad_high <- subset(loads, loadtype == "high" & age == "adult")

gerp <- subset(loads, loadtype=="gerp")
high <- subset(loads, loadtype=="high")

load_wide <- left_join(gerp, high, by = "id")
load_wide$age.x <- gsub("adult", "Adults and yearlings", load_wide$age.x)
load_wide$age.x <- gsub("chick", "Chicks", load_wide$age.x)
ggplot(load_wide, aes(x = total_load.x, y = total_load.y, fill = age.x), col = "black") + 
  geom_point(shape=21, size = 3) + 
  geom_smooth(method="lm", aes(col = age.x), linewidth = 2) + 
  labs(x = "Total GERP load", y = "Total SnpEff load", col = "Age category", fill = "Age category") + 
  scale_fill_manual(values=alpha(c(clrs_hunting[4], clrs_hunting[2]), 0.7))+
  scale_color_manual(values=alpha(c(clrs_hunting[4], clrs_hunting[2]), 0.7)) +
  theme(legend.position = "bottom")-> plot_cor_loads

ggsave(plot_cor_loads, file = "plots/compare_loads_gerp_high.png", width = 10, height = 8)


cor.test(load_wide$total_load.x[which(load_wide$age.x=="Chicks")], load_wide$total_load.y[which(load_wide$age.x=="Chicks")])


### loads over time ###

loads$age <- factor(loads$age, levels = c("adult", "chick"))
load("data/phenotypic/phenotypes_lifetime.RData")
loads <- left_join(loads, pheno_wide[,c("id", "lifespan")], by = "id")
loads <- loads %>% mutate(lifespan_cat = as.factor(case_when(
  lifespan == 1 ~ "yearling",
  lifespan > 1 ~ "adult"
)))
loads$lifespan_cat <- factor(loads$lifespan_cat, levels = c("adult", "yearling"))

ggplot(loads, aes(x = year, y = total_load, fill = age)) + geom_point() + geom_line() + 
  facet_grid(~loadtype)

loads_sum <- loads %>% group_by(age, loadtype, year) %>% summarise(
  n = n(),
  mean = mean(total_load, na.rm=TRUE),
  sd = sd(total_load))

loads_sum$age <- gsub("adult", "Adults and yearlings", loads_sum$age)
loads_sum$age <- gsub("chick", "Chicks", loads_sum$age)
loads_sum$loadtype <- gsub("gerp", "GERP", loads_sum$loadtype)
loads_sum$loadtype <- gsub("high", "SnpEff", loads_sum$loadtype)

loads_sum$age <- factor(loads_sum$age, levels = c("Chicks", "Adults and yearlings"))
ggplot(subset(loads_sum, loadtype=="GERP"), aes(x = year, y = mean, col = age, fill = age)) + 
  geom_point(shape=21, aes(size = n), position=position_dodge(.9)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
 # facet_wrap(~loadtype, ncol=1, scale="free") +
  theme(legend.position="bottom")+
  scale_fill_manual(values=alpha(c(clrs_hunting[2], clrs_hunting[4]), 0.7))+
  scale_color_manual(values=alpha(c(clrs_hunting[2], clrs_hunting[4]), 0.7)) +
  labs(x = "Birth year", y = "Total load", fill = "Age class", col = "Age class", size = "Sample size") +
  theme( plot.margin = margin(2,1,1,1, "cm"),
         panel.spacing = unit(3,"lines"))+
  guides(fill  = guide_legend(order = 1),
         col = guide_legend(order=1),
         size = guide_legend(order = 2))-> time_size_gerp

time_size_gerp

ggplot(subset(loads_sum, loadtype=="SnpEff"), aes(x = year, y = mean, col = age, fill = age)) + 
  geom_point(shape=21, aes(size = n), position=position_dodge(.9)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  # facet_wrap(~loadtype, ncol=1, scale="free") +
  theme(legend.position="bottom")+
  scale_fill_manual(values=alpha(c(clrs_hunting[2], clrs_hunting[4]), 0.7))+
  scale_color_manual(values=alpha(c(clrs_hunting[2], clrs_hunting[4]), 0.7)) +
  labs(x = "Birth year", y = "Total load", fill = "Age class", col = "Age class", size = "Sample size") +
  theme( panel.spacing = unit(3,"lines"),
         plot.margin = margin(2,1,1,1, "cm"))+
  guides(fill  = guide_legend(order = 1),
         col = guide_legend(order=1),
         size = guide_legend(order = 2))-> time_size_high

time_size_high

plot_grid(time_size_gerp, time_size_high,
          ncol = 1, 
          labels = c("A) GERP", "B) SnpEff"), label_fontface = "plain", label_size = 22) -> sup_time
sup_time

ggsave(sup_time, file = "plots/plot_loads_time.png", width = 12, height = 10)

# model
summary(lm(total_load ~ year, data = subset(loads, loadtype=="gerp"& age == "chick")))
summary(lm(total_load ~ year, data = subset(loads, loadtype=="gerp"& age == "adult")))
summary(lm(total_load ~ year, data = subset(loads, loadtype=="gerp"& lifespan_cat == "yearling")))
summary(lm(total_load ~ year, data = subset(loads, loadtype=="gerp"& lifespan_cat == "adult")))

summary(lm(total_load ~ site, data = subset(loads, loadtype=="gerp"& age == "chick")))
summary(lm(total_load ~ site, data = subset(loads, loadtype=="gerp"& age == "adult")))
summary(lm(total_load ~ site, data = subset(loads, loadtype=="gerp"& lifespan_cat == "yearling")))
summary(lm(total_load ~ site, data = subset(loads, loadtype=="gerp"& lifespan_cat == "adult")))

#bayesian
## parameters
iter = 1000000
burn = 500000
thin = 1000

## total
brm_gerp_year <- brm(total_load ~ year, data = subset(loads, loadtype == "gerp"),
                family = "gaussian",
                prior = prior(normal(0,1), class = b),
                cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                iter = iter, thin = thin, warmup = burn, seed = 1908)
save(brm_gerp_year, file = "output/5_models/brm_gerp_year.RData")

mcmc_intervals_data(brm_gerp_year, prob =0.8, prob_outer = 0.95)

brm_snpeff_year <- brm(total_load ~ year, data = subset(loads, loadtype == "high"),
                     family = "gaussian",
                     prior = prior(normal(0,1), class = b),
                     cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                     iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_snpeff_year, file = "output/5_models/brm_high_year.RData")

mcmc_intervals_data(brm_snpeff_year, prob =0.8, prob_outer = 0.95)

brm_gerp_site <- brm(total_load ~ site, data = subset(loads, loadtype == "gerp"),
                     family = "gaussian",
                     prior = prior(normal(0,1), class = b),
                     cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                     iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_gerp_site, file = "output/5_models/brm_gerp_site.RData")

mcmc_intervals_data(brm_gerp_site, prob =0.8, prob_outer = 0.95)

brm_high_site <- brm(total_load ~ site, data = subset(loads, loadtype == "high"),
                     family = "gaussian",
                     prior = prior(normal(0,1), class = b),
                     cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                     iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_high_site, file = "output/5_models/brm_high_site.RData")

mcmc_intervals_data(brm_high_site, prob =0.8, prob_outer = 0.95)

### differences among leks  

loads_sum_lek <- loads %>% group_by(age, loadtype, site) %>% summarise(
  n = n(),
  mean = mean(total_load, na.rm=TRUE),
  sd = sd(total_load))

loads_sum_lek$age <- gsub("adult", "Adults and yearlings", loads_sum_lek$age)
loads_sum_lek$age <- gsub("chick", "Chicks", loads_sum_lek$age)
loads_sum_lek$loadtype <- gsub("gerp", "GERP", loads_sum_lek$loadtype)
loads_sum_lek$loadtype <- gsub("high", "SnpEff", loads_sum_lek$loadtype)
loads_sum_lek$age <- factor(loads_sum_lek$age, levels = c("Chicks", "Adults and yearlings"))

ggplot(subset(loads_sum_lek, loadtype=="GERP"), aes(x = site, y = mean, col = age, fill = age)) + 
  geom_point(shape=21, aes(size = n), position=position_dodge(.9)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  theme(legend.position="bottom")+
  scale_fill_manual(values=alpha(c(clrs_hunting[4], clrs_hunting[2]), 0.7))+
  scale_color_manual(values=alpha(c(clrs_hunting[4], clrs_hunting[2]), 0.7)) +
  labs(x = "Lek", y = "Total load", fill = "Age class", col = "Age class", size = "Sample size")+
  theme( plot.margin = margin(2,1,1,1, "cm"),panel.spacing = unit(3,"lines"))+
  guides(fill  = guide_legend(order = 1),
         col = guide_legend(order=1),
         size = guide_legend(order = 2))-> lek_size_gerp

lek_size_gerp

ggplot(subset(loads_sum_lek, loadtype=="SnpEff"), aes(x = site, y = mean, col = age, fill = age)) + 
  geom_point(shape=21, aes(size = n), position=position_dodge(.9)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  theme(legend.position="bottom")+
  scale_fill_manual(values=alpha(c(clrs_hunting[4], clrs_hunting[2]), 0.7))+
  scale_color_manual(values=alpha(c(clrs_hunting[4], clrs_hunting[2]), 0.7)) +
  labs(x = "Lek", y = "Total load", fill = "Age class", col = "Age class", size = "Sample size")+
  theme( plot.margin = margin(2,1,1,1, "cm"),panel.spacing = unit(3,"lines"))+
  guides(fill  = guide_legend(order = 1),
         col = guide_legend(order=1),
         size = guide_legend(order = 2))-> lek_size_high

lek_size_high

plot_grid(lek_size_gerp, lek_size_high,
          ncol = 1, 
          labels = c("A) GERP", "B) SnpEff"), label_fontface = "plain", label_size = 22) -> sup_lek
sup_lek

ggsave(sup_lek, file = "plots/plot_loads_lek.png", width = 12, height = 10)

#### test if effect is still there while excluding oldest individuals to have overlapping time frames ####

summary(lmerTest::lmer(scale(total_load) ~ age + (1|site), data = subset(loads, loadtype == "gerp" & year >= 2002)))
summary(lmerTest::lmer(scale(total_load) ~ age + (1|site), data = subset(loads, loadtype == "gerp")))
