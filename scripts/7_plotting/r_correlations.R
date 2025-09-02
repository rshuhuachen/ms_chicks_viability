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

ggplot(loads_sum, aes(x = year, y = mean, col = age, fill = age)) + 
  geom_point(shape=21, aes(size = n), position=position_dodge(.9)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~loadtype, ncol=1, scale="free") +
  theme(legend.position="bottom")+
  scale_fill_manual(values=alpha(c(clrs_hunting[4], clrs_hunting[2]), 0.7))+
  scale_color_manual(values=alpha(c(clrs_hunting[4], clrs_hunting[2]), 0.7)) +
  labs(x = "Birth year", y = "Total load", fill = "Age class", col = "Age class", size = "Sample size") -> time_size

time_size

ggsave(time_size, file = "plots/plot_loads_time.png", width = 10, height = 8)

### differences among leks  


loads_sum_lek <- loads %>% group_by(age, loadtype, site) %>% summarise(
  n = n(),
  mean = mean(total_load, na.rm=TRUE),
  sd = sd(total_load))

loads_sum_lek$age <- gsub("adult", "Adults and yearlings", loads_sum_lek$age)
loads_sum_lek$age <- gsub("chick", "Chicks", loads_sum_lek$age)
loads_sum_lek$loadtype <- gsub("gerp", "GERP", loads_sum_lek$loadtype)
loads_sum_lek$loadtype <- gsub("high", "SnpEff", loads_sum_lek$loadtype)

ggplot(loads_sum_lek, aes(x = site, y = mean, col = age, fill = age)) + 
  geom_point(shape=21, aes(size = n), position=position_dodge(.9)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~loadtype, ncol=1, scale="free") +
  theme(legend.position="bottom")+
  scale_fill_manual(values=alpha(c(clrs_hunting[4], clrs_hunting[2]), 0.7))+
  scale_color_manual(values=alpha(c(clrs_hunting[4], clrs_hunting[2]), 0.7)) +
  labs(x = "Lek", y = "Total load", fill = "Age class", col = "Age class", size = "Sample size") -> lek_size

lek_size

ggsave(lek_size, file = "plots/plot_loads_lek.png", width = 10, height = 8)

#### test if effect is still there while excluding oldest individuals to have overlapping time frames ####

summary(lmerTest::lmer(scale(total_load) ~ age + (1|site), data = subset(loads, loadtype == "gerp" & year >= 2002)))
summary(lmerTest::lmer(scale(total_load) ~ age + (1|site), data = subset(loads, loadtype == "gerp")))
