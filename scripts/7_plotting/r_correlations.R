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
