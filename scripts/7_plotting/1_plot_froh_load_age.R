#### packages ####
extrafont::loadfonts(device="all")
pacman::p_load(brms, bayesplot, tidyverse, data.table, ggridges, performance)

#### theme ####
source("scripts/theme_ggplot.R")
source("scripts/function_diagnose_brms.R")

#### inbreeding depression ####
load(file = "output/5_models/brm_froh_chicks.RData")

diagnose_froh <- diagnose(fit = brm_froh_chick, modelname = "froh")

# get intervals
interval_froh <- mcmc_intervals_data(brm_froh_chick, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
interval_froh$model <- "froh"

# get areas
area_froh <- mcmc_areas_data(brm_froh_chick, pars = "b_ageadult")
area_froh$model <- "froh"

### plot

# split by interval
brms_froh <- split(area_froh, area_froh$interval)

brms_froh$bottom <- brms_froh$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

ggplot(data = brms_froh$outer) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density), fill = clr_froh, col = clr_froh)+
  geom_segment(data=interval_froh, aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=interval_froh, aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=interval_froh, aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression(beta~"estimate for adults compared to chicks"), y = "Density", title = "Inbreeding")+
  scale_y_discrete(labels = c(expression(italic(F)[ROH])))+
  # xlim(-0.6, 0.6)+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        #axis.title.y = element_blank(),
       # axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) -> froh_chick

froh_chick

#### loads: total, hom and het ####
# total
load(file = "output/5_models/brm_gerp_total_chicks.RData")
load(file = "output/5_models/brm_high_total_chicks.RData")
brm_high <- brm_high_noroh
# hom
load(file = "output/5_models/brm_gerp_hom_chicks.RData")
load(file = "output/5_models/brm_high_hom_chicks.RData")
# het
load(file = "output/5_models/brm_gerp_het_chicks.RData")
load(file = "output/5_models/brm_high_het_chicks.RData")

##### diagnose ####
source("scripts/function_diagnose_brms.R")
diagnose_gerp <- diagnose(fit = brm_gerp, modelname = "gerp_total")
diagnose_high <- diagnose(fit = brm_high, modelname = "high_total")
diagnose_gerp_hom <- diagnose(fit = brm_gerp_hom, modelname = "gerp_hom")
diagnose_high_hom <- diagnose(fit = brm_high_hom, modelname = "high_hom")
diagnose_gerp_het <- diagnose(fit = brm_gerp_het, modelname = "gerp_het")
diagnose_high_het <- diagnose(fit = brm_high_het, modelname = "high_het")

### plot ###

# get intervals
gerp_interval <- mcmc_intervals_data(brm_gerp, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
high_interval <-  mcmc_intervals_data(brm_high, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
gerp_hom_interval <- mcmc_intervals_data(brm_gerp_hom, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
high_hom_interval <-  mcmc_intervals_data(brm_high_hom, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
gerp_het_interval <- mcmc_intervals_data(brm_gerp_het, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
high_het_interval <-  mcmc_intervals_data(brm_high_het, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")

intervals <- rbind(gerp_interval,
                   high_interval,
                   gerp_hom_interval,
                   high_hom_interval,
                   gerp_het_interval,
                   high_het_interval)

intervals$model <- c("GERP", "SnpEff", "GERP", "SnpEff", "GERP", "SnpEff")
intervals$load <- c("Total", "Total", "Hom", "Hom", "Het", "Het")

# get areas
gerp_area <- mcmc_areas_data(brm_gerp, pars = "b_ageadult")
high_area <- mcmc_areas_data(brm_high, pars = "b_ageadult")
gerp_hom_area <- mcmc_areas_data(brm_gerp_hom, pars = "b_ageadult")
high_hom_area <- mcmc_areas_data(brm_high_hom, pars = "b_ageadult")
gerp_het_area <- mcmc_areas_data(brm_gerp_het, pars = "b_ageadult")
high_het_area <- mcmc_areas_data(brm_high_het, pars = "b_ageadult")

areas <- rbind(gerp_area,
               high_area,
               gerp_hom_area,
               high_hom_area,
               gerp_het_area,
               high_het_area)

areas$model <- c(rep("GERP", nrow(gerp_area)),
                 rep("SnpEff", nrow(high_area)),
                 rep("GERP", nrow(gerp_hom_area)),
                 rep("SnpEff", nrow(high_hom_area)),
                 rep("GERP", nrow(gerp_het_area)),
                 rep("SnpEff", nrow(high_het_area)))

areas$load <- c(rep("Total", nrow(gerp_area)),
                rep("Total", nrow(high_area)),
                rep("Hom", nrow(gerp_hom_area)),
                rep("Hom", nrow(high_hom_area)),
                rep("Het", nrow(gerp_het_area)),
                rep("Het", nrow(high_het_area)))

#rearrange order for visualization
intervals$model  <- factor(as.factor(intervals$model),
                           levels= c("SnpEff", "GERP"))

intervals$load  <- factor(as.factor(intervals$load),
                          levels= c("Het", "Hom", "Total"))

areas$model  <- factor(as.factor(areas$model),
                       levels= c("SnpEff", "GERP"))

areas$load  <- factor(as.factor(areas$load),
                      levels= c("Het", "Hom", "Total"))

### plot

# split by interval
brms <- split(areas, areas$interval)

brms$bottom <- brms$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

##### plot total load ####

ggplot(data = subset(brms$outer, load == "Total")) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  geom_segment(data=subset(intervals, load == "Total"), aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=subset(intervals, load == "Total"), aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=subset(intervals, load == "Total"), aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression(beta~"estimate for adults compared to chicks"), y = "Density", title = "Total load")+
  scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
  scale_color_manual(values =c(clr_high, clr_gerp)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> total_load

total_load

##### plot hom load ####

ggplot(data = subset(brms$outer, load == "Hom")) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  geom_segment(data=subset(intervals, load == "Hom"), aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=subset(intervals, load == "Hom"), aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=subset(intervals, load == "Hom"), aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression(beta~"estimate for adults compared to chicks"), y = "Density", title = "Homozygous load")+
  scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
  scale_color_manual(values =c(clr_high, clr_gerp)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> hom_load

hom_load


##### plot het load ####

ggplot(data = subset(brms$outer, load == "Het")) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  geom_segment(data=subset(intervals, load == "Het"), aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=subset(intervals, load == "Het"), aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=subset(intervals, load == "Het"), aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression(beta~"estimate for adults compared to chicks"), y = "Density", title = "Heterozygous load")+
  scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
  scale_color_manual(values =c(clr_high, clr_gerp)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> het_load

het_load

#### Combine figs ####
cowplot::plot_grid(froh_chick,  total_load, hom_load,  het_load, 
                   ncol = 2, align = "hv", axis = "lb",
                   labels = "auto", 
                   label_fontface = "plain", label_size = 22) -> fig

png("plots/froh_loads_chicks.png", height = 1000, width = 1000)
fig
dev.off()


### export intervals
intervals_clean <- data.frame(parameter = intervals$parameter,
                              model = intervals$model, 
                              load = intervals$load,
                              median = round(intervals$m,2),
                              ci_95 = paste0(round(intervals$ll, 2),",", round(intervals$hh,2)))
intervals_clean

write_tsv(intervals_clean, file = "output/intervals_loads.tsv")
