#### packages ####
extrafont::loadfonts(device="all")
pacman::p_load(brms, bayesplot, tidyverse, data.table, ggridges, performance, scales)

#### theme ####
source("scripts/theme_ggplot.R")
source("scripts/function_diagnose_brms.R")

#### chick vs adult ####
#### inbreeding depression ####
load(file = "output/5_models/brm_froh_chicks.RData")

diagnose_froh_chick <- diagnose(fit = brm_froh_chick, modelname = "froh")

# get intervals
interval_froh_chick <- mcmc_intervals_data(brm_froh_chick, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
interval_froh_chick$model <- "froh"

# get areas
area_froh_chick <- mcmc_areas_data(brm_froh_chick, pars = "b_ageadult")
area_froh_chick$model <- "froh"

### plot

# split by interval
brms_froh_chick <- split(area_froh_chick, area_froh_chick$interval)

brms_froh_chick$bottom <- brms_froh_chick$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

ggplot(data = brms_froh_chick$outer) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density), fill = clr_froh, col = clr_froh)+
  geom_segment(data=interval_froh_chick, aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=interval_froh_chick, aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=interval_froh_chick, aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
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
brm_gerp_chick <- brm_gerp
load(file = "output/5_models/brm_high_total_chicks.RData")
brm_high_chick <- brm_high_noroh
# hom
# load(file = "output/5_models/brm_gerp_hom_chicks.RData")
# load(file = "output/5_models/brm_high_hom_chicks.RData")
# # het
# load(file = "output/5_models/brm_gerp_het_chicks.RData")
# load(file = "output/5_models/brm_high_het_chicks.RData")

##### diagnose ####
source("scripts/function_diagnose_brms.R")
diagnose_gerp_chick <- diagnose(fit = brm_gerp_chick, modelname = "gerp_total")
diagnose_high_chick <- diagnose(fit = brm_high_chick, modelname = "high_total")
# diagnose_gerp_hom <- diagnose(fit = brm_gerp_hom, modelname = "gerp_hom")
# diagnose_high_hom <- diagnose(fit = brm_high_hom, modelname = "high_hom")
# diagnose_gerp_het <- diagnose(fit = brm_gerp_het, modelname = "gerp_het")
# diagnose_high_het <- diagnose(fit = brm_high_het, modelname = "high_het")

### plot ###

# get intervals
gerp_interval_chick <- mcmc_intervals_data(brm_gerp_chick, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
high_interval_chick <-  mcmc_intervals_data(brm_high_chick, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
# gerp_hom_interval <- mcmc_intervals_data(brm_gerp_hom, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
# high_hom_interval <-  mcmc_intervals_data(brm_high_hom, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
# gerp_het_interval <- mcmc_intervals_data(brm_gerp_het, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
# high_het_interval <-  mcmc_intervals_data(brm_high_het, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")

intervals_load_chick <- rbind(gerp_interval_chick,
                   high_interval_chick
                   # gerp_hom_interval,
                   # high_hom_interval,
                   # gerp_het_interval,
                   # high_het_interval
                   )

intervals_load_chick$model <- c("GERP", "SnpEff")#, "GERP", "SnpEff", "GERP", "SnpEff")
intervals_load_chick$load <- c("Total", "Total")#, "Hom", "Hom", "Het", "Het")

# get areas
gerp_area_chick <- mcmc_areas_data(brm_gerp_chick, pars = "b_ageadult")
high_area_chick<- mcmc_areas_data(brm_high_chick, pars = "b_ageadult")
# gerp_hom_area <- mcmc_areas_data(brm_gerp_hom, pars = "b_ageadult")
# high_hom_area <- mcmc_areas_data(brm_high_hom, pars = "b_ageadult")
# gerp_het_area <- mcmc_areas_data(brm_gerp_het, pars = "b_ageadult")
# high_het_area <- mcmc_areas_data(brm_high_het, pars = "b_ageadult")

areas_load_chick <- rbind(gerp_area_chick,
               high_area_chick
               # gerp_hom_area,
               # high_hom_area,
               # gerp_het_area,
               # high_het_area
               )

areas_load_chick$model <- c(rep("GERP", nrow(gerp_area_chick)),
                 rep("SnpEff", nrow(high_area_chick))
                 # rep("GERP", nrow(gerp_hom_area)),
                 # rep("SnpEff", nrow(high_hom_area)),
                 # rep("GERP", nrow(gerp_het_area)),
                 # rep("SnpEff", nrow(high_het_area))
                 )

areas_load_chick$load <- c(rep("Total", nrow(gerp_area_chick)),
                rep("Total", nrow(high_area_chick))
                # rep("Hom", nrow(gerp_hom_area)),
                # rep("Hom", nrow(high_hom_area)),
                # rep("Het", nrow(gerp_het_area)),
                # rep("Het", nrow(high_het_area))
                )

#rearrange order for visualization
intervals_load_chick$model  <- factor(as.factor(intervals_load_chick$model),
                           levels= c("SnpEff", "GERP"))

# intervals_chick$load  <- factor(as.factor(intervals_chick$load),
#                           levels= c("Het", "Hom", "Total"))

areas_load_chick$model  <- factor(as.factor(areas_load_chick$model),
                       levels= c("SnpEff", "GERP"))

# areas_load_chick$load  <- factor(as.factor(areas_load_chick$load),
#                       levels= c("Het", "Hom", "Total"))

### plot

# split by interval
brms_loads_chick <- split(areas_load_chick, areas_load_chick$interval)

brms_loads_chick$bottom <- brms_loads_chick$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

##### plot total load ####

ggplot(data = subset(brms_loads_chick$outer, load == "Total")) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  geom_segment(data=subset(intervals_load_chick, load == "Total"), aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=subset(intervals_load_chick, load == "Total"), aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=subset(intervals_load_chick, load == "Total"), aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression(beta~"estimate for adults compared to chicks"), y = "Density", title = "Total load")+
  scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
  scale_color_manual(values =c(clr_high, clr_gerp)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> total_load_chicks

total_load_chicks

##### plot hom load ####
# 
# ggplot(data = subset(brms$outer, load == "Hom")) +  
#   aes(x = .data$x, y = .data$model) + 
#   geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
#   geom_segment(data=subset(intervals, load == "Hom"), aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
#   geom_segment(data=subset(intervals, load == "Hom"), aes(x = ll, xend = hh, yend = model), col = "black")+
#   geom_point(data=subset(intervals, load == "Hom"), aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
#   geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
#   labs(x = expression(beta~"estimate for adults compared to chicks"), y = "Density", title = "Homozygous load")+
#   scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
#   scale_color_manual(values =c(clr_high, clr_gerp)) +
#   theme(panel.border = element_blank(),
#         panel.grid = element_blank(),
#         strip.background = element_blank(),
#         legend.position = "none",
#         axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> hom_load
# 
# hom_load
# 
# 
# ##### plot het load ####
# 
# ggplot(data = subset(brms$outer, load == "Het")) +  
#   aes(x = .data$x, y = .data$model) + 
#   geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
#   geom_segment(data=subset(intervals, load == "Het"), aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
#   geom_segment(data=subset(intervals, load == "Het"), aes(x = ll, xend = hh, yend = model), col = "black")+
#   geom_point(data=subset(intervals, load == "Het"), aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
#   geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
#   labs(x = expression(beta~"estimate for adults compared to chicks"), y = "Density", title = "Heterozygous load")+
#   scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
#   scale_color_manual(values =c(clr_high, clr_gerp)) +
#   theme(panel.border = element_blank(),
#         panel.grid = element_blank(),
#         strip.background = element_blank(),
#         legend.position = "none",
#         axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> het_load
# 
# het_load
# 
# #### Combine figs ####
# cowplot::plot_grid(froh_chick,  total_load, hom_load,  het_load, 
#                    ncol = 2, align = "hv", axis = "lb",
#                    labels = "auto", 
#                    label_fontface = "plain", label_size = 22) -> fig
# 
# png("plots/froh_loads_chicks.png", height = 1000, width = 1000)
# fig
# dev.off()


### export intervals
intervals_load_chick_clean <- data.frame(parameter = intervals_load_chick$parameter,
                              model = intervals_load_chick$model, 
                              load = intervals_load_chick$load,
                              median = round(intervals_load_chick$m,2),
                              ci_95 = paste0(round(intervals_load_chick$ll, 2),", ", round(intervals_load_chick$hh,2)),
                              ci_80 = paste0(round(intervals_load_chick$l, 2),", ", round(intervals_load_chick$h,2)))
intervals_load_chick_clean

write_tsv(intervals_load_chick_clean, file = "output/intervals_chick_loads.tsv")

#### yearling vs adult #####
#### inbreeding depression ####
load(file = "output/5_models/brm_froh_yearling.RData")

diagnose_froh_yearling <- diagnose(fit = brm_froh_yearling, modelname = "froh")

# get intervals
interval_froh_yearling <- mcmc_intervals_data(brm_froh_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catadult")
interval_froh_yearling$model <- "froh"

# get areas
area_froh_yearling <- mcmc_areas_data(brm_froh_yearling, pars = "b_lifespan_catadult")
area_froh_yearling$model <- "froh"

### plot

# split by interval
brms_froh_yearling <- split(area_froh_yearling, area_froh_yearling$interval)

brms_froh_yearling$bottom <- brms_froh_yearling$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

ggplot(data = brms_froh_yearling$outer) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density), fill = clr_froh, col = clr_froh)+
  geom_segment(data=interval_froh_yearling, aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=interval_froh_yearling, aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=interval_froh_yearling, aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression(beta~"estimate for long- compared to short-lived males"), y = "Density", title = "Inbreeding")+
  scale_y_discrete(labels = c(expression(italic(F)[ROH])))+
  # xlim(-0.6, 0.6)+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        #axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) -> froh_yearling

froh_yearling

#### loads: total, hom and het ####
# total
load(file = "output/5_models/brm_gerp_total_yearlings.RData")
load(file = "output/5_models/brm_high_total_yearlings.RData")

# # hom
# load(file = "output/5_models/brm_gerp_hom_yearlings.RData")
# load(file = "output/5_models/brm_high_hom_yearlings.RData")
# # het
# load(file = "output/5_models/brm_gerp_het_yearlings.RData")
# load(file = "output/5_models/brm_high_het_yearlings.RData")

##### diagnose ####
source("scripts/function_diagnose_brms.R")
diagnose_gerp_yearling <- diagnose(fit = brm_gerp_yearling, modelname = "gerp_total_yearling")
diagnose_high_yearling <- diagnose(fit = brm_high_yearling, modelname = "high_total_yearling")
# diagnose_gerp_hom_yearling <- diagnose(fit = brm_gerp_hom_yearling, modelname = "gerp_hom_yearling")
# diagnose_high_hom_yearling <- diagnose(fit = brm_high_hom_yearling, modelname = "high_hom_yearling")
# diagnose_gerp_het_yearling <- diagnose(fit = brm_gerp_het_yearling, modelname = "gerp_het_yearling")
# diagnose_high_het_yearling <- diagnose(fit = brm_high_het_yearling, modelname = "high_het_yearling")

### plot ###

# get intervals
gerp_interval_yearling <- mcmc_intervals_data(brm_gerp_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catadult")
high_interval_yearling <-  mcmc_intervals_data(brm_high_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catadult")
# gerp_hom_interval_yearling <- mcmc_intervals_data(brm_gerp_hom_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catadult")
# high_hom_interval_yearling <-  mcmc_intervals_data(brm_high_hom_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catadult")
# gerp_het_interval_yearling <- mcmc_intervals_data(brm_gerp_het_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catadult")
# high_het_interval_yearling <-  mcmc_intervals_data(brm_high_het_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catadult")

intervals_yearling_load <- rbind(gerp_interval_yearling,
                   high_interval_yearling
                   # gerp_hom_interval_yearling,
                   # high_hom_interval_yearling,
                   # gerp_het_interval_yearling,
                   # high_het_interval_yearling
                   )

intervals_yearling_load$model <- c("GERP", "SnpEff")#, "GERP", "SnpEff", "GERP", "SnpEff")
intervals_yearling_load$load <- c("Total", "Total")#, "Hom", "Hom", "Het", "Het")

# get areas
gerp_area_yearling <- mcmc_areas_data(brm_gerp_yearling, pars = "b_lifespan_catadult")
high_area_yearling <- mcmc_areas_data(brm_high_yearling, pars = "b_lifespan_catadult")
# gerp_hom_area_yearling <- mcmc_areas_data(brm_gerp_hom_yearling, pars = "b_lifespan_catadult")
# high_hom_area_yearling <- mcmc_areas_data(brm_high_hom_yearling, pars = "b_lifespan_catadult")
# gerp_het_area_yearling <- mcmc_areas_data(brm_gerp_het_yearling, pars = "b_lifespan_catadult")
# high_het_area_yearling <- mcmc_areas_data(brm_high_het_yearling, pars = "b_lifespan_catadult")

areas_yearling_load <- rbind(gerp_area_yearling,
               high_area_yearling
               # gerp_hom_area_yearling,
               # high_hom_area_yearling,
               # gerp_het_area_yearling,
               # high_het_area_yearling
               )

areas_yearling_load$model <- c(rep("GERP", nrow(gerp_area_yearling)),
                 rep("SnpEff", nrow(high_area_yearling))
                 # rep("GERP", nrow(gerp_hom_area_yearling)),
                 # rep("SnpEff", nrow(high_hom_area_yearling)),
                 # rep("GERP", nrow(gerp_het_area_yearling)),
                 # rep("SnpEff", nrow(high_het_area_yearling))
                 )

areas_yearling_load$load <- c(rep("Total", nrow(gerp_area_yearling)),
                rep("Total", nrow(high_area_yearling))
                # rep("Hom", nrow(gerp_hom_area_yearling)),
                # rep("Hom", nrow(high_hom_area_yearling)),
                # rep("Het", nrow(gerp_het_area_yearling)),
                # rep("Het", nrow(high_het_area_yearling))
                )

#rearrange order for visualization
intervals_yearling_load$model  <- factor(as.factor(intervals_yearling_load$model),
                           levels= c("SnpEff", "GERP"))

# intervals_yearling_load$load  <- factor(as.factor(intervals_yearling_load$load),
#                           levels= c("Het", "Hom", "Total"))

areas_yearling_load$model  <- factor(as.factor(areas_yearling_load$model),
                       levels= c("SnpEff", "GERP"))

# areas_yearling_load$load  <- factor(as.factor(areas_yearling_load$load),
#                       levels= c("Het", "Hom", "Total"))

### plot

# split by interval
brms_yearling_load <- split(areas_yearling_load, areas_yearling_load$interval)

brms_yearling_load$bottom <- brms_yearling_load$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

##### plot total load ####

ggplot(data = subset(brms_yearling_load$outer, load == "Total")) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  geom_segment(data=subset(intervals_yearling_load, load == "Total"), aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=subset(intervals_yearling_load, load == "Total"), aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=subset(intervals_yearling_load, load == "Total"), aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression(beta~"estimate for long- compared to short-lived males"), y = "Density", title = "Total load")+
  scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
  scale_color_manual(values =c(clr_high, clr_gerp)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> total_load_yearling

total_load_yearling
# 
# ##### plot hom load ####
# 
# ggplot(data = subset(brms_yearling$outer, load == "Hom")) +  
#   aes(x = .data$x, y = .data$model) + 
#   geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
#   geom_segment(data=subset(intervals_yearling, load == "Hom"), aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
#   geom_segment(data=subset(intervals_yearling, load == "Hom"), aes(x = ll, xend = hh, yend = model), col = "black")+
#   geom_point(data=subset(intervals_yearling, load == "Hom"), aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
#   geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
#   labs(x = expression(beta~"estimate for long- compared to short-lived males"), y = "Density", title = "Homozygous load")+
#   scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
#   scale_color_manual(values =c(clr_high, clr_gerp)) +
#   theme(panel.border = element_blank(),
#         panel.grid = element_blank(),
#         strip.background = element_blank(),
#         legend.position = "none",
#         axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> hom_load
# 
# hom_load
# 
# 
# ##### plot het load ####
# 
# ggplot(data = subset(brms_yearling$outer, load == "Het")) +  
#   aes(x = .data$x, y = .data$model) + 
#   geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
#   geom_segment(data=subset(intervals_yearling, load == "Het"), aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
#   geom_segment(data=subset(intervals_yearling, load == "Het"), aes(x = ll, xend = hh, yend = model), col = "black")+
#   geom_point(data=subset(intervals_yearling, load == "Het"), aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
#   geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
#   labs(x = expression(beta~"estimate for long- compared to short-lived males"), y = "Density", title = "Heterozygous load")+
#   scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
#   scale_color_manual(values =c(clr_high, clr_gerp)) +
#   theme(panel.border = element_blank(),
#         panel.grid = element_blank(),
#         strip.background = element_blank(),
#         legend.position = "none",
#         axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> het_load
# 
# het_load
# 
# #### Combine figs ####
# cowplot::plot_grid(froh_yearling,  total_load, hom_load,  het_load, 
#                    ncol = 2, align = "hv", axis = "lb",
#                    labels = "auto", 
#                    label_fontface = "plain", label_size = 22) -> fig
# 
# png("plots/froh_loads_yearlings.png", height = 1000, width = 1000)
# fig
# dev.off()

### export intervals
intervals_yearling_load_clean <- data.frame(parameter = intervals_yearling_load$parameter,
                              model = intervals_yearling_load$model, 
                              load = intervals_yearling_load$load,
                              median = round(intervals_yearling_load$m,2),
                              ci_95 = paste0(round(intervals_yearling_load$ll, 2),", ", round(intervals_yearling_load$hh,2)),
                              ci_80 = paste0(round(intervals_yearling_load$l, 2),", ", round(intervals_yearling_load$h,2)))

intervals_yearling_load_clean

write_tsv(intervals_yearling_load_clean, file = "output/intervals_loads_yearling.tsv")


#### Create figure: separated by froh, total gerp load, total snpeff load but put chick and yearling togehter ####
intervals_load_chick$load <- NULL
intervals_yearling_load$load <- NULL

all_intervals <- rbind(interval_froh_chick,
                       interval_froh_yearling,
                       intervals_load_chick,
                       intervals_yearling_load)

areas_load_chick$load <- NULL
areas_yearling_load$load <- NULL
all_areas <- rbind(area_froh_chick,
                   area_froh_yearling,
                   areas_load_chick,
                   areas_yearling_load)

# rename
all_intervals$parameter <- gsub("b_ageadult", "Post-juveniles compared to chicks", all_intervals$parameter)
all_intervals$parameter <- gsub ("b_lifespan_catadult", "Adults compared to yearlings", all_intervals$parameter)

all_intervals$model <- gsub("froh", "Genomic inbreeding", all_intervals$model)
all_intervals$model <- gsub("GERP", "Total GERP load", all_intervals$model)
all_intervals$model <- gsub("SnpEff", "Total SnpEff load", all_intervals$model)


all_areas$parameter <- gsub("b_ageadult", "Post-juveniles compared to chicks", all_areas$parameter)
all_areas$parameter <- gsub ("b_lifespan_catadult", "Adults compared to yearlings", all_areas$parameter)

all_areas$model <- gsub("froh", "Genomic inbreeding", all_areas$model)
all_areas$model <- gsub("GERP", "Total GERP load", all_areas$model)
all_areas$model <- gsub("SnpEff", "Total SnpEff load", all_areas$model)

# prep figure
brms_all <- split(all_areas, all_areas$interval)

brms_all$bottom <- brms_all$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

froh_both <- ggplot(data = subset(brms_all$outer, model == "Genomic inbreeding")) +  
  aes(x = .data$x, y = .data$parameter) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = parameter, col = parameter))+
  geom_segment(data=subset(all_intervals, model == "Genomic inbreeding"), aes(x = l, xend = h, yend = parameter), col = "black", linewidth=3)+
  geom_segment(data=subset(all_intervals, model == "Genomic inbreeding"), aes(x = ll, xend = hh, yend = parameter), col = "black")+
  geom_point(data=subset(all_intervals, model == "Genomic inbreeding"), aes(x = m, y = parameter), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression("Standardised"~beta~"estimate"), title = expression(italic(F)[ROH]))+
  scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
  scale_color_manual(values =c(clr_high, clr_gerp)) +
  scale_y_discrete(labels = label_wrap(12)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank()) # element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) 

froh_both

gerp_both <- ggplot(data = subset(brms_all$outer, model == "Total GERP load")) +  
  aes(x = .data$x, y = .data$parameter) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = parameter, col = parameter))+
  geom_segment(data=subset(all_intervals, model == "Total GERP load"), aes(x = l, xend = h, yend = parameter), col = "black", linewidth=3)+
  geom_segment(data=subset(all_intervals, model == "Total GERP load"), aes(x = ll, xend = hh, yend = parameter), col = "black")+
  geom_point(data=subset(all_intervals, model == "Total GERP load"), aes(x = m, y = parameter), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression("Standardised"~beta~"estimate"), title = "Total GERP load")+
  scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
  scale_color_manual(values =c(clr_high, clr_gerp)) +
  scale_y_discrete(labels = label_wrap(12)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank()) # element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) 

gerp_both

snpeff_both <- ggplot(data = subset(brms_all$outer, model == "Total SnpEff load")) +  
  aes(x = .data$x, y = .data$parameter) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = parameter, col = parameter))+
  geom_segment(data=subset(all_intervals, model == "Total SnpEff load"), aes(x = l, xend = h, yend = parameter), col = "black", linewidth=3)+
  geom_segment(data=subset(all_intervals, model == "Total SnpEff load"), aes(x = ll, xend = hh, yend = parameter), col = "black")+
  geom_point(data=subset(all_intervals, model == "Total SnpEff load"), aes(x = m, y = parameter), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression("Standardised"~beta~"estimate"), title = "Total SnpEff load")+
  scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
  scale_color_manual(values =c(clr_high, clr_gerp)) +
  scale_y_discrete(labels = label_wrap(12)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank()) # element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) 

snpeff_both

combined_plot <- cowplot::plot_grid(froh_both, gerp_both, snpeff_both, 
                                    ncol = 1, align = "hv", axis = "lb",
                                    labels = "auto", 
                                    label_fontface = "plain", label_size = 22) -> fig
                                     
png("plots/figure_1_compare.png", height = 1000, width = 600)
fig
dev.off()
