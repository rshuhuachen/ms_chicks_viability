### packages ###
extrafont::loadfonts(device="all")
pacman::p_load(tidyverse, data.table, brms, bayesplot, cowplot, ggsignif, ggpubr, reshape2, ggridges)

#### theme ####
source("scripts/theme_ggplot.R")
source("scripts/function_diagnose_brms.R")

#### Chick vs (sub)-adult ####
##### load model outputs ####
load(file = "output/5_models/per_region/brms_total_chicks_gerp_promoters.RData")
gerp_promo <- fit
load(file = "output/5_models/per_region/brms_total_chicks_gerp_exons.RData")
gerp_exon <- fit
load(file = "output/5_models/per_region/brms_total_chicks_gerp_introns.RData")
gerp_intron <- fit

load(file = "output/5_models/per_region/brms_total_chicks_high_promoters.RData")
high_promo <- fit
load(file = "output/5_models/per_region/brms_total_chicks_high_exons.RData")
high_exon <- fit
load(file = "output/5_models/per_region/brms_total_chicks_high_introns.RData")
high_intron <- fit

rm(fit)

##### diagnose ####
diagnose_gerp_promo <- diagnose(fit = gerp_promo, modelname = "gerp_promo")
diagnose_gerp_exon <- diagnose(fit = gerp_exon, modelname = "gerp_exon")
diagnose_gerp_intron <- diagnose(fit = gerp_intron, modelname = "gerp_intron")
diagnose_high_promo <- diagnose(fit = high_promo, modelname = "high_promo")
diagnose_high_exon <- diagnose(fit = high_exon, modelname = "high_exon")
diagnose_high_intron <- diagnose(fit = high_intron, modelname = "high_intron")

### plot ###

# get intervals
gerp_promo_interval <- mcmc_intervals_data(gerp_promo, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
gerp_exon_interval <-  mcmc_intervals_data(gerp_exon, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
gerp_intron_interval <- mcmc_intervals_data(gerp_intron, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
high_promo_interval <-  mcmc_intervals_data(high_promo, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
high_exon_interval <- mcmc_intervals_data(high_exon, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")
high_intron_interval <-  mcmc_intervals_data(high_intron, prob =0.8, prob_outer = 0.95, pars = "b_ageadult")

intervals <- rbind(gerp_promo_interval,
                   gerp_exon_interval,
                   gerp_intron_interval,
                   high_promo_interval,
                   high_exon_interval,
                   high_intron_interval)

intervals$model <- c("GERP", "GERP", "GERP", "SnpEff", "SnpEff", "SnpEff")
intervals$region <- c("Promoter", "Exon", "Intron", "Promoter", "Exon", "Intron")

# get areas
gerp_promo_area <- mcmc_areas_data(gerp_promo, pars = "b_ageadult")
gerp_exon_area <- mcmc_areas_data(gerp_exon, pars = "b_ageadult")
gerp_intron_area <- mcmc_areas_data(gerp_intron, pars = "b_ageadult")
high_promo_area <- mcmc_areas_data(high_promo, pars = "b_ageadult")
high_exon_area <- mcmc_areas_data(high_exon, pars = "b_ageadult")
high_intron_area <- mcmc_areas_data(high_intron, pars = "b_ageadult")

areas <- rbind(gerp_promo_area,
               gerp_exon_area,
               gerp_intron_area,
               high_promo_area,
               high_exon_area,
               high_intron_area)

areas$model <- c(rep("GERP", nrow(gerp_promo_area)),
                 rep("GERP", nrow(gerp_exon_area)),
                 rep("GERP", nrow(gerp_intron_area)),
                 rep("SnpEff", nrow(high_promo_area)),
                 rep("SnpEff", nrow(high_exon_area)),
                 rep("SnpEff", nrow(high_intron_area)))

areas$region <- c(rep("Promoter", nrow(gerp_promo_area)),
                rep("Exon", nrow(gerp_exon_area)),
                rep("Intron", nrow(gerp_intron_area)),
                rep("Promoter", nrow(high_promo_area)),
                rep("Exon", nrow(high_exon_area)),
                rep("Intron", nrow(high_intron_area)))

#rearrange order for visualization
intervals$model  <- factor(as.factor(intervals$model),
                           levels= c("SnpEff", "GERP"))

intervals$region  <- factor(as.factor(intervals$region),
                          levels= c("Exon", "Promoter","Intron"))

areas$model  <- factor(as.factor(areas$model),
                       levels= c("SnpEff", "GERP"))

areas$region  <- factor(as.factor(areas$region),
                        levels= c("Exon", "Promoter","Intron"))

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

ggplot(data = brms$outer) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  geom_segment(data=intervals, aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=intervals, aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=intervals, aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression(beta~" for yearlings and adults compared to chicks"), y = "Density")+
  scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
  scale_color_manual(values =c(clr_high, clr_gerp)) +
  facet_wrap(~region, scales="free")+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> total_load_region

total_load_region

ggsave(total_load_region, file = "plots/load/posterior_load_per_region.png", width = 14, height = 8)

### export intervals
intervals_clean <- data.frame(parameter = intervals$parameter,
                              model = intervals$model, 
                              region = intervals$region,
                              median = round(intervals$m,2),
                              ci_95 = paste0(round(intervals$ll, 2),",", round(intervals$hh,2)),
                              ci_80 = paste0(round(intervals$l, 2),",", round(intervals$h,2)),
                              conditional_r2 = c(paste0(round(as.data.frame(r2_bayes(gerp_promo))$R2[1], 2), " [", 
                                                        round(as.data.frame(r2_bayes(gerp_promo))$CI_low[1], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(gerp_promo))$CI_high[1], 2), "]"),
                                                 paste0(round(as.data.frame(r2_bayes(gerp_exon))$R2[1], 2), " [", 
                                                        round(as.data.frame(r2_bayes(gerp_exon))$CI_low[1], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(gerp_exon))$CI_high[1], 2), "]"),
                                                 paste0(round(as.data.frame(r2_bayes(gerp_intron))$R2[1], 2), " [", 
                                                        round(as.data.frame(r2_bayes(gerp_intron))$CI_low[1], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(gerp_intron))$CI_high[1], 2), "]"),
                                                 paste0(round(as.data.frame(r2_bayes(high_promo))$R2[1], 2), " [", 
                                                        round(as.data.frame(r2_bayes(high_promo))$CI_low[1], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(high_promo))$CI_high[1], 2), "]"),
                                                 paste0(round(as.data.frame(r2_bayes(high_exon))$R2[1], 2), " [", 
                                                        round(as.data.frame(r2_bayes(high_exon))$CI_low[1], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(high_exon))$CI_high[1], 2), "]"),
                                                 paste0(round(as.data.frame(r2_bayes(high_intron))$R2[1], 2), " [", 
                                                        round(as.data.frame(r2_bayes(high_intron))$CI_low[1], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(high_intron))$CI_high[1], 2), "]")),
                              marginal_r2 = c(paste0(round(as.data.frame(r2_bayes(gerp_promo))$R2[2], 2), " [", 
                                                        round(as.data.frame(r2_bayes(gerp_promo))$CI_low[2], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(gerp_promo))$CI_high[2], 2), "]"),
                                                 paste0(round(as.data.frame(r2_bayes(gerp_exon))$R2[2], 2), " [", 
                                                        round(as.data.frame(r2_bayes(gerp_exon))$CI_low[2], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(gerp_exon))$CI_high[2], 2), "]"),
                                                 paste0(round(as.data.frame(r2_bayes(gerp_intron))$R2[2], 2), " [", 
                                                        round(as.data.frame(r2_bayes(gerp_intron))$CI_low[2], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(gerp_intron))$CI_high[2], 2), "]"),
                                                 paste0(round(as.data.frame(r2_bayes(high_promo))$R2[2], 2), " [", 
                                                        round(as.data.frame(r2_bayes(high_promo))$CI_low[2], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(high_promo))$CI_high[2], 2), "]"),
                                                 paste0(round(as.data.frame(r2_bayes(high_exon))$R2[2], 2), " [", 
                                                        round(as.data.frame(r2_bayes(high_exon))$CI_low[2], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(high_exon))$CI_high[2], 2), "]"),
                                                 paste0(round(as.data.frame(r2_bayes(high_intron))$R2[2], 2), " [", 
                                                        round(as.data.frame(r2_bayes(high_intron))$CI_low[2], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(high_intron))$CI_high[2], 2), "]")))


intervals_clean

write_tsv(intervals_clean, file = "output/intervals_loads_region_chick.tsv")

#### Yearling vs adult ####
load(file = "output/5_models/per_region/brms_total_subadult_gerp_promoters.RData")
gerp_promo_yearling <- fit
load(file = "output/5_models/per_region/brms_total_subadult_gerp_exons.RData")
gerp_exon_yearling <- fit
load(file = "output/5_models/per_region/brms_total_subadult_gerp_introns.RData")
gerp_intron_yearling <- fit

load(file = "output/5_models/per_region/brms_total_subadult_high_promoters.RData")
high_promo_yearling <- fit
load(file = "output/5_models/per_region/brms_total_subadult_high_exons.RData")
high_exon_yearling <- fit
load(file = "output/5_models/per_region/brms_total_subadult_high_introns.RData")
high_intron_yearling <- fit

rm(fit)

#### diagnose ####
diagnose_gerp_promo_yearling <- diagnose(fit = gerp_promo_yearling, modelname = "gerp_promo_sa")
diagnose_gerp_exon_yearling <- diagnose(fit = gerp_exon_yearling, modelname = "gerp_exon_sa")
diagnose_gerp_intron_yearling <- diagnose(fit = gerp_intron_yearling, modelname = "gerp_intron_sa")
diagnose_high_promo_yearling <- diagnose(fit = high_promo_yearling, modelname = "high_promo_sa")
diagnose_high_exon_yearling <- diagnose(fit = high_exon_yearling, modelname = "high_exon_sa")
diagnose_high_intron_yearling <- diagnose(fit = high_intron_yearling, modelname = "high_intron_sa")

### plot ###

# get intervals
gerp_promo_interval_yearling <- mcmc_intervals_data(gerp_promo_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catYearling")
gerp_exon_interval_yearling <-  mcmc_intervals_data(gerp_exon_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catYearling")
gerp_intron_interval_yearling <- mcmc_intervals_data(gerp_intron_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catYearling")
high_promo_interval_yearling <-  mcmc_intervals_data(high_promo_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catYearling")
high_exon_interval_yearling <- mcmc_intervals_data(high_exon_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catYearling")
high_intron_interval_yearling <-  mcmc_intervals_data(high_intron_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catYearling")

intervals_yearling <- rbind(gerp_promo_interval_yearling,
                   gerp_exon_interval_yearling,
                   gerp_intron_interval_yearling,
                   high_promo_interval_yearling,
                   high_exon_interval_yearling,
                   high_intron_interval_yearling)

intervals_yearling$model <- c("GERP", "GERP", "GERP", "SnpEff", "SnpEff", "SnpEff")
intervals_yearling$region <- c("Promoter", "Exon", "Intron", "Promoter", "Exon", "Intron")

# get areas
gerp_promo_area_yearling <- mcmc_areas_data(gerp_promo_yearling, pars = "b_lifespan_catYearling")
gerp_exon_area_yearling <- mcmc_areas_data(gerp_exon_yearling, pars = "b_lifespan_catYearling")
gerp_intron_area_yearling <- mcmc_areas_data(gerp_intron_yearling, pars = "b_lifespan_catYearling")
high_promo_area_yearling <- mcmc_areas_data(high_promo_yearling, pars = "b_lifespan_catYearling")
high_exon_area_yearling <- mcmc_areas_data(high_exon_yearling, pars = "b_lifespan_catYearling")
high_intron_area_yearling <- mcmc_areas_data(high_intron_yearling, pars = "b_lifespan_catYearling")

areas_yearling <- rbind(gerp_promo_area_yearling,
               gerp_exon_area_yearling,
               gerp_intron_area_yearling,
               high_promo_area_yearling,
               high_exon_area_yearling,
               high_intron_area_yearling)

areas_yearling$model <- c(rep("GERP", nrow(gerp_promo_area_yearling)),
                 rep("GERP", nrow(gerp_exon_area_yearling)),
                 rep("GERP", nrow(gerp_intron_area_yearling)),
                 rep("SnpEff", nrow(high_promo_area_yearling)),
                 rep("SnpEff", nrow(high_exon_area_yearling)),
                 rep("SnpEff", nrow(high_intron_area_yearling)))

areas_yearling$region <- c(rep("Promoter", nrow(gerp_promo_area_yearling)),
                  rep("Exon", nrow(gerp_exon_area_yearling)),
                  rep("Intron", nrow(gerp_intron_area_yearling)),
                  rep("Promoter", nrow(high_promo_area_yearling)),
                  rep("Exon", nrow(high_exon_area_yearling)),
                  rep("Intron", nrow(high_intron_area_yearling)))

#rearrange order for visualization
intervals_yearling$model  <- factor(as.factor(intervals_yearling$model),
                           levels= c("SnpEff", "GERP"))

intervals_yearling$region  <- factor(as.factor(intervals_yearling$region),
                            levels= c("Exon", "Promoter","Intron"))

areas_yearling$model  <- factor(as.factor(areas_yearling$model),
                       levels= c("SnpEff", "GERP"))

areas_yearling$region  <- factor(as.factor(areas_yearling$region),
                        levels= c("Exon", "Promoter","Intron"))

### plot

# split by interval
brms_yearling <- split(areas_yearling, areas_yearling$interval)

brms_yearling$bottom <- brms_yearling$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

##### plot total load ####

ggplot(data = brms_yearling$outer) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  geom_segment(data=intervals_yearling, aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=intervals_yearling, aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=intervals_yearling, aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression(beta~" for adults compared to yearlings"), y = "Density")+
  scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
  scale_color_manual(values =c(clr_high, clr_gerp)) +
  facet_wrap(~region, scales="free")+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> total_load_region_sa

total_load_region_sa

ggsave(total_load_region_sa, file = "plots/load/posterior_load_per_region_subadult.png", width = 14, height = 8)

### export intervals
intervals_clean_yearling <- data.frame(parameter = intervals_yearling$parameter,
                              model = intervals_yearling$model, 
                              region = intervals_yearling$region,
                              median = round(intervals_yearling$m,2),
                              ci_95 = paste0(round(intervals_yearling$ll, 2),",", round(intervals_yearling$hh,2)),
                              ci_80 = paste0(round(intervals_yearling$l, 2),",", round(intervals_yearling$h,2)),
                              conditional_r2 = c(paste0(round(as.data.frame(r2_bayes(gerp_promo_yearling))$R2[1], 2), " [", 
                                                        round(as.data.frame(r2_bayes(gerp_promo_yearling))$CI_low[1], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(gerp_promo_yearling))$CI_high[1], 2), "]"),
                                                 paste0(round(as.data.frame(r2_bayes(gerp_exon_yearling))$R2[1], 2), " [", 
                                                        round(as.data.frame(r2_bayes(gerp_exon_yearling))$CI_low[1], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(gerp_exon_yearling))$CI_high[1], 2), "]"),
                                                 paste0(round(as.data.frame(r2_bayes(gerp_intron_yearling))$R2[1], 2), " [", 
                                                        round(as.data.frame(r2_bayes(gerp_intron_yearling))$CI_low[1], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(gerp_intron_yearling))$CI_high[1], 2), "]"),
                                                 paste0(round(as.data.frame(r2_bayes(high_promo_yearling))$R2[1], 2), " [", 
                                                        round(as.data.frame(r2_bayes(high_promo_yearling))$CI_low[1], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(high_promo_yearling))$CI_high[1], 2), "]"),
                                                 paste0(round(as.data.frame(r2_bayes(high_exon_yearling))$R2[1], 2), " [", 
                                                        round(as.data.frame(r2_bayes(high_exon_yearling))$CI_low[1], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(high_exon_yearling))$CI_high[1], 2), "]"),
                                                 paste0(round(as.data.frame(r2_bayes(high_intron_yearling))$R2[1], 2), " [", 
                                                        round(as.data.frame(r2_bayes(high_intron_yearling))$CI_low[1], 2), ", ", 
                                                        round(as.data.frame(r2_bayes(high_intron_yearling))$CI_high[1], 2), "]")),
                              marginal_r2 = c(paste0(round(as.data.frame(r2_bayes(gerp_promo_yearling))$R2[2], 2), " [", 
                                                     round(as.data.frame(r2_bayes(gerp_promo_yearling))$CI_low[2], 2), ", ", 
                                                     round(as.data.frame(r2_bayes(gerp_promo_yearling))$CI_high[2], 2), "]"),
                                              paste0(round(as.data.frame(r2_bayes(gerp_exon_yearling))$R2[2], 2), " [", 
                                                     round(as.data.frame(r2_bayes(gerp_exon_yearling))$CI_low[2], 2), ", ", 
                                                     round(as.data.frame(r2_bayes(gerp_exon_yearling))$CI_high[2], 2), "]"),
                                              paste0(round(as.data.frame(r2_bayes(gerp_intron_yearling))$R2[2], 2), " [", 
                                                     round(as.data.frame(r2_bayes(gerp_intron_yearling))$CI_low[2], 2), ", ", 
                                                     round(as.data.frame(r2_bayes(gerp_intron_yearling))$CI_high[2], 2), "]"),
                                              paste0(round(as.data.frame(r2_bayes(high_promo_yearling))$R2[2], 2), " [", 
                                                     round(as.data.frame(r2_bayes(high_promo_yearling))$CI_low[2], 2), ", ", 
                                                     round(as.data.frame(r2_bayes(high_promo_yearling))$CI_high[2], 2), "]"),
                                              paste0(round(as.data.frame(r2_bayes(high_exon_yearling))$R2[2], 2), " [", 
                                                     round(as.data.frame(r2_bayes(high_exon_yearling))$CI_low[2], 2), ", ", 
                                                     round(as.data.frame(r2_bayes(high_exon_yearling))$CI_high[2], 2), "]"),
                                              paste0(round(as.data.frame(r2_bayes(high_intron_yearling))$R2[2], 2), " [", 
                                                     round(as.data.frame(r2_bayes(high_intron_yearling))$CI_low[2], 2), ", ", 
                                                     round(as.data.frame(r2_bayes(high_intron_yearling))$CI_high[2], 2), "]")))
intervals_clean_yearling

write_tsv(intervals_clean_yearling, file = "output/intervals_loads_region_yearling.tsv")


#### Combine chick and yearling GERP posteriors ####

intervals_both <- rbind(intervals, intervals_yearling)
areas_both <- rbind(areas, areas_yearling)

### plot

# split by interval
brms_both <- split(areas_both, areas_both$interval)

brms_both$bottom <- brms_both$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

##### plot total GERP load only
brms_both$outer$parameter <- gsub ("b_ageadult", "Post-juvenile compared to chick", brms_both$outer$parameter )
brms_both$outer$parameter <- gsub ("b_lifespan_catAdult", "Adult compared to yearling", brms_both$outer$parameter )

intervals_both$parameter <- gsub ("b_ageadult", "Post-juvenile compared to chick", intervals_both$parameter)
intervals_both$parameter <- gsub ("b_lifespan_catAdult", "Adult compared to yearling", intervals_both$parameter)
# relevel
intervals_both$parameter <- factor(intervals_both$parameter, levels = c("Post-juvenile compared to chick", "Adult compared to yearling"))
brms_both$outer$parameter <- factor(brms_both$outer$parameter, levels = c("Post-juvenile compared to chick", "Adult compared to yearling"))

intervals_both$region <- factor(intervals_both$region, levels = c("Intron", "Promoter", "Exon"))
brms_both$outer$region <- factor(brms_both$outer$region, levels = c("Intron", "Promoter", "Exon"))


ggplot(data = subset(brms_both$outer, brms_both$outer$model == "GERP")) +  
  aes(x = .data$x, y = .data$region) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = parameter, col = parameter))+
  geom_segment(data=subset(intervals_both, model == "GERP"), aes(x = l, xend = h, yend = region), col = "black", linewidth=3)+
  geom_segment(data=subset(intervals_both, model == "GERP"), aes(x = ll, xend = hh, yend = region), col = "black")+
  geom_point(data=subset(intervals_both, model == "GERP"), aes(x = m, y = region), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression("Standardised"~beta~" estimate"))+
  scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
  scale_color_manual(values =c(clr_high, clr_gerp)) +
  facet_grid(~parameter, labeller = label_wrap_gen())+
  xlim(-1.5,1)+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank()) -> gerp_per_region_both
gerp_per_region_both
ggsave(gerp_per_region_both, file = "plots/figure_2_posteriors_gerp_region.png", width = 12, height = 8)

#### Boxplots ####
#### Load loads per region ####
load("output/loads_per_region.RData")

#### Load metadata ####
load("data/metadata/metadata_adult_chick.RData")

meta <- meta %>% mutate(age = case_when(
  grepl("C", meta$id) ~ "Chick",
  grepl("D", meta$id) ~ "Adult"
))

meta$age <- factor(meta$age, levels = c("Chick", "Adult"))

### Merge
load_per_region <- left_join(load_per_region, meta, by = "id")

summary(lmerTest::lmer(total_load ~ age + (1|year), data = subset(load_per_region, loadtype == "gerp_exons")))
summary(lmerTest::lmer(total_load ~ age + (1|year), data = subset(load_per_region, loadtype == "gerp_promoters")))
summary(lmerTest::lmer(total_load ~ age + (1|year), data = subset(load_per_region, loadtype == "gerp_introns")))


summary(lmerTest::lmer(total_load ~ age + (1|year), data = subset(load_per_region, loadtype == "high_exons")))
summary(lmerTest::lmer(total_load ~ age + (1|year), data = subset(load_per_region, loadtype == "high_promoters")))
summary(lmerTest::lmer(total_load ~ age + (1|year), data = subset(load_per_region, loadtype == "high_introns")))

#### Merge all and plot ####

load_per_region <- load_per_region %>% mutate(
  method = case_when(grepl("gerp", load_per_region$loadtype)~ "GERP ≥ 4",
                     grepl("high", load_per_region$loadtype)~ "High impact SnpEff"))

load_per_region <- load_per_region %>% mutate(
  region = case_when(grepl("introns", load_per_region$loadtype)~ "Introns",
                     grepl("promoters", load_per_region$loadtype)~ "Promoters",
                     grepl("exon", load_per_region$loadtype)~ "Exons"))

#### Plot ####
ggplot(load_per_region, aes(x = age, y = total_load)) + 
  geom_point(position="jitter", aes(fill = method, color = method))+
  geom_boxplot(fill = alpha("white", 0.7), outlier.shape=NA) + 
  theme(plot.title = element_text(margin=margin(0,0,1,0, "cm")))+
  scale_color_manual(values = c(clr_gerp, clr_high))+
  scale_fill_manual(values = c(clr_gerp, clr_high))+
  labs(x = "Age", y = "Total load") +
  facet_wrap(method ~ region, scales="free") +
  theme(legend.position="none") -> compare_loads

compare_loads
ggsave(compare_loads, file = "plots/load/boxplots_region_gerp_snpeff.png", width=14,height=18)

### Figure 2 alternative: posteriors and raw data ####

# posteriors
ggplot(data = subset(brms_both$outer, brms_both$outer$model == "GERP" & parameter == "Post-juvenile compared to chick")) +  
  aes(x = .data$x, y = .data$region) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = parameter, col = parameter))+
  geom_segment(data=subset(intervals_both, model == "GERP"& parameter == "Post-juvenile compared to chick"), aes(x = l, xend = h, yend = region), col = "black", linewidth=3)+
  geom_segment(data=subset(intervals_both, model == "GERP"& parameter == "Post-juvenile compared to chick"), aes(x = ll, xend = hh, yend = region), col = "black")+
  geom_point(data=subset(intervals_both, model == "GERP"& parameter == "Post-juvenile compared to chick"), aes(x = m, y = region), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression("Standardised"~beta~" estimate"))+
  scale_fill_manual(values =alpha(c(clr_high), 0.7)) +
  scale_color_manual(values =c(clr_high)) +
  xlim(-1.5,1)+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank()) -> gerp_per_region_yearling

# raw data
load_per_region$region <- gsub("Exons", "Exon", load_per_region$region)
load_per_region$region <- gsub("Promoters", "Promoter", load_per_region$region)
load_per_region$region <- gsub("Introns", "Intron", load_per_region$region)
load_per_region$region <- factor(load_per_region$region, levels = c("Exon", "Promoter", "Intron"))

load_per_region$age <- gsub("Adult", "Yearling", load_per_region$age)

subset(load_per_region, method == "GERP ≥ 4") %>% 
  ggplot(aes(x = age, y = total_load)) + 
  geom_point(position="jitter", aes(fill = method), color="black", shape=21, size = 3)+
  geom_violin(fill = alpha("white", 0.7)) + 
  theme(plot.title = element_text(margin=margin(0,0,1,0, "cm")))+
  geom_boxplot(width=0.2, fill = alpha("#8DAA91", 0.7), outlier.shape = NA)+
#  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "black")+
  scale_color_manual(values = c(clr_high))+
  scale_fill_manual(values = c(clr_high))+
  labs(x = "Age", y = "Total GERP load") +
  facet_grid(region~., scales="free") +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank()) -> rawdata_gerp_chick 

plot_grid(gerp_per_region_yearling, rawdata_gerp_chick, 
          ncol = 2, align = "hv", axis = "lb",
          labels = c("A", "B"), label_fontface = "plain", label_size = 22) -> fig2
fig2
ggsave(fig2, file = "plots/figure_2_alternative_with_raw.png", width=14,height=10)

