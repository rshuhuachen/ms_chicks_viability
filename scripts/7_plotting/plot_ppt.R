#### packages ####
extrafont::loadfonts(device="all")
pacman::p_load(brms, bayesplot, tidyverse, data.table, ggridges, performance, scales)

#### theme ####
source("scripts/theme_ggplot.R")

#### chick vs adult ####
#### inbreeding depression ####
load(file = "output/5_models/brm_froh_chicks.RData")

# get intervals
interval_froh_chick <- mcmc_intervals_data(brm_froh_chick, prob =0.8, prob_outer = 0.95, pars = "b_agechick")
interval_froh_chick$model <- "froh"

# get areas
area_froh_chick <- mcmc_areas_data(brm_froh_chick, pars = "b_agechick")
area_froh_chick$model <- "froh"

#### loads: total, hom and het ####
# total
load(file = "output/5_models/brm_gerp_total_chicks.RData")
brm_gerp_chick <- brm_gerp
load(file = "output/5_models/brm_high_total_chicks.RData")
brm_high_chick <- brm_high

### plot ###

# get intervals
gerp_interval_chick <- mcmc_intervals_data(brm_gerp_chick, prob =0.8, prob_outer = 0.95, pars = "b_agechick")
high_interval_chick <-  mcmc_intervals_data(brm_high_chick, prob =0.8, prob_outer = 0.95, pars = "b_agechick")

intervals_load_chick <- rbind(gerp_interval_chick,
                              high_interval_chick)

intervals_load_chick$model <- c("GERP", "SnpEff")
intervals_load_chick$load <- c("Total", "Total")

# get areas
gerp_area_chick <- mcmc_areas_data(brm_gerp_chick, pars = "b_agechick")
high_area_chick<- mcmc_areas_data(brm_high_chick, pars = "b_agechick")

areas_load_chick <- rbind(gerp_area_chick,
                          high_area_chick)

areas_load_chick$model <- c(rep("GERP", nrow(gerp_area_chick)),
                            rep("SnpEff", nrow(high_area_chick)))

areas_load_chick$load <- c(rep("Total", nrow(gerp_area_chick)),
                           rep("Total", nrow(high_area_chick)))

#rearrange order for visualization
intervals_load_chick$model  <- factor(as.factor(intervals_load_chick$model),
                                      levels= c("SnpEff", "GERP"))

areas_load_chick$model  <- factor(as.factor(areas_load_chick$model),
                                  levels= c("SnpEff", "GERP"))


#### yearling vs adult #####
#### inbreeding depression ####
load(file = "output/5_models/brm_froh_yearling.RData")

#diagnose_froh_yearling <- diagnose(fit = brm_froh_yearling, modelname = "froh")

# get intervals
interval_froh_yearling <- mcmc_intervals_data(brm_froh_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catyearling")
interval_froh_yearling$model <- "froh"

# get areas
area_froh_yearling <- mcmc_areas_data(brm_froh_yearling, pars = "b_lifespan_catyearling")
area_froh_yearling$model <- "froh"


#### loads ####
load(file = "output/5_models/brm_gerp_total_yearlings.RData")
load(file = "output/5_models/brm_high_total_yearlings.RData")

### plot ###

# get intervals
gerp_interval_yearling <- mcmc_intervals_data(brm_gerp_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catyearling")
high_interval_yearling <-  mcmc_intervals_data(brm_high_yearling, prob =0.8, prob_outer = 0.95, pars = "b_lifespan_catyearling")

intervals_yearling_load <- rbind(gerp_interval_yearling,
                                 high_interval_yearling)

intervals_yearling_load$model <- c("GERP", "SnpEff")
intervals_yearling_load$load <- c("Total", "Total")

# get areas
gerp_area_yearling <- mcmc_areas_data(brm_gerp_yearling, pars = "b_lifespan_catyearling")
high_area_yearling <- mcmc_areas_data(brm_high_yearling, pars = "b_lifespan_catyearling")

areas_yearling_load <- rbind(gerp_area_yearling,
                             high_area_yearling)

areas_yearling_load$model <- c(rep("GERP", nrow(gerp_area_yearling)),
                               rep("SnpEff", nrow(high_area_yearling)))

areas_yearling_load$load <- c(rep("Total", nrow(gerp_area_yearling)),
                              rep("Total", nrow(high_area_yearling)))

#rearrange order for visualization
intervals_yearling_load$model  <- factor(as.factor(intervals_yearling_load$model),
                                         levels= c("SnpEff", "GERP"))

areas_yearling_load$model  <- factor(as.factor(areas_yearling_load$model),
                                     levels= c("SnpEff", "GERP"))


#### Create figure: separated by froh, total gerp load, total snpeff load but put chick and yearling togehter ####
intervals_load_chick$load <- NULL
intervals_yearling_load$load <- NULL

all_intervals <- rbind(interval_froh_chick,
                       interval_froh_yearling,
                       intervals_load_chick,
                       intervals_yearling_load)

all_intervals <- all_intervals %>% arrange(parameter)


areas_load_chick$load <- NULL
areas_yearling_load$load <- NULL
all_areas <- rbind(area_froh_chick,
                   area_froh_yearling,
                   areas_load_chick,
                   areas_yearling_load)

# rename
all_intervals$parameter <- gsub("b_agechick", "Chicks compared to yearlings/adults", all_intervals$parameter)
all_intervals$parameter <- gsub ("b_lifespan_catyearling", "Yearlings compared to adults", all_intervals$parameter)

all_intervals$model <- gsub("froh", "Genomic inbreeding", all_intervals$model)
all_intervals$model <- gsub("GERP", "Total GERP load", all_intervals$model)
all_intervals$model <- gsub("SnpEff", "Total SnpEff load", all_intervals$model)


all_areas$parameter <- gsub("b_agechick", "Chicks compared to yearlings/adults", all_areas$parameter)
all_areas$parameter <- gsub ("b_lifespan_catyearling", "Yearlings compared to adults", all_areas$parameter)

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

pacman::p_load(scales)


### plot for ppt ####
brms_all$outer$model <- factor(brms_all$outer$model, levels = c("Genomic inbreeding", "Total SnpEff load", "Total GERP load"))
all_intervals$model <- factor(all_intervals$model, levels = c("Genomic inbreeding", "Total SnpEff load", "Total GERP load"))

brms_all$outer$parameter <- factor(brms_all$outer$parameter, levels = c("Chicks compared to yearlings/adults", "Yearlings compared to adults"))
all_intervals$parameter <- factor(all_intervals$parameter, levels = c("Chicks compared to yearlings/adults", "Yearlings compared to adults"))

all_posteriors_a <- ggplot(data = subset(brms_all$outer, parameter == "Chicks compared to yearlings/adults")) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  # geom_segment(data=subset(all_intervals, parameter == "Chicks compared to yearlings/adults"),  aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=subset(all_intervals, parameter == "Chicks compared to yearlings/adults"),  
               aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=subset(all_intervals, parameter == "Chicks compared to yearlings/adults"),  
               aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=subset(all_intervals, parameter == "Chicks compared to yearlings/adults"),  
             aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression("Standardised"~beta))+
  scale_fill_manual(values =alpha(c(clr_froh, clr_high, clr_gerp), 0.7)) +
  scale_color_manual(values =c(clr_froh, clr_high, clr_gerp)) +
  xlim(-1.5,1.5)+
  scale_y_discrete(labels = c(expression(italic(F)[ROH]), "Total SnpEff load", "Total GERP load"))+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill='transparent'),
        strip.background = element_blank(),
        legend.position = "none",
        plot.title = element_blank(),
        strip.text.x = element_text(size = 22, margin = margin(15,15, 200, 15)),
        axis.title.y = element_blank()) 

all_posteriors_a

png("plots/figure_1a.png", width=500, height=400)
all_posteriors_a
dev.off()

all_posteriors_b <- ggplot(data = subset(brms_all$outer, parameter == "Yearlings compared to adults")) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = parameter, col = parameter))+
  # geom_segment(data=subset(all_intervals, parameter == "Yearlings compared to adults"),  aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=subset(all_intervals, parameter == "Yearlings compared to adults"), 
               aes(x = ll, xend = hh, yend = model), col = "black", linewidth=1.5, 
               position=position_nudge(y = -0.1))+
  geom_point(data=subset(all_intervals, parameter == "Yearlings compared to adults"), aes(x = m, y = model), 
             fill="white",  col = "black", shape=21, size = 6,
             position=position_nudge(y = -0.1)) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression("Standardised"~beta~"estimate"))+
  scale_fill_manual(values =alpha(c(clr_gerp), 0.7)) +
  scale_color_manual(values =c(clr_gerp)) +
  xlim(-1.5,1.5)+
  scale_y_discrete(labels = c("Total SnpEff load", "Total GERP load", expression(italic(F)[ROH])))+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        plot.title = element_blank(),
        strip.text.x = element_text(size = 22, margin = margin(15,15, 200, 15)),
        axis.title.y = element_blank()) 
all_posteriors_b

all_posteriors <- cowplot::plot_grid(all_posteriors_a + theme(plot.margin = margin(10,1,1,1, "cm")), 
                                     all_posteriors_b + theme(plot.margin = margin(10,1,1,1, "cm")), 
                                     ncol = 2, align = "hv", axis = "lb",
                                     labels = c("A", "B"), 
                                     label_fontface = "bold", label_size = 22)

png("plots/figure_1_all.png", height = 800, width = 1000)
all_posteriors
dev.off()

ggsave(all_posteriors, file = "plots/figures_1_all.pdf", device=cairo_pdf, height=10,width=16)
ggsave(all_posteriors, file = "plots/figures_1_all.eps", device=cairo_pdf, height=10,width=16)
