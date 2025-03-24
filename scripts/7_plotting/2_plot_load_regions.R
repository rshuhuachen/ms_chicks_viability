
### Load loads per region ####
load("output/4_load/loads_per_region.RData")

### Load metadata ####
load("metadata/metadata_adult_chick.RData")

meta <- meta %>% mutate(age = case_when(
  grepl("C", meta$id) ~ "Chick",
  grepl("D", meta$id) ~ "Adult"
))

meta$age <- factor(meta$age, levels = c("Chick", "Adult"))

### Load phenotypes adults ####
load("metadata/phenotypes_wide_extra.RData")

### Merge
loads_regions <- left_join(loads_regions, meta, by = "id")

summary(lmerTest::lmer(total_load ~ age + (1|year), data = subset(loads_regions, loadtype == "exons")))
summary(lmerTest::lmer(total_load ~ age + (1|year), data = subset(loads_regions, loadtype == "promoters")))
summary(lmerTest::lmer(total_load ~ age + (1|year), data = subset(loads_regions, loadtype == "introns")))


summary(lmerTest::lmer(total_load ~ age + (1|year), data = subset(loads_regions_high, loadtype == "exons")))
summary(lmerTest::lmer(total_load ~ age + (1|year), data = subset(loads_regions_high, loadtype == "promoters")))
summary(lmerTest::lmer(total_load ~ age + (1|year), data = subset(loads_regions_high, loadtype == "introns")))

#### Merge all and plot ####
loads_regions_high$method <- "high"
loads_regions$method <- "gerp"

all_loads <- rbind(loads_regions_high, loads_regions)
save(all_loads, file = "output/loads_all_per_region.RData")

source("scripts/theme_ggplot.R")
pacman::p_load(ggsignif,ggpubr)

all_loads$loadtype <- gsub("exons", "Exons", all_loads$loadtype)
all_loads$loadtype <- gsub("introns", "Introns", all_loads$loadtype)
all_loads$loadtype <- gsub("promoters", "Promoters", all_loads$loadtype)
all_loads$method <- gsub("gerp", "GERP ≥ 4", all_loads$method)
all_loads$method <- gsub("high", "High impact SnpEff", all_loads$method)


### gerp - exon
aov_gerp_exon <- summary(aov(total_load ~ age, data = subset(all_loads, loadtype == "Exons" & method == "GERP ≥ 4")))
ggplot(subset(all_loads, loadtype == "Exons" & method == "GERP ≥ 4"), aes(x = age, y = total_load)) + 
  geom_point(position="jitter")+
  geom_boxplot(fill = alpha("white", 0.7), outlier.shape=NA) + 
  theme(plot.title = element_text(margin=margin(0,0,1,0, "cm")))+
  labs(x = "Age", y = "Total load", title = "GERP ≥ 4 - exons")+
  annotate("text", x = 1.5, y = 0.12, 
           label = "ANOVA: p-value < 0.001", size=6)-> compare_loads_gerp_exon
compare_loads_gerp_exon

### gerp - intron
aov_gerp_intron <- summary(aov(total_load ~ age, data = subset(all_loads, loadtype == "Introns" & method == "GERP ≥ 4")))
ggplot(subset(all_loads, loadtype == "Introns" & method == "GERP ≥ 4"), aes(x = age, y = total_load)) + 
  geom_point(position="jitter")+
  geom_boxplot(fill = alpha("white", 0.7), outlier.shape=NA) + 
  theme(plot.title = element_text(margin=margin(0,0,1,0, "cm")))+
  labs(x = "Age", y = "Total load", title = "GERP ≥ 4 - introns")+
  annotate("text", x = 1.5, y = 0.165, 
           label = paste0("ANOVA: p-value = ", round(aov_gerp_intron[[1]]$`Pr(>F)`[1], 3)), size=6)-> compare_loads_gerp_intron

compare_loads_gerp_intron

### gerp - promo
aov_gerp_Promoter <- summary(aov(total_load ~ age, data = subset(all_loads, loadtype == "Promoters" & method == "GERP ≥ 4")))
ggplot(subset(all_loads, loadtype == "Promoters" & method == "GERP ≥ 4"), aes(x = age, y = total_load)) + 
  geom_point(position="jitter")+
  geom_boxplot(fill = alpha("white", 0.7), outlier.shape=NA) + 
  theme(plot.title = element_text(margin=margin(0,0,1,0, "cm")))+
  labs(x = "Age", y = "Total load", title = "GERP ≥ 4 - promoters")+
  annotate("text", x = 1.5, y = 0.16, 
           label = paste0("ANOVA: p-value = ", round(aov_gerp_Promoter[[1]]$`Pr(>F)`[1], 3)), size=6)-> compare_loads_gerp_Promoter

compare_loads_gerp_Promoter

### high - exon
aov_high_exon <- summary(aov(total_load ~ age, data = subset(all_loads, loadtype == "Exons" & method == "High impact SnpEff")))
ggplot(subset(all_loads, loadtype == "Exons" & method == "High impact SnpEff"), aes(x = age, y = total_load)) + 
  geom_point(position="jitter")+
  geom_boxplot(fill = alpha("white", 0.7), outlier.shape=NA) + 
  theme(plot.title = element_text(margin=margin(0,0,1,0, "cm")))+
  labs(x = "Age", y = "Total load", title = "High impact SnpEff - exons")+
  annotate("text", x = 1.5, y = 0.161, 
           label = paste0("ANOVA: p-value = ", round(aov_high_exon[[1]]$`Pr(>F)`[1], 3)), size=6)-> compare_loads_high_exon
compare_loads_high_exon

### high - intron
aov_high_intron <- summary(aov(total_load ~ age, data = subset(all_loads, loadtype == "Introns" & method == "High impact SnpEff")))
ggplot(subset(all_loads, loadtype == "Introns" & method == "High impact SnpEff"), aes(x = age, y = total_load)) + 
  geom_point(position="jitter")+
  geom_boxplot(fill = alpha("white", 0.7), outlier.shape=NA) + 
  theme(plot.title = element_text(margin=margin(0,0,1,0, "cm")))+
  labs(x = "Age", y = "Total load", title = "High impact SnpEff - introns")+
  annotate("text", x = 1.5, y = 0.20, 
           label = paste0("ANOVA: p-value = ", round(aov_high_intron[[1]]$`Pr(>F)`[1], 3)), size=6)-> compare_loads_high_intron

compare_loads_high_intron

### high - promo
aov_high_Promoter <- summary(aov(total_load ~ age, data = subset(all_loads, loadtype == "Promoters" & method == "High impact SnpEff")))
ggplot(subset(all_loads, loadtype == "Promoters" & method == "High impact SnpEff"), aes(x = age, y = total_load)) + 
  geom_point(position="jitter")+
  geom_boxplot(fill = alpha("white", 0.7), outlier.shape=NA) + 
  theme(plot.title = element_text(margin=margin(0,0,1,0, "cm")))+
  labs(x = "Age", y = "Total load", title = "High impact SnpEff - promoters")+
  annotate("text", x = 1.5, y = 0.171, 
           label = paste0("ANOVA: p-value = ", round(aov_high_Promoter[[1]]$`Pr(>F)`[1], 3)), size=6)-> compare_loads_high_Promoter

compare_loads_high_Promoter


cowplot::plot_grid(compare_loads_gerp_exon, compare_loads_high_exon, compare_loads_gerp_intron, 
                   compare_loads_high_intron, compare_loads_gerp_Promoter,
                   compare_loads_high_Promoter,
                   ncol = 2, align = "hv", axis = "lb",
                   labels = "auto", label_fontface = "plain", label_size = 22) -> compare_loads
compare_loads
ggsave(compare_loads, file = "plots/compare_total_load_chick_adult_regions.png", width=14,height=18)

ggplot(subset(all_loads, method == "gerp" & loadtype == "exons"), aes(x = age, y = total_load)) + 
  geom_boxplot() + geom_point(position="jitter")
