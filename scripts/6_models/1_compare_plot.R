### Load packages ####
pacman::p_load(tidyverse, data.table, cowplot, ggnewscale, lmerTest, lme4, DHARMa)
source("scripts/theme_ggplot.R")

### Load metadata ####
load("metadata/metadata_adult_chick.RData")

meta <- meta %>% mutate(age = case_when(
  grepl("C", meta$id) ~ "Chick",
  grepl("D", meta$id) ~ "Adult"
))

meta$age <- factor(meta$age, levels = c("Chick", "Adult"))

### Load phenotypes adults ####
load("metadata/phenotypes_wide_extra.RData")

### Load function posteriors ####
source("scripts/function_plot_posterior.R")

### Inbreeding ####
#### Load data, filter for autosome, quality filter ####
roh <- read.csv("output/inbreeding/bcftools_roh_chick_adult_clean.txt")

scaf <- fread("/vol/cluster-data/rchen/geneticload/gerp/analyses/git/cactus_insert_ltet_take4/data/genomes/30_largest.scafs.tsv")
scaf$scaf <- gsub(":", ";", scaf$scaf)
scaf$scaf <- gsub("\\.", "=", scaf$scaf)

roh_clean <- subset(roh, qual > 30 & length >= 100000 & nsnp >= 100 & chr != scaf$scaf[which(scaf$scaf_no==4)]) # same as load manuscript

roh_bp_perid <- roh_clean %>% group_by(id) %>% summarise(total_length_bp = sum(length))
roh_bp_perid$froh <- roh_bp_perid$total_length_bp/(1003452484-scaf$size_base[which(scaf$scaf_no==4)])
summary(roh_bp_perid$froh)
roh_bp_perid_n <- roh_clean %>% group_by(id) %>% count()
summary(roh_bp_perid_n$n)

froh <- left_join(roh_bp_perid, roh_bp_perid_n, by = "id")
froh <- froh %>% rename(n_rohs = n)
froh <- left_join(froh, meta, by = c("id"))

write.csv(froh, "output/inbreeding/froh_chick_adult_autosomal.txt", quote=F, row.names=F)

#### Model differences ####
ggplot(froh, aes(froh)) + geom_histogram() # normal enough, gamma doesn't improve fit
ggplot(froh, aes(site, froh)) + geom_boxplot() 
ggplot(froh, aes(year, froh)) + geom_boxplot() 
summary(lm(froh ~ site, data = froh))# leh lower
summary(lm(froh ~ year, data = froh))

null_froh  <- lmerTest::lmer(froh ~ 1 + (1|site) , data = froh)
model_froh <- lmerTest::lmer(froh ~ age + (1|site) , data = froh)
anova(null_froh, model_froh)
summary(model_froh)

simulation_gaussian <- simulateResiduals(model_froh, plot=T) # deviation from normality

#### Plot differences raw data ####

ggplot(froh, aes(x = age, y = froh)) + 
  geom_point(position="jitter", aes(col = age), size = 2, alpha = 0.6) +
  geom_boxplot(aes(fill = age), alpha = 0.5, outlier.shape=NA) + 
  theme(legend.position = "none") +
  geom_signif(comparisons = list(c("Chick", "Adult")), annotations="*",
              map_signif_level=TRUE, text = 10)+
  ylim(0.1, 0.4)+
  scale_color_manual(values=c(clrs_hunting[2], clrs_hunting[3]))+
  scale_fill_manual(values=c(clrs_hunting[2], clrs_hunting[3])) + 
  labs(x = "Age", y = expression(F[ROH])) -> plot_froh

ggsave(plot_froh, file = "plots/compare_froh_chick_adult.png", width = 10, height = 10)

#### Plot posteriors ####
load("results/bayes_models/age_froh_out.RData")
fit_froh <- fit
post_froh <- get_posterior_data_age(posteriors = fit_froh, response = "FROH", predictor = "ageAdult", name = "froh_chick_adult")

post_froh$diagnostics$autocor
post_froh$diagnostics$trace
post_froh$diagnostics$rhat
post_froh$diagnostics$neff

post_froh$plot +  scale_y_discrete(labels = c(expression(italic(F)[ROH])))

# ### FROH categories ####
# 
# #### Transform data ####
# 
# cats_kb <- c(100, 1000, 2000, 5000, 10000, 100000, 200000)
# cats_kb <- cats_kb*1000
# froh_per_cat <- data.frame()
# for (i in 2:length(cats_kb)){
#   roh_cat <- subset(roh_clean, qual>30 & length < cats_kb[i] & length >= cats_kb[i-1])
#   froh_cat <- roh_cat %>% group_by(id) %>% summarise(total_length_bp = sum(length))
#   froh_cat$froh <- froh_cat$total_length_bp/1011624072
#   froh_cat$cat <- paste0("froh_", as.character(cats_kb[i]))
#   froh_per_cat <- rbind(froh_per_cat, froh_cat)}
# 
# froh_per_cat <- froh_per_cat %>% select(-c(total_length_bp))
# 
# froh_per_cat <- froh_per_cat %>% mutate(pretty_cat = as.factor(case_when(
#   cat == "froh_1e+06" ~ "100 - 1000 kb",
#   cat == "froh_1e+07" ~ "5 - 10 Mb",
#   cat == "froh_1e+08" ~ "> 10 Mb",
#   cat == "froh_2e+06" ~ "1 - 2 Mb",
#   cat == "froh_5e+06" ~ "2 - 5 Mb")))
# 
# froh_per_cat$pretty_cat <- factor(froh_per_cat$pretty_cat, levels = c("100 - 1000 kb",
#                                                                       "1 - 2 Mb", "2 - 5 Mb","5 - 10 Mb","> 10 Mb"))
# froh_per_cat <- froh_per_cat %>% select(-c(cat))
# 
# froh_per_cat_wide <- spread(froh_per_cat, pretty_cat, froh)
# froh_per_cat_wide <- left_join(froh_per_cat_wide, meta, by = "id")
# froh_per_cat_wide$`100 - 1000 kb`[which(is.na(froh_per_cat_wide$`100 - 1000 kb`))] <- 0 #if NA in froh then 0
# froh_per_cat_wide$`1 - 2 Mb`[which(is.na(froh_per_cat_wide$`1 - 2 Mb`))] <- 0 #if NA in froh then 0
# froh_per_cat_wide$`5 - 10 Mb`[which(is.na(froh_per_cat_wide$`5 - 10 Mb`))] <- 0 #if NA in froh then 0
# froh_per_cat_wide$`> 10 Mb`[which(is.na(froh_per_cat_wide$`> 10 Mb`))] <- 0 #if NA in froh then 0
# froh_per_cat_wide$`2 - 5 Mb`[which(is.na(froh_per_cat_wide$`2 - 5 Mb`))] <- 0 #if NA in froh then 0
# 
# #### Model differences ####
# 
# null_froh_cat1  <- lmerTest::lmer(`100 - 1000 kb` ~ 1 + (1|site) , data = froh_per_cat_wide)
# model_froh_cat1 <- lmerTest::lmer(`100 - 1000 kb` ~ age + (1|site) , data = froh_per_cat_wide)
# anova(null_froh_cat1, model_froh_cat1) # sig, less in adults
# summary(model_froh_cat1)
# 
# null_froh_cat2  <- lmerTest::lmer(`1 - 2 Mb` ~ 1 + (1|site) , data = froh_per_cat_wide)
# model_froh_cat2 <- lmerTest::lmer(`1 - 2 Mb` ~ age + (1|site) , data = froh_per_cat_wide)
# anova(null_froh_cat2, model_froh_cat2) #ns
# summary(model_froh_cat2)
# 
# null_froh_cat3  <- lmerTest::lmer(`2 - 5 Mb` ~ 1 + (1|site) , data = froh_per_cat_wide)
# model_froh_cat3 <- lmerTest::lmer(`2 - 5 Mb` ~ age + (1|site) , data = froh_per_cat_wide)
# anova(null_froh_cat3, model_froh_cat3)
# summary(model_froh_cat3) #ns
# 
# null_froh_cat4  <- lmerTest::lmer(`5 - 10 Mb` ~ 1 + (1|site) , data = froh_per_cat_wide)
# model_froh_cat4 <- lmerTest::lmer(`5 - 10 Mb` ~ age + (1|site) , data = froh_per_cat_wide)
# anova(null_froh_cat4, model_froh_cat4)
# summary(model_froh_cat4) #ns
# 
# null_froh_cat5  <- lmerTest::lmer(`> 10 Mb` ~ 1 + (1|site) , data = froh_per_cat_wide)
# model_froh_cat5 <- lmerTest::lmer(`> 10 Mb` ~ age + (1|site) , data = froh_per_cat_wide)
# anova(null_froh_cat5, model_froh_cat5) 
# summary(model_froh_cat5) #ns

### GERP load ####
#### Load data ####
load(file="results/gerp_loads_chick_adult_no_exon.RData") #gerp_load

gerp_load <- gerp_load %>% mutate(
  gerp_load_hetero = gerp_count_p / n_genotyped,
  gerp_load_homo = gerp_count_r / n_genotyped,
  gerp_load_total_add = gerp_count_t_add / n_genotyped,
  gerp_load_total_codom = gerp_count_t_cum / n_genotyped)

gerp_load <- gerp_load %>% select(c(id, gerp_load_hetero:gerp_load_total_codom))
gerp_load <- left_join(gerp_load, meta, by = "id")

#### Model differences ####

### total additive
ggplot(gerp_load, aes(gerp_load_total_add)) + geom_histogram() # normal enough
ggplot(gerp_load, aes(site, gerp_load_total_add)) + geom_boxplot() 
ggplot(gerp_load, aes(year, gerp_load_total_add)) + geom_boxplot() 
summary(lm(gerp_load_total_add ~ site, data = gerp_load))# ns 
summary(lm(gerp_load_total_add ~ year, data = gerp_load))# ns

null_gerp_total_add  <- lmerTest::lmer(gerp_load_total_add ~ 1 + (1|site) , data = gerp_load)
model_gerp_total_add <- lmerTest::lmer(gerp_load_total_add ~ age + (1|site) , data = gerp_load)
anova(null_gerp_total_add, model_gerp_total_add) # sig
summary(model_gerp_total_add)

### homozygous
ggplot(gerp_load, aes(gerp_load_homo)) + geom_histogram() # normal enough
ggplot(gerp_load, aes(site, gerp_load_homo)) + geom_boxplot() 
ggplot(gerp_load, aes(year, gerp_load_homo)) + geom_boxplot() 
summary(lm(gerp_load_homo ~ site, data = gerp_load))# ns 
summary(lm(gerp_load_homo ~ year, data = gerp_load))# ns

null_gerp_homo  <- lmerTest::lmer(gerp_load_homo ~ 1 + (1|site) , data = gerp_load)
model_gerp_homo <- lmerTest::lmer(gerp_load_homo ~ age + (1|site) , data = gerp_load)
anova(null_gerp_homo, model_gerp_homo) # ns
summary(model_gerp_homo)

### heterozygous
ggplot(gerp_load, aes(gerp_load_hetero)) + geom_histogram() # normal enough
ggplot(gerp_load, aes(site, gerp_load_hetero)) + geom_boxplot() 
ggplot(gerp_load, aes(year, gerp_load_hetero)) + geom_boxplot() 
summary(lm(gerp_load_hetero ~ site, data = gerp_load))# ns 
summary(lm(gerp_load_hetero ~ year, data = gerp_load))# ns

null_gerp_hetero  <- lmerTest::lmer(gerp_load_hetero ~ 1 + (1|site) , data = gerp_load)
model_gerp_hetero <- lmerTest::lmer(gerp_load_hetero ~ age + (1|site) , data = gerp_load)
anova(null_gerp_hetero, model_gerp_hetero) # ns
summary(model_gerp_hetero)

#### Plot differences raw data ####

ggplot(gerp_load, aes(x = age, y = gerp_load_total_add)) + 
  geom_point(position="jitter", aes(col = age), size = 2, alpha = 0.6) +
  geom_boxplot(aes(fill = age), alpha = 0.5, outlier.shape=NA) + 
  theme(legend.position = "none") +
  geom_signif(comparisons = list(c("Chick", "Adult")), annotations="*",
              map_signif_level=TRUE, text = 10)+
  ylim(0.1445, 0.149)+
  scale_color_manual(values=c(clrs_hunting[2], clrs_hunting[3]))+
  scale_fill_manual(values=c(clrs_hunting[2], clrs_hunting[3])) + 
  labs(x = "Age", y = "Total GERP load") -> plot_gerp_total

ggplot(gerp_load, aes(x = age, y = gerp_load_homo)) + 
  geom_point(position="jitter", aes(col = age), size = 2, alpha = 0.6) +
  geom_boxplot(aes(fill = age), alpha = 0.5, outlier.shape=NA) + 
  theme(legend.position = "none") +
  scale_color_manual(values=c(clrs_hunting[2], clrs_hunting[3]))+
  scale_fill_manual(values=c(clrs_hunting[2], clrs_hunting[3])) + 
  labs(x = "Age", y = "Homozygous GERP load") -> plot_gerp_homo

ggplot(gerp_load, aes(x = age, y = gerp_load_hetero)) + 
  geom_point(position="jitter", aes(col = age), size = 2, alpha = 0.6) +
  geom_boxplot(aes(fill = age), alpha = 0.5, outlier.shape=NA) + 
  theme(legend.position = "none") +
  scale_color_manual(values=c(clrs_hunting[2], clrs_hunting[3]))+
  scale_fill_manual(values=c(clrs_hunting[2], clrs_hunting[3])) + 
  labs(x = "Age", y = "Heterozygous GERP load") -> plot_gerp_hetero

ggsave(plot_gerp_total, file = "plots/compare_gerp_total_chick_adult.png", width = 10, height = 10)
ggsave(plot_gerp_homo, file = "plots/compare_gerp_homo_chick_adult.png", width = 10, height = 10)
ggsave(plot_gerp_hetero, file = "plots/compare_gerp_hetero_chick_adult.png", width = 10, height = 10)

#### Plot posteriors ####

## total
load("results/bayes_models/age_gerp_load_total_add_out.RData")
fit_gerp_total <- fit
post_gerp_total <- get_posterior_data_age(posteriors = fit_gerp_total, response = "Total GERP load", predictor = "ageAdult", name = "gerp_total_chick_adult")

post_gerp_total$diagnostics$autocor
post_gerp_total$diagnostics$trace
post_gerp_total$diagnostics$rhat
post_gerp_total$diagnostics$neff

post_gerp_total$plot 

## homo
load("results/bayes_models/age_gerp_load_homo_out.RData")
fit_gerp_homo <- fit
post_gerp_homo <- get_posterior_data_age(posteriors = fit_gerp_homo, response = "Homozygous GERP load", predictor = "ageAdult", name = "gerp_homo_chick_adult")

post_gerp_homo$diagnostics$autocor
post_gerp_homo$diagnostics$trace
post_gerp_homo$diagnostics$rhat
post_gerp_homo$diagnostics$neff

post_gerp_homo$plot 

## hetero
load("results/bayes_models/age_gerp_load_hetero_out.RData")
fit_gerp_hetero <- fit
post_gerp_hetero <- get_posterior_data_age(posteriors = fit_gerp_hetero, response = "Heterozygous GERP load", predictor = "ageAdult", name = "gerp_hetero_chick_adult")

post_gerp_hetero$diagnostics$autocor
post_gerp_hetero$diagnostics$trace
post_gerp_hetero$diagnostics$rhat
post_gerp_hetero$diagnostics$neff

post_gerp_hetero$plot 

### combine in plot
intervals_gerp <- rbind(post_gerp_total$interval,
                        post_gerp_homo$interval,
                        post_gerp_hetero$interval)

intervals_gerp$response <- gsub(" GERP load", "", intervals_gerp$response)

intervals_gerp$response <- factor(intervals_gerp$response, levels = c("Heterozygous", "Homozygous", "Total"))

areas_gerp <- rbind(post_gerp_total$area,
                        post_gerp_homo$area,
                        post_gerp_hetero$area)
areas_gerp$response <- gsub(" GERP load", "", areas_gerp$response)
areas_gerp$response <- factor(areas_gerp$response, levels = c("Heterozygous", "Homozygous", "Total"))

brms_gerp <- split(areas_gerp, areas_gerp$interval)

brms_gerp$bottom <- brms_gerp$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

ggplot(data = brms_gerp$outer) +  
  aes(x = .data$x, y = .data$response) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = response, col = response))+
  geom_segment(data=intervals_gerp, aes(x = l, xend = h, yend = response), col = "black", linewidth=3)+
  geom_segment(data=intervals_gerp, aes(x = ll, xend = hh, yend = response), col = "black")+
  geom_point(data=intervals_gerp, aes(x = m, y = response), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Beta coefficient (sub) adult", title= "GERP", y = "Load type")+
  scale_fill_manual(values =alpha(c(clrs_hunting[1], clrs_hunting[2], clrs_hunting[3]), 0.7)) +
  scale_color_manual(values =c(clrs_hunting[1], clrs_hunting[2], clrs_hunting[3])) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> posterior_plot_gerp 

### SnpEff load ####
#### Load data ####
high <- fread("output/snpeff/genetic_load_da_nosex_30scaf_HIGH.tsv")
high_load <- left_join(high[,c("id", "load_p", "load_r", "load_t_add", "load_t_codom")], meta, by = "id")
names(high_load) <- c("id",  "high_load_hetero", "high_load_homo", "high_load_total_add", "high_load_total_codom", "site", "father", "year", "age")
write.csv(high_load, "output/high_load_chick_adult.txt", quote=F, row.names=F)

#### Model differences ####

### total additive
ggplot(high_load, aes(high_load_total_add)) + geom_histogram() # normal 
ggplot(high_load, aes(site, high_load_total_add)) + geom_boxplot() 
ggplot(high_load, aes(year, high_load_total_add)) + geom_boxplot() 
summary(lm(high_load_total_add ~ site, data = high_load))# ns 
summary(lm(high_load_total_add ~ year, data = high_load))# ns

null_high_total_add  <- lmerTest::lmer(high_load_total_add ~ 1 + (1|site) , data = high_load)
model_high_total_add <- lmerTest::lmer(high_load_total_add ~ age + (1|site) , data = high_load)
anova(null_high_total_add, model_high_total_add) # ns
summary(model_high_total_add)

### homozygous
ggplot(high_load, aes(high_load_homo)) + geom_histogram() # normal enough
ggplot(high_load, aes(site, high_load_homo)) + geom_boxplot() 
ggplot(high_load, aes(year, high_load_homo)) + geom_boxplot() 
summary(lm(high_load_homo ~ site, data = high_load))# ns 
summary(lm(high_load_homo ~ year, data = high_load))# sig

null_high_homo  <- lmerTest::lmer(high_load_homo ~ 1 + (1|year) + (1|site) , data = high_load)
model_high_homo <- lmerTest::lmer(high_load_homo ~ age + (1|year) + (1|site) , data = high_load)
anova(null_high_homo, model_high_homo) 
summary(model_high_homo)

### heterozygous
ggplot(high_load, aes(high_load_hetero)) + geom_histogram() # normal enough
ggplot(high_load, aes(site, high_load_hetero)) + geom_boxplot() 
ggplot(high_load, aes(year, high_load_hetero)) + geom_boxplot() 
summary(lm(high_load_hetero ~ site, data = high_load))# ns 
summary(lm(high_load_hetero ~ year, data = high_load))# sig

null_high_hetero  <- lmerTest::lmer(high_load_hetero ~ 1 + (1|year) + (1|site) , data = high_load)
model_high_hetero <- lmerTest::lmer(high_load_hetero ~ age + (1|year) + (1|site) , data = high_load)
anova(null_high_hetero, model_high_hetero) # sig
summary(model_high_hetero)

#### Plot differences raw data ####

ggplot(high_load, aes(x = age, y = high_load_total_add)) + 
  geom_point(position="jitter", aes(col = age), size = 2, alpha = 0.6) +
  geom_boxplot(aes(fill = age), alpha = 0.5, outlier.shape=NA) + 
  theme(legend.position = "none") +
  scale_color_manual(values=c(clrs_hunting[2], clrs_hunting[3]))+
  scale_fill_manual(values=c(clrs_hunting[2], clrs_hunting[3])) + 
  labs(x = "Age", y = "Total SnpEff load") -> plot_high_total

ggplot(high_load, aes(x = age, y = high_load_homo)) + 
  geom_point(position="jitter", aes(col = age), size = 2, alpha = 0.6) +
  geom_boxplot(aes(fill = age), alpha = 0.5, outlier.shape=NA) + 
  theme(legend.position = "none") +
  scale_color_manual(values=c(clrs_hunting[2], clrs_hunting[3]))+
  scale_fill_manual(values=c(clrs_hunting[2], clrs_hunting[3])) + 
  labs(x = "Age", y = "Homozygous SnpEff load") -> plot_high_homo

ggplot(high_load, aes(x = age, y = high_load_hetero)) + 
  geom_point(position="jitter", aes(col = age), size = 2, alpha = 0.6) +
  geom_boxplot(aes(fill = age), alpha = 0.5, outlier.shape=NA) + 
  geom_signif(comparisons = list(c("Chick", "Adult")), annotations="*",
              map_signif_level=TRUE, text = 10)+
  ylim(0.09, 0.16)+
  theme(legend.position = "none") +
  scale_color_manual(values=c(clrs_hunting[2], clrs_hunting[3]))+
  scale_fill_manual(values=c(clrs_hunting[2], clrs_hunting[3])) + 
  labs(x = "Age", y = "Heterozygous SnpEff load") -> plot_high_hetero

ggsave(plot_high_total, file = "plots/compare_high_total_chick_adult.png", width = 10, height = 10)
ggsave(plot_high_homo, file = "plots/compare_high_homo_chick_adult.png", width = 10, height = 10)
ggsave(plot_high_hetero, file = "plots/compare_high_hetero_chick_adult.png", width = 10, height = 10)

#### Plot posteriors ####

## total
load("results/bayes_models/age_high_load_total_add_out.RData")
fit_high_total <- fit
post_high_total <- get_posterior_data_age(posteriors = fit_high_total, response = "Total SnpEff load", predictor = "ageAdult", name = "high_total_chick_adult")

post_high_total$diagnostics$autocor
post_high_total$diagnostics$trace
post_high_total$diagnostics$rhat
post_high_total$diagnostics$neff

post_high_total$plot 

## homo
load("results/bayes_models/age_high_load_homo_out.RData")
fit_high_homo <- fit
post_high_homo <- get_posterior_data_age(posteriors = fit_high_homo, response = "Homozygous SnpEff load", predictor = "ageAdult", name = "high_homo_chick_adult")

post_high_homo$diagnostics$autocor
post_high_homo$diagnostics$trace
post_high_homo$diagnostics$rhat
post_high_homo$diagnostics$neff

post_high_homo$plot 

## hetero
load("results/bayes_models/age_high_load_hetero_out.RData")
fit_high_hetero <- fit
post_high_hetero <- get_posterior_data_age(posteriors = fit_high_hetero, response = "Heterozygous SnpEff load", predictor = "ageAdult", name = "high_hetero_chick_adult")

post_high_hetero$diagnostics$autocor
post_high_hetero$diagnostics$trace
post_high_hetero$diagnostics$rhat
post_high_hetero$diagnostics$neff

post_high_hetero$plot 

### combine in plot
intervals_high <- rbind(post_high_total$interval,
                        post_high_homo$interval,
                        post_high_hetero$interval)
intervals_high$response <- gsub(" SnpEff load", "", intervals_high$response)
intervals_high$response <- factor(intervals_high$response, levels = c("Heterozygous", "Homozygous", "Total"))

areas_high <- rbind(post_high_total$area,
                    post_high_homo$area,
                    post_high_hetero$area)
areas_high$response <- gsub(" SnpEff load", "", areas_high$response)
areas_high$response <- factor(areas_high$response, levels = c("Heterozygous", "Homozygous", "Total"))

brms_high <- split(areas_high, areas_high$interval)

brms_high$bottom <- brms_high$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

ggplot(data = brms_high$outer) +  
  aes(x = .data$x, y = .data$response) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = response, col = response))+
  geom_segment(data=intervals_high, aes(x = l, xend = h, yend = response), col = "black", linewidth=3)+
  geom_segment(data=intervals_high, aes(x = ll, xend = hh, yend = response), col = "black")+
  geom_point(data=intervals_high, aes(x = m, y = response), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Beta coefficient (sub) adult", title= "SnpEff", y = "Load type")+
  scale_fill_manual(values =alpha(c(clrs_hunting[1], clrs_hunting[2], clrs_hunting[3]), 0.7)) +
  scale_color_manual(values =c(clrs_hunting[1], clrs_hunting[2], clrs_hunting[3])) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> posterior_plot_high 


### Calculating a ratio compared to 'neutral' expectations including exons ####
#### Load data ####
gerps_data <- list.files(path="output/scores/", pattern = "load_per_id_cat", full.names=T)

gerps_cat <- data.frame()
for (i in 1:length(gerps_data)){
  scaf <- fread(gerps_data[i])
  gerps_cat <- rbind(gerps_cat, scaf)
}

sum_gerp_cat <- gerps_cat %>% select(-c(scafno)) %>% group_by(id, gerp_cat) %>%
  summarise(across(n_total:gerp_count_t_cum, sum))

sum_gerp_cat <- sum_gerp_cat %>% mutate(
  gerp_load_hetero = gerp_count_p / n_genotyped,
  gerp_load_homo = gerp_count_r / n_genotyped,
  gerp_load_total_add = gerp_count_t_add / n_genotyped,
  gerp_load_total_codom = gerp_count_t_cum / n_genotyped)

sum_gerp_cat <- as.data.frame(sum_gerp_cat)

# gerp_load <- left_join(subset(sum_gerp_cat[,c("id", "gerp_cat",  "gerp_load_hetero", "gerp_load_homo","gerp_load_total_add", "gerp_load_total_codom")], gerp_cat == "4-5"), 
#                        meta, by = "id")
# 
# write.csv(gerp_load, "output/gerp_load_chick_adult.txt", quote=F, row.names=F)

gerp_load_neutral <- left_join(subset(sum_gerp_cat[,c("id", "gerp_cat",  "gerp_load_hetero", "gerp_load_homo","gerp_load_total_add", 
                                                      "gerp_load_total_codom")], gerp_cat == "0"), meta, by = "id")

gerp_ratio <- left_join(gerp_load[,c("id", "age", "site", "year", "gerp_load_total_add")], gerp_load_neutral[,c("id", "age", "gerp_load_total_add")], 
                        suffix = c("_deleterious", "_neutral"), by = c("id", "age"))

gerp_ratio <- gerp_ratio %>% mutate(ratio = gerp_load_total_add_deleterious / gerp_load_total_add_neutral)

ggplot(gerp_ratio, aes(x = ratio)) + geom_histogram()

#### Model differences ####
null_gerp_age_ratio  <- lmerTest::lmer(ratio ~ 1 + (1|year) + (1|site) , data = gerp_ratio)
model_gerp_age_ratio <- lmerTest::lmer(ratio ~ age + (1|year) + (1|site) , data = gerp_ratio)
anova(null_gerp_age_ratio, model_gerp_age_ratio)
summary(model_gerp_age_ratio) #very sig

#### Plot differences raw data ####

ggplot(gerp_ratio, aes(x = age, y = ratio)) + 
  geom_point(position="jitter", aes(col = age), size = 2, alpha = 0.6) +
  geom_boxplot(aes(fill = age), alpha = 0.5, outlier.shape=NA) + 
  theme(legend.position = "none") +
  geom_signif(comparisons = list(c("Chick", "Adult")), annotations="*",
              map_signif_level=TRUE, text = 10)+
  ylim(0.775, 0.805)+
  scale_color_manual(values=c(clrs_hunting[2], clrs_hunting[3]))+
  scale_fill_manual(values=c(clrs_hunting[2], clrs_hunting[3])) + 
  labs(x = "Age", y = "Ratio total GERP ≥ 4 to total GERP < 0") -> plot_gerp_ratio

ggsave(plot_gerp_ratio, file = "plots/compare_gerp_ratio_total_chick_adult.png", width = 10, height = 10)

#### Plot posteriors ####
load("results/bayes_models/age_ratio_out.RData")
fit_ratio <- fit
post_ratio <- get_posterior_data_age(posteriors = fit_ratio, response = "Ratio GERP to neutral", predictor = "ageAdult", name = "ratio_chick_adult")

post_ratio$diagnostics$autocor
post_ratio$diagnostics$trace
post_ratio$diagnostics$rhat
post_ratio$diagnostics$neff

post_ratio$plot 

### Combine data ####

data <- left_join(meta, froh[c("id", "froh", "n_rohs")], by= "id") %>% 
  left_join(gerp_load[,c("id", "gerp_load_total_add", "gerp_load_total_codom", "gerp_load_homo", "gerp_load_hetero")], by = "id") %>%
  left_join(high_load[,c("id", "high_load_total_add", "high_load_total_codom", "high_load_homo", "high_load_hetero")], by = "id") %>%
  left_join(gerp_ratio[,c("id", "ratio")])

### Same comparisons but yearling vs adults lifespans ####

pheno <- pheno %>% mutate(lifespan_bi = as.factor(case_when(
  lifespan == 1 ~ "1 year",
  lifespan > 1 ~ "> 1 year"
)))
pheno$lifespan_bi <- factor(pheno$lifespan_bi, levels = c("1 year", "> 1 year"))

data_pheno <- left_join(data, pheno[,c("id", "lifespan_bi", "core")], by = "id")

##### Inbreeding ####
null_froh_lifespan  <- lmerTest::lmer(froh ~ 1 + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
model_froh_lifespan <- lmerTest::lmer(froh ~ lifespan_bi + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
anova(null_froh_lifespan, model_froh_lifespan)
summary(model_froh_lifespan) #ns

ggplot(subset(data_pheno, core == "core"), aes(x = lifespan_bi, y = froh)) + 
  geom_point(position="jitter", aes(col = lifespan_bi), size = 2, alpha = 0.6) +
  geom_boxplot(aes(fill = lifespan_bi), alpha = 0.5, outlier.shape=NA) + 
  theme(legend.position = "none") +
  scale_color_manual(values=c(clrs_hunting[4], clrs_hunting[1]))+
  scale_fill_manual(values=c(clrs_hunting[4], clrs_hunting[1])) + 
  labs(x = "Lifespan", y = expression(F[ROH])) -> plot_froh_lifespan

ggsave(plot_froh_lifespan, file = "plots/compare_froh_lifespan.png", width = 10, height = 10)

#### Plot posteriors ###
load("results/bayes_models/lifespan_froh_out.RData")
fit_froh_ls <- fit
post_froh_ls <- get_posterior_data_lifespan(posteriors = fit_froh_ls, response = "FROH", predictor = "lifespan_bi>1year", name = "froh_lifespan")

post_froh_ls$diagnostics$autocor
post_froh_ls$diagnostics$trace
post_froh_ls$diagnostics$rhat
post_froh_ls$diagnostics$neff

post_froh_ls$plot +  scale_y_discrete(labels = c(expression(italic(F)[ROH])))

##### GERP total ####

null_gerp_lifespan  <- lmerTest::lmer(gerp_load_total_add ~ 1 + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
model_gerp_lifespan <- lmerTest::lmer(gerp_load_total_add ~ lifespan_bi + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
anova(null_gerp_lifespan, model_gerp_lifespan)
summary(model_gerp_lifespan) #ns

ggplot(subset(data_pheno, core == "core"), aes(x = lifespan_bi, y = gerp_load_total_add)) + 
  geom_point(position="jitter", aes(col = lifespan_bi), size = 2, alpha = 0.6) +
  geom_boxplot(aes(fill = lifespan_bi), alpha = 0.5, outlier.shape=NA) + 
  theme(legend.position = "none") +
  scale_color_manual(values=c(clrs_hunting[4], clrs_hunting[1]))+
  scale_fill_manual(values=c(clrs_hunting[4], clrs_hunting[1])) + 
  labs(x = "Lifespan", y = "Total GERP load") -> plot_gerp_total_lifespan

ggsave(plot_gerp_total_lifespan, file = "plots/compare_gerp_total_lifespan.png", width = 10, height = 10)

##### GERP homozygosity ####
null_gerp_lifespan_homo  <- lmerTest::lmer(gerp_load_homo ~ 1 + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
model_gerp_lifespan_homo <- lmerTest::lmer(gerp_load_homo ~ lifespan_bi + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
anova(null_gerp_lifespan_homo, model_gerp_lifespan_homo)
summary(model_gerp_lifespan_homo) #ns

##### GERP heterozygosity ####
null_gerp_lifespan_hetero  <- lmerTest::lmer(gerp_load_hetero ~ 1 + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
model_gerp_lifespan_hetero <- lmerTest::lmer(gerp_load_hetero ~ lifespan_bi + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
anova(null_gerp_lifespan_hetero, model_gerp_lifespan_hetero)
summary(model_gerp_lifespan_hetero) #ns

#### Plot posteriors ####

## total
load("results/bayes_models/lifespan_gerp_load_total_add_out.RData")
fit_gerp_total_ls <- fit
post_gerp_total_ls <- get_posterior_data_lifespan(posteriors = fit_gerp_total_ls, response = "Total GERP load", predictor = "lifespan_bi>1year", 
                                          name = "gerp_total_lifespan")

post_gerp_total_ls$diagnostics$autocor
post_gerp_total_ls$diagnostics$trace
post_gerp_total_ls$diagnostics$rhat
post_gerp_total_ls$diagnostics$neff

post_gerp_total_ls$plot 

## homo
load("results/bayes_models/lifespan_gerp_load_homo_out.RData")
fit_gerp_homo_ls <- fit
post_gerp_homo_ls <- get_posterior_data_lifespan(posteriors = fit_gerp_homo_ls, response = "Homozygous GERP load", predictor = "lifespan_bi>1year", 
                                         name = "gerp_homo_lifespan")

post_gerp_homo_ls$diagnostics$autocor
post_gerp_homo_ls$diagnostics$trace
post_gerp_homo_ls$diagnostics$rhat
post_gerp_homo_ls$diagnostics$neff

post_gerp_homo_ls$plot 

## hetero
load("results/bayes_models/lifespan_gerp_load_hetero_out.RData")
fit_gerp_hetero_ls <- fit
post_gerp_hetero_ls <- get_posterior_data_lifespan(posteriors = fit_gerp_hetero_ls, response = "Heterozygous GERP load", predictor = "lifespan_bi>1year", 
                                           name = "gerp_hetero_lifespan")

post_gerp_hetero_ls$diagnostics$autocor
post_gerp_hetero_ls$diagnostics$trace
post_gerp_hetero_ls$diagnostics$rhat
post_gerp_hetero_ls$diagnostics$neff

post_gerp_hetero_ls$plot 

### combine in plot
intervals_gerp_ls <- rbind(post_gerp_total_ls$interval,
                        post_gerp_homo_ls$interval,
                        post_gerp_hetero_ls$interval)

intervals_gerp_ls$response <- gsub(" GERP load", "", intervals_gerp_ls$response)

intervals_gerp_ls$response <- factor(intervals_gerp_ls$response, levels = c("Heterozygous", "Homozygous", "Total"))

areas_gerp_ls <- rbind(post_gerp_total_ls$area,
                    post_gerp_homo_ls$area,
                    post_gerp_hetero_ls$area)
areas_gerp_ls$response <- gsub(" GERP load", "", areas_gerp_ls$response)
areas_gerp_ls$response <- factor(areas_gerp_ls$response, levels = c("Heterozygous", "Homozygous", "Total"))

brms_gerp_ls <- split(areas_gerp_ls, areas_gerp_ls$interval)

brms_gerp_ls$bottom <- brms_gerp_ls$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

ggplot(data = brms_gerp_ls$outer) +  
  aes(x = .data$x, y = .data$response) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = response, col = response))+
  geom_segment(data=intervals_gerp_ls, aes(x = l, xend = h, yend = response), col = "black", linewidth=3)+
  geom_segment(data=intervals_gerp_ls, aes(x = ll, xend = hh, yend = response), col = "black")+
  geom_point(data=intervals_gerp_ls, aes(x = m, y = response), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Beta coefficient lifespan > 1 year", title= "GERP", y = "Load type")+
  scale_fill_manual(values =alpha(c(clrs_hunting[1], clrs_hunting[2], clrs_hunting[3]), 0.7)) +
  scale_color_manual(values =c(clrs_hunting[1], clrs_hunting[2], clrs_hunting[3])) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> posterior_plot_gerp 

##### GERP ratio ####
null_gerp_lifespan_ratio  <- lmerTest::lmer(ratio ~ 1 + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
model_gerp_lifespan_ratio <- lmerTest::lmer(ratio ~ lifespan_bi + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
anova(null_gerp_lifespan_ratio, model_gerp_lifespan_ratio)
summary(model_gerp_lifespan_ratio) #ns

ggplot(subset(data_pheno, core == "core"), aes(x = lifespan_bi, y = ratio)) + 
  geom_point(position="jitter", aes(col = lifespan_bi), size = 2, alpha = 0.6) +
  geom_boxplot(aes(fill = lifespan_bi), alpha = 0.5, outlier.shape=NA) + 
  theme(legend.position = "none") +
  scale_color_manual(values=c(clrs_hunting[4], clrs_hunting[1]))+
  scale_fill_manual(values=c(clrs_hunting[4], clrs_hunting[1])) + 
  labs(x = "Lifespan", y = "Ratio total GERP ≥ 4 to total GERP < 0") -> plot_gerp_ratio_lifespan

ggsave(plot_gerp_ratio_lifespan, file = "plots/compare_gerp_ratio_lifespan.png", width = 10, height = 10)

### posteriors
load("results/bayes_models/lifespan_ratio_out.RData")
fit_ratio_ls <- fit
post_ratio_ls <- get_posterior_data_lifespan(posteriors = fit_ratio_ls, response = "Ratio GERP to neutral", predictor = "lifespan_bi>1year", name = "ratio_lifespan")

post_ratio_ls$diagnostics$autocor
post_ratio_ls$diagnostics$trace
post_ratio_ls$diagnostics$rhat
post_ratio_ls$diagnostics$neff

post_ratio_ls$plot 

##### SnpEff total ####
null_high_lifespan  <- lmerTest::lmer(high_load_total_add ~ 1 + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
model_high_lifespan <- lmerTest::lmer(high_load_total_add ~ lifespan_bi + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
anova(null_high_lifespan, model_high_lifespan)
summary(model_high_lifespan) #ns

ggplot(subset(data_pheno, core == "core"), aes(x = lifespan_bi, y = high_load_total_add)) + 
  geom_point(position="jitter", aes(col = lifespan_bi), size = 2, alpha = 0.6) +
  geom_boxplot(aes(fill = lifespan_bi), alpha = 0.5, outlier.shape=NA) + 
  theme(legend.position = "none") +
  scale_color_manual(values=c(clrs_hunting[4], clrs_hunting[1]))+
  scale_fill_manual(values=c(clrs_hunting[4], clrs_hunting[1])) + 
  labs(x = "Lifespan", y = "Total SnpEff load") -> plot_high_total_lifespan

ggsave(plot_high_total_lifespan, file = "plots/compare_high_total_lifespan.png", width = 10, height = 10)

##### SnpEff homozygosity ####
null_high_lifespan_homo  <- lmerTest::lmer(high_load_homo ~ 1 + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
model_high_lifespan_homo <- lmerTest::lmer(high_load_homo ~ lifespan_bi + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
anova(null_high_lifespan_homo, model_high_lifespan_homo)
summary(model_high_lifespan_homo) #ns

##### SnpEff heterozygosity ####
null_high_lifespan_hetero  <- lmerTest::lmer(high_load_hetero ~ 1 + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
model_high_lifespan_hetero <- lmerTest::lmer(high_load_hetero ~ lifespan_bi + (1|year) + (1|site) , data = subset(data_pheno, core == "core"))
anova(null_high_lifespan_hetero, model_high_lifespan_hetero)
summary(model_high_lifespan_hetero) #ns

#### Plot posteriors ####

## total
load("results/bayes_models/lifespan_high_load_total_add_out.RData")
fit_high_total_ls <- fit
post_high_total_ls <- get_posterior_data_lifespan(posteriors = fit_high_total_ls, response = "Total SnpEff load", predictor = "lifespan_bi>1year", 
                                          name = "high_total_lifespan")

post_high_total_ls$diagnostics$autocor
post_high_total_ls$diagnostics$trace
post_high_total_ls$diagnostics$rhat
post_high_total_ls$diagnostics$neff

post_high_total_ls$plot 

## homo
load("results/bayes_models/lifespan_high_load_homo_out.RData")
fit_high_homo_ls <- fit
post_high_homo_ls <- get_posterior_data_lifespan(posteriors = fit_high_homo_ls, response = "Homozygous SnpEff load", predictor = "lifespan_bi>1year", 
                                         name = "high_homo_lifespan")

post_high_homo_ls$diagnostics$autocor
post_high_homo_ls$diagnostics$trace
post_high_homo_ls$diagnostics$rhat
post_high_homo_ls$diagnostics$neff

post_high_homo_ls$plot 

## hetero
load("results/bayes_models/lifespan_high_load_hetero_out.RData")
fit_high_hetero_ls <- fit
post_high_hetero_ls <- get_posterior_data_lifespan(posteriors = fit_high_hetero_ls, response = "Heterozygous SnpEff load", predictor = "lifespan_bi>1year", 
                                           name = "high_hetero_lifespan")

post_high_hetero_ls$diagnostics$autocor
post_high_hetero_ls$diagnostics$trace
post_high_hetero_ls$diagnostics$rhat
post_high_hetero_ls$diagnostics$neff

post_high_hetero$plot 

### combine in plot
intervals_high_ls <- rbind(post_high_total_ls$interval,
                        post_high_homo_ls$interval,
                        post_high_hetero_ls$interval)
intervals_high_ls$response <- gsub(" SnpEff load", "", intervals_high_ls$response)
intervals_high_ls$response <- factor(intervals_high_ls$response, levels = c("Heterozygous", "Homozygous", "Total"))

areas_high_ls <- rbind(post_high_total_ls$area,
                    post_high_homo_ls$area,
                    post_high_hetero_ls$area)
areas_high_ls$response <- gsub(" SnpEff load", "", areas_high_ls$response)
areas_high_ls$response <- factor(areas_high_ls$response, levels = c("Heterozygous", "Homozygous", "Total"))

brms_high_ls <- split(areas_high_ls, areas_high_ls$interval)

brms_high_ls$bottom <- brms_high_ls$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

ggplot(data = brms_high_ls$outer) +  
  aes(x = .data$x, y = .data$response) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = response, col = response))+
  geom_segment(data=intervals_high_ls, aes(x = l, xend = h, yend = response), col = "black", linewidth=3)+
  geom_segment(data=intervals_high_ls, aes(x = ll, xend = hh, yend = response), col = "black")+
  geom_point(data=intervals_high_ls, aes(x = m, y = response), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Beta coefficient lifespan > 1 year", title= "SnpEff", y = "Load type")+
  scale_fill_manual(values =alpha(c(clrs_hunting[1], clrs_hunting[2], clrs_hunting[3]), 0.7)) +
  scale_color_manual(values =c(clrs_hunting[1], clrs_hunting[2], clrs_hunting[3])) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> posterior_plot_high_ls 


### Save data ####
save(data, file = "results/froh_loads_chick_adult.RData")

### Purging plot GERP ####
# make df long
longdata <- gather(data, loadtype, load, gerp_load_total_add:high_load_hetero, factor_key=T)

# split up load type and impact cat
longdata <- longdata %>% mutate(impact = as.factor(case_when(
  grepl("high", longdata$loadtype) ~ "SnpEff",
  grepl("gerp", longdata$loadtype) ~ "GERP")))

longdata$loadtype <- gsub("high_", "", longdata$loadtype)
longdata$loadtype <- gsub("gerp_", "", longdata$loadtype)

longdata$loadtype <- gsub("load_hetero", "Heterozygous", longdata$loadtype)
longdata$loadtype <- gsub("load_homo", "Homozygous", longdata$loadtype)
longdata$loadtype <- gsub("load_total_add", "Total (weighted)", longdata$loadtype)
longdata$loadtype <- gsub("load_total_codom", "Total (unweighted)", longdata$loadtype)

### Weighted (additive)
# define function to calculate hline without purging

lm_adult_high_w <- lm(load ~ froh, data = subset(longdata, impact=="SnpEff" & loadtype=="Total (weighted)" & age == "Adult"))
lm_chick_high_w <- lm(load ~ froh, data = subset(longdata, impact=="SnpEff" & loadtype=="Total (weighted)" & age == "Chick"))
lm_adult_gerp_w <- lm(load ~ froh, data = subset(longdata, impact=="GERP" & loadtype=="Total (weighted)" & age == "Adult"))
lm_chick_gerp_w <- lm(load ~ froh, data = subset(longdata, impact=="GERP" & loadtype=="Total (weighted)" & age == "Chick"))

calc_nopurge <- function(froh, model){
  return(as.numeric(model$coefficients[1]) + froh*as.numeric(model$coefficients[2]))
}

nopurge_w <- data.frame(impact= c("GERP", "GERP", "SnpEff", "SnpEff"),
                        age = c("Adult", "Chick", "Adult", "Chick"),
                        yline = c(calc_nopurge(froh = min(longdata$froh[which(longdata$age=="Adult")], na.rm=T), model = lm_adult_gerp_w),
                                  calc_nopurge(froh = min(longdata$froh[which(longdata$age=="Chick")], na.rm=T), model = lm_chick_gerp_w),
                                  calc_nopurge(froh = min(longdata$froh[which(longdata$age=="Adult")], na.rm=T), model = lm_adult_high_w),
                                  calc_nopurge(froh = min(longdata$froh[which(longdata$age=="Chick")], na.rm=T), model = lm_chick_high_w)))



ggplot(subset(longdata, loadtype != "Total (unweighted)"), 
       aes(x = froh, y = load, col = loadtype)) + 
  geom_point() + geom_smooth(method="lm", fullrange = T) +
  facet_grid(impact~factor(age, levels = c("Chick", "Adult")), scales="free_y") +
  geom_hline(data=nopurge_w, aes(yintercept=yline), linetype = "dotted", col = "grey30")+
  scale_color_manual(values=c(clrs_hunting[1],
                              clrs_hunting[3],
                              clrs_hunting[5]))+
  theme(legend.position = "top")+
  labs(x = expression(F[ROH]), y = "Load", col = "Load type") -> froh_load_total_weighted

ggsave(froh_load_total_weighted, file="plots/froh_load_adult_chick_weighted.png", width=12, height=10)

### Codominant (unweighted)

lm_adult_high_uw <- lm(load ~ froh, data = subset(longdata, impact=="SnpEff" & loadtype=="Total (unweighted)" & age == "Adult"))
lm_chick_high_uw <- lm(load ~ froh, data = subset(longdata, impact=="SnpEff" & loadtype=="Total (unweighted)" & age == "Chick"))
lm_adult_gerp_uw <- lm(load ~ froh, data = subset(longdata, impact=="GERP" & loadtype=="Total (unweighted)" & age == "Adult"))
lm_chick_gerp_uw <- lm(load ~ froh, data = subset(longdata, impact=="GERP" & loadtype=="Total (unweighted)" & age == "Chick"))

nopurge_uw <- data.frame(impact= c("GERP", "GERP", "SnpEff", "SnpEff"),
                         age = c("Adult", "Chick", "Adult", "Chick"),
                         yline = c(calc_nopurge(froh = min(longdata$froh[which(longdata$age=="Adult")], na.rm=T), model = lm_adult_gerp_uw),
                                   calc_nopurge(froh = min(longdata$froh[which(longdata$age=="Chick")], na.rm=T), model = lm_chick_gerp_uw),
                                   calc_nopurge(froh = min(longdata$froh[which(longdata$age=="Adult")], na.rm=T), model = lm_adult_high_uw),
                                   calc_nopurge(froh = min(longdata$froh[which(longdata$age=="Chick")], na.rm=T), model = lm_chick_high_uw)))


ggplot(subset(longdata, loadtype != "Total (weighted)"), 
       aes(x = froh, y = load, col = loadtype)) + 
  geom_point() + geom_smooth(method="lm", fullrange=T) +
  facet_grid(impact~factor(age, levels = c("Chick", "Adult")), scales="free_y") +
  geom_hline(data=nopurge_uw, aes(yintercept=yline), linetype = "dotted", col = "grey30")+
  scale_color_manual(values=c(clrs_hunting[1],
                              clrs_hunting[3],
                              clrs_hunting[5]))+
  theme(legend.position = "top")+
  labs(x = expression(F[ROH]), y = "Load", col = "Load type") -> froh_load_total_unweighted

ggsave(froh_load_total_unweighted, file="plots/froh_load_adult_chick_unweighted.png", width=12, height=10)

cowplot::plot_grid(froh_load_total_weighted, froh_load_total_unweighted, ncol = 1, align = "hv", axis = "lb",
                   labels = c("a) Weighted", "b) Unweighted"), label_fontface = "plain", label_size = 22) -> froh_load_total

ggsave(froh_load_total, file="plots/froh_load_adult_chick.png", width=12, height=18)

### Combine plots ####

# raw data totals
cowplot::plot_grid(plot_froh, plot_froh_lifespan,
                   plot_gerp_total, plot_gerp_total_lifespan, 
                   plot_high_total, plot_high_total_lifespan,
                   plot_gerp_ratio, plot_gerp_ratio_lifespan, ncol = 2, align = "hv", axis = "lb",
                   labels = "auto", label_fontface = "plain", label_size = 22) -> compare_total

ggsave(compare_total, file="plots/total_compare_age_lifespan.png", width=14, height=22)

# posteriors (sub)adult compared to chick
cowplot::plot_grid(post_froh$plot+  scale_y_discrete(labels = c(expression(italic(F)[ROH]))) + 
                     theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), color = "black")), 
                   posterior_plot_gerp,
                   posterior_plot_high, 
                   rel_heights = c(0.5, 1, 1, 0.5),
                   post_ratio$plot + theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), color = "black")), ncol = 1, align = "hv", axis = "lb",
                   labels = "auto", label_fontface = "plain", label_size = 22) -> chick_adult_posteriors

ggsave(chick_adult_posteriors, file="plots/posteriors_chick_adult.png", width=12, height=20)

# posteriors lifespan
cowplot::plot_grid(post_froh_ls$plot+  scale_y_discrete(labels = c(expression(italic(F)[ROH]))) + 
                     theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), color = "black")), 
                   posterior_plot_gerp_ls,
                   posterior_plot_high_ls, 
                   rel_heights = c(0.5, 1, 1, 0.5),
                   post_ratio_ls$plot + theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0), color = "black")), ncol = 1, align = "hv", axis = "lb",
                   labels = "auto", label_fontface = "plain", label_size = 22) -> lifespan_posteriors

ggsave(lifespan_posteriors, file="plots/posteriors_lifespan.png", width=12, height=20)
