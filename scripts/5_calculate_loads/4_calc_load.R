pacman::p_load(dplyr, data.table, tidyverse)
### load names ####
load(file = "output/nam")

#### high impact snpeff mutations ####
high <- read.table("output/4_load/snpeff/ltet_filtered_ann_aa_chick_high.vcf.gz")

## load headers
names <- fread("output/2_inbreeding/chicks_adults_samples.txt", header = F)
names(names) <- "file_id"

## load id info
# add real id
ids <- read.csv("data/metadata/file_list_all_bgi_clean.csv")
ids$file_id <- gsub(".sorted.bam", "", ids$file)
names <- left_join(names, ids[,c("file_id", "id")], by = c("file_id"))
names <- names %>% mutate(id = case_when(
  is.na(id) ~ file_id,
  TRUE ~ id
))

names(high)[1:9] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT")
names(high)[c(10:ncol(high))] <- names$id

source("scripts/5_calculate_loads/function_calculate_load.R")
source("scripts/theme_ggplot.R")

## calculate load
high_load <- calculate_load_snpeff(vcf = high, loadtype = "high", output_vcf = F)

## load metadata
load("data/metadata/metadata_adult_chick.RData")

high_load <- left_join(high_load, meta, by = "id")
high_load <- high_load %>% mutate(age = as.factor(case_when(
  grepl("C", high_load$id) ~ "chick",
  grepl("D", high_load$id) ~ "adult"
)))

## plot raw data 
compare_high <- ggplot(high_load, aes(x = age, y = total_load)) + 
  geom_boxplot(aes(fill = age), outlier.shape = NA) + 
  geom_point(position = "jitter", aes(col = age), size = 2) + labs(x = "Age", y = "Total SnpEff load")+
  scale_color_manual(values = c(clr_high, clr_highlight)) +
  scale_fill_manual(values = alpha(c(clr_high, clr_highlight), 0.4)) +
  theme(legend.position = "none")

compare_high

png("plots/load/compare_high_total.png", width=600, height = 600)
compare_high
dev.off()

summary(lmerTest::lmer(total_load ~ age + (1|site), data = high_load)) # ns

compare_high_hom <- ggplot(high_load, aes(x = age, y = hom_load)) + 
  geom_boxplot(aes(fill = age), outlier.shape = NA) + 
  geom_point(position = "jitter", aes(col = age), size = 2) + labs(x = "Age", y = "Homozygous SnpEff load")+
  scale_color_manual(values = c(clr_high, clr_highlight)) +
  scale_fill_manual(values = alpha(c(clr_high, clr_highlight), 0.4)) +
  theme(legend.position = "none")

compare_high_hom

png("plots/load/compare_high_hom.png", width=600, height = 600)
compare_high_hom
dev.off()


### GERP ####

## load data
load(file = "output/4_load/gerp/gerps_all.RData")

## calculate load
gerp_load <- calculate_load_gerp(vcf = gerp, output_vcf = F, loadtype = "gerp")

## load metadata
load("data/metadata/metadata_adult_chick.RData")

gerp_load <- left_join(gerp_load, meta, by = "id")
gerp_load <- gerp_load %>% mutate(age = as.factor(case_when(
  grepl("C", gerp_load$id) ~ "chick",
  grepl("D", gerp_load$id) ~ "adult"
)))

## plot raw data
compare_gerp <- ggplot(gerp_load, aes(x = age, y = total_load)) + 
  geom_boxplot(aes(fill = age), outlier.shape = NA) + 
  geom_point(position = "jitter", aes(col = age), size = 2) + labs(x = "Age", y = "Total GERP load")+
  scale_color_manual(values = c(clr_gerp, clr_highlight)) +
  scale_fill_manual(values = alpha(c(clr_gerp, clr_highlight), 0.4)) +
  theme(legend.position = "none")

compare_gerp

png("plots/load/compare_gerp_total.png", width=600, height = 600)
compare_gerp
dev.off()

summary(lmerTest::lmer(total_load ~ age + (1|site), data = gerp_load)) # sig
summary(lmerTest::lmer(hom_load ~ age + (1|site), data = gerp_load)) # almost sig
summary(lmerTest::lmer(het_load ~ age + (1|site), data = gerp_load)) # non sig

compare_gerp_hom <- ggplot(gerp_load, aes(x = age, y = hom_load)) + 
  geom_boxplot(aes(fill = age), outlier.shape = NA) + 
  geom_point(position = "jitter", aes(col = age), size = 2) + labs(x = "Age", y = "Homozygous GERP load")+
  scale_color_manual(values = c(clr_gerp, clr_highlight)) +
  scale_fill_manual(values = alpha(c(clr_gerp, clr_highlight), 0.4)) +
  theme(legend.position = "none")

compare_gerp_hom

png("plots/load/compare_gerp_hom.png", width=600, height = 600)
compare_gerp_hom
dev.off()

### Save loads ####
loads <- rbind(high_load, gerp_load)
save(loads, file = "output/4_load/loads.RData")


