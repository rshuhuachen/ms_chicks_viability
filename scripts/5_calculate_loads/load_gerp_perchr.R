#### load packages ####
pacman::p_load(dplyr, data.table, tidyverse)

#### load gerps ####

load(file = "output/load/gerps_all.RData")

#### load functions #####
source("scripts/5_del_mutations/function_calculate_load.R")
source("scripts/theme_ggplot.R")

#### scaf name ####
load("/vol/cluster-data/rchen/git/scaffold_names_dovetail.RData")
gerp <- left_join(gerp, genome[,c("contig", "scaf_nr")], by = c("chr" = "contig"))

#### load per chr ####
list_gerps <- data.frame()
for (i in 1:length(unique(gerp$chr))){
  sub_chr <- subset(gerp, chr == unique(gerp$chr)[i])
  load_chr <- calculate_load_gerp(vcf = sub_chr, output_vcf = F, loadtype = unique(gerp$scaf_nr)[i])
  list_gerps <- rbind(list_gerps, load_chr)
  
}

#### get metadata ####
### load metadata
load("data/metadata/metadata_adult_chick.RData")

gerp_load <- left_join(list_gerps, meta, by = "id")
gerp_load <- gerp_load %>% mutate(age = as.factor(case_when(
  grepl("C", gerp_load$id) ~ "chick",
  grepl("D", gerp_load$id) ~ "adult"
)))


#### loop over chr to get output of model ####
n_mut <- unique(list_gerps[,c("loadtype", "n_total")])
summary_models <- data.frame()

for (i in 1:length(unique(gerp_load$loadtype))){
  
  ### total load
  model_total <- lm(total_load ~ age, data = subset(gerp_load, loadtype == unique(gerp_load$loadtype)[i]))
  sum_total <- summary(model_total)
  coef_total <- as.data.frame(sum_total$coefficients)
  
  summary_total <- data.frame(chr = unique(gerp_load$loadtype)[i],
                    est = coef_total$Estimate[2],
                    se = coef_total$`Std. Error`[2],
                     tval = coef_total$`t value`[2],
                    pval = coef_total$`Pr(>|t|)`[2],
                    load= "Total")
  
  ### hom load
  model_hom <- lm(hom_load ~ age, data = subset(gerp_load, loadtype == unique(gerp_load$loadtype)[i]))
  sum_hom <- summary(model_hom)
  coef_hom <- as.data.frame(sum_hom$coefficients)
  
  summary_hom <- data.frame(chr = unique(gerp_load$loadtype)[i],
                              est = coef_hom$Estimate[2],
                              se = coef_hom$`Std. Error`[2],
                              tval = coef_hom$`t value`[2],
                              pval = coef_hom$`Pr(>|t|)`[2],
                              load= "Homozygous")
  
  ### het load
  model_het <- lm(het_load ~ age, data = subset(gerp_load, loadtype == unique(gerp_load$loadtype)[i]))
  sum_het <- summary(model_het)
  coef_het <- as.data.frame(sum_het$coefficients)
  
  summary_het <- data.frame(chr = unique(gerp_load$loadtype)[i],
                              est = coef_het$Estimate[2],
                              se = coef_het$`Std. Error`[2],
                              tval = coef_het$`t value`[2],
                              pval = coef_het$`Pr(>|t|)`[2],
                              load= "Heterozygous")
  
  summary <- rbind(summary_total, summary_hom, summary_het)
  summary <- summary %>% mutate(sig = case_when(pval < 0.05 ~ "Significant", TRUE ~ "Non significant"))
  summary_models <- rbind(summary_models, summary)
}

summary_models <- left_join(summary_models, n_mut, by = c("chr" = "loadtype"))
save(summary_models, file = "output/load/summary_load_per_chr.RData")

#### plot ####
ggplot(summary_models, aes(x = est*1000, y = as.factor(chr), col = sig)) + 
  geom_point(size=4) + 
  geom_segment(aes(x = (est-se)*1000, xend = (est+se)*1000, yend = as.factor(chr)), col = "black", linewidth=1) +
  geom_vline(xintercept = 0, col = "darkred", linetype = "dotted") +
  facet_grid(~load, scales="free")+
  scale_color_manual(values=c(clr_grey, clr_gerp)) +
  labs(col = "Significant", x = expression("Beta estimate (chicks vs adults) * "~10^3), y = "Chromosome") -> gerp_per_chr

png("plots/load/gerp_per_chr_chicks_vs_adults.png", width=1000, height = 800)
gerp_per_chr
dev.off()


## nr of snps per chr ##
load("/vol/cluster-data/rchen/git/scaffold_names_dovetail.RData")

system('zcat data/processed/annotated/ltet_filtered_ann_aa_chick_correct.vcf.gz | grep -v "^#" | cut -f 1 | sort | uniq -c > output/n_snps_per_chr.txt')

n_snp <- read.table("output/n_snps_per_chr.txt")
names(n_snp) <- c("n_snp", "contig")

genome <- left_join(genome, n_snp, by = "contig")

ggplot(subset(genome[c(1:30),]), aes(x = n_bases, y = n_snp)) + geom_point() 
save(genome, file = "output/chr_mutations.RData")

## nr of gerp mutations per chr
n_gerp <- gerp %>% group_by(chr) %>% count()
genome <- left_join(genome, n_gerp, by = c("contig" = "chr"))

ggplot(subset(genome[c(1:30),]), aes(x = n_bases, y = n)) + geom_point() 

summary_models <- left_join(summary_models, genome, by = c("chr" = "scaf_nr"))

summary(lm(est ~ n, subset(summary_models, load == "Total"))) #nothing
summary(lm(est ~ n, subset(summary_models, load == "Homozygous"))) #nothing
summary(lm(est ~ n, subset(summary_models, load == "Heterozygous"))) #nothing

summary(lm(est ~ n_bases, subset(summary_models, load == "Total"))) #nothing
summary(lm(est ~ n_bases, subset(summary_models, load == "Homozygous"))) #nothing
summary(lm(est ~ n_bases, subset(summary_models, load == "Heterozygous"))) #nothing
