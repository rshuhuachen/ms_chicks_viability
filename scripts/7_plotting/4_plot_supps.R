##### Supplementary descriptive figure ####
pacman::p_load(tidyverse, data.table)
source("scripts/theme_ggplot.R")

### a: froh distribution ####
load(file = "output/froh_chick_adult.RData")
froh$age <- gsub("chick", "Chick", froh$age)

#add longevity data 
load("data/phenotypic/phenotypes_lifetime.RData")
froh <- left_join(froh, pheno_wide[,c("id", "lifespan")], by = "id")
froh <- froh %>% mutate(age = as.factor(case_when(
  age == "Chick" ~ "Chick",
  lifespan == 1 ~ "Yearling",
  lifespan > 1 ~ "Adult"
)))

ggplot(froh, aes(x = froh, fill = age)) + geom_histogram() +
  labs(x = expression(italic(F)[ROH]), y = "Count", fill = "Age") +
  scale_fill_manual(values = c(clrs_hunting[2:4])) -> froh_dist

### b: snpeff distribution ####
# scaf info
scaf <- fread("data/metadata/30_largest.scafs.tsv")
scaf$scaf <- gsub(":", ";", scaf$scaf)
scaf$scaf <- gsub("\\.", "=", scaf$scaf)
scaf <- subset(scaf, scaf_no != 4)

vcf <- read.table("output/3_annotated_genome/ltet_filtered_ann_aa_chick_correct.vcf.gz")
vcf_29 <- subset(vcf, V1 %in% scaf$scaf)
vcf_29_nowarn <- subset(vcf_29, !grepl("WARNING", V8))

#subset regions
high <- subset(vcf_29_nowarn, grepl("HIGH", vcf_29_nowarn$V8))
mod <- subset(vcf_29_nowarn, grepl("MODERATE", vcf_29_nowarn$V8))
low <- subset(vcf_29_nowarn, grepl("LOW", vcf_29_nowarn$V8))
modifier <- subset(vcf_29_nowarn, grepl("MODIFIER", vcf_29_nowarn$V8))

sum <- data.frame(type = c("total_29scaf", "total_high", "total_mod", "total_low", "total_modifier",
                           "downstream_gene_variant","five_prime_UTR_premature_start_codon_gain_variant",
                           "five_prime_UTR_variant", "initiator_codon_variant",  
                           "intron_variant", "LOF", "missense_variant" ,"NMD","splice_acceptor_variant"   , 
                           "splice_donor_variant" ,  "splice_region_variant" , "start_lost" ,   
                           "stop_gained" , "stop_lost",
                           "stop_retained_variant" ,"synonymous_variant",
                           "three_prime_UTR_variant", "upstream_gene_variant", 
                           "intergenic_region"),
                  n_mutations = c(nrow(vcf_29_nowarn), nrow(high), nrow(mod), nrow(low), nrow(modifier),
                                  nrow(subset(vcf_29_nowarn, grepl("downstream_gene_variant", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("5_prime_UTR_premature_start_codon_gain_variant", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("5_prime_UTR_variant", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("initiator_codon_variant", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("intron_variant", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("LOF", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("missense_variant", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("NMD", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("splice_acceptor_variant", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("splice_donor_variant", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("splice_region_variant", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("start_lost", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("stop_gained", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("stop_lost", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("stop_retained_variant", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("synonymous_variant", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("3_prime_UTR_variant", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("upstream_gene_variant", vcf_29_nowarn$V8))),
                                  nrow(subset(vcf_29_nowarn, grepl("intergenic_region", vcf_29_nowarn$V8)))),
                  impact = c("NA", "NA", "NA", "NA", "NA", "modify","low", "modify", "low", "modify", "high", "moderate", "high", "high","high",    
                             "low","high", "high", "high","low", "low", "modify","modify", "low"))

write.csv(sum, file = "output/load/snpeff/n_mutations_per_type", quote=F, row.names = F)

