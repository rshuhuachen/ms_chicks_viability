### load packages ####
pacman::p_load(data.table, tidyverse)
source("scripts/theme_ggplot.R")

### run BCFtools
### Set the directories of both your pre-processed VCF file as well as the absolute path of the github repository 

VCF= paste0(getwd(), "/output/genomic/genomes/processed/ltet_snps_filtered_chicks_adults.vcf.gz")
ROH_OUT = paste0(getwd(), "/output/2_inbreeding/bcftools_roh_chicks_adult.txt")

# run bcftools
system(paste0("bcftools roh -G30 --AF-dflt 0.4 ", VCF, " -o ", ROH_OUT))

# filter bcftools to only include the ROHs
system(paste0("grep \"^RG\" ", ROH_OUT, "> ", getwd(), "/output/2_inbreeding/bcftools_roh_chick_adult_rg.txt")) 

roh <- fread("output/2_inbreeding/bcftools_roh_chick_adult_rg.txt")
names(roh) <- c("state", "file_id", "chr", "start", "end", "length", "nsnp", "qual")

#summarise
summary(roh$qual)

#add real id
ids <- read.csv("data/metadata/file_list_all_bgi_clean.csv")
ids$file <- gsub(".sorted.bam", "", ids$file)

roh <- left_join(roh, ids[,c("file", "id")], by = c("file_id" = "file"))

roh <- roh %>% mutate(id = case_when(
  is.na(id) ~ file_id,
  TRUE ~ id
))

roh <- roh %>% mutate(age = as.factor(case_when(
  grepl("C", roh$id) ~ "chick",
  grepl("D", roh$id) ~ "adult"
)))

write.csv(roh, "output/2_inbreeding/bcftools_roh_chick_adult.txt", quote=F, row.names=F)

#### Load in clean and filter sex scaf ####
roh <- read.csv("output/2_inbreeding/bcftools_roh_chick_adult.txt")

scaf <- fread("data/metadata/30_largest.scafs.tsv")
scaf$scaf <- gsub(":", ";", scaf$scaf)
scaf$scaf <- gsub("\\.", "=", scaf$scaf)

roh <- subset(roh, qual>30 & chr != scaf$scaf[which(scaf$scaf_no==4)]) #filter out sex scaf 4
roh_bp_perid_raw <- roh %>% group_by(id) %>% summarise(total_length_bp = sum(length))
roh_bp_perid_raw$froh <- roh_bp_perid_raw$total_length_bp/(1003452484-scaf$size_base[which(scaf$scaf_no==4)])
summary(roh_bp_perid_raw$froh)

roh_bp_perid_n_raw <- roh %>% group_by(id) %>% count()
summary(roh_bp_perid_n_raw$n)

froh_raw <- left_join(roh_bp_perid_raw, unique(roh[,c("id", "age")]), by = "id")
ggplot(froh_raw, aes(x = age, y = froh))+ geom_boxplot() + theme_classic()
summary(lm(froh ~ age, data = froh_raw)) #ns

# filtering 
roh_clean <- subset(roh, qual > 30 & length >= 100000 & nsnp >= 100 & chr != scaf$scaf[which(scaf$scaf_no==4)]) 

roh_bp_perid <- roh_clean %>% group_by(id) %>% summarise(total_length_bp = sum(length))
roh_bp_perid$froh <- roh_bp_perid$total_length_bp/(1003452484-scaf$size_base[which(scaf$scaf_no==4)])
summary(roh_bp_perid$froh)
roh_bp_perid_n <- roh_clean %>% group_by(id) %>% count()
summary(roh_bp_perid_n$n)

froh <- left_join(roh_bp_perid, unique(roh[,c("id", "age")]), by = "id")
ggplot(froh, aes(x = age, y = froh))+ geom_boxplot() 

compare_froh <- ggplot(froh, aes(x = age, y = froh)) + 
  geom_boxplot(aes(fill = age), outlier.shape = NA) + 
  geom_point(position = "jitter", aes(col = age), size = 2) + labs(x = "Age", y = expression(F[ROH]))+
  scale_color_manual(values = c(clr_froh, clr_highlight)) +
  scale_fill_manual(values = alpha(c(clr_froh, clr_highlight), 0.4)) +
  theme(legend.position = "none")

compare_froh

png("plots/inbreeding/compare_froh.png", width=600, height = 600)
compare_froh
dev.off()

froh$age <- factor(froh$age, levels = c("chick", "adult"))
summary(lm(froh ~ age, data = subset(froh, id != "C09"))) #sig, so only when excluding small rohs (< 100kb & < 100) do chicks have higher froh
summary(lm(froh ~ age, data = subset(froh, froh < 0.3 & id != "C09")))

write.csv(froh, "output/2_inbreeding/froh_chick_adult.txt", quote=F, row.names=F)
save(froh, file = "output/2_inbreeding/froh_chick_adult.RData")

### summary
summary(froh$froh[which(froh$age == "chick" & froh$id!="C09")])
summary(froh$froh[which(froh$age == "adult")])
