#### Load packages ####

pacman::p_load(dplyr, data.table, readr)

#source("scripts/theme_ggplot.R")

#### Load data #####

## gerp
load(file = "output/4_load/gerp/gerps_all.RData")

## snpeff
high <- read.table("output/4_load/snpeff/ltet_filtered_ann_aa_chick_high.vcf.gz")

# load headers
names <- fread("output/2_inbreeding/chicks_adults_samples.txt", header = F)
names(names) <- "file_id"

# load id info
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

#### Unique mutations GERP, have to calculate allele frequencies based on equal sample sizes of adults and chicks ####

#### Prepare files ####

### get a file just noting down the chr and pos of deleterious mutations for GERP and SnpEff separately

# gerp > 4
#write_tsv(data.frame(chr = gerp$chr, pos = gerp$pos), file = "output/4_load/gerp/chr_pos_gerp4.txt", col_names = F)

# high impact
#write_tsv(data.frame(chr = high$CHROM, pos = high$POS), file = "output/4_load/snpeff/chr_pos_high.txt", col_names = F)

#### Function to get AF ####

get_unique_mutations <- function(n_random_draws, n_ids, vcf, out){
  # load adult and chick ids
  adult_ids <- fread("data/metadata/ids_adults.txt", header=F)
  chicks_ids <- fread("data/metadata/ids_chicks.txt", header=F)
  summary_gerp <- data.frame()
  summary_high <- data.frame()
  
  for (i in 1:n_random_draws){
    ## take X number of random adult ids
    adult_ids_subset <- sample_n(adult_ids, size = n_ids, replace = F)
    adult_ids_subset$V1 <- paste0("/vol/cluster-data/rchen/wgr/data/processed/alignments_dovetailrefgenome/sorted_bamfiles/", adult_ids_subset$V1, ".sorted.bam")
    write_tsv(as.data.frame(adult_ids_subset$V1), file = paste0("output/3_annotated_genome/unique_mut/adults_ids_subset_",i,".txt"), col_names = F)
    
    ## take X number of random chick ids
    chicks_ids <- chicks_ids[which(chicks_ids$V1 != "C09")]
    chick_ids_subset <- sample_n(chicks_ids, size = n_ids, replace = F)
    chick_ids_subset$V1 <- paste0("/vol/cluster-data/rchen/wgr/chicks_genomes/output/processed/", chick_ids_subset$V1, ".sorted.bam")
    
    write_tsv(as.data.frame(chick_ids_subset$V1), file = paste0("output/3_annotated_genome/unique_mut/chicks_ids_subset_",i,".txt"), col_names = F)
    
    ## calculate allele frequencies based on the two random subsets for adults and chicks separately
    
    ### gerp ####
    # gerp adults
    system(paste0("vcftools --vcf ", vcf, " --positions output/4_load/gerp/chr_pos_gerp4.txt --keep output/3_annotated_genome/unique_mut/adults_ids_subset_",
    i,".txt --freq --out output/3_annotated_genome/unique_mut/allele_freq_gerp4_adults_subset_", i))
    
    # gerp chicks
    system(paste0("vcftools --vcf ", vcf, " --positions output/4_load/gerp/chr_pos_gerp4.txt --keep output/3_annotated_genome/unique_mut/chicks_ids_subset_",
                  i,".txt --freq --out output/3_annotated_genome/unique_mut/allele_freq_gerp4_chicks_subset_", i))
    
    
    ## load in data, merge and clean up
    gerp4_ad <- fread(paste0("output/3_annotated_genome/unique_mut/allele_freq_gerp4_adults_subset_", i, ".frq"), 
                      skip = 1, col.names = c("scaf",  "pos",  "n_alleles",  "n_chr",  "af_1", "af_2"))
    
    gerp4_ch <- fread(paste0("output/3_annotated_genome/unique_mut/allele_freq_gerp4_chicks_subset_", i, ".frq"), 
                      skip = 1, col.names = c("scaf",  "pos",  "n_alleles",  "n_chr",  "af_1", "af_2"))
    
    # merge adult and chick af
    gerp4_af <- left_join(gerp4_ad[,c("scaf", "pos", "af_1", "af_2")], gerp4_ch[,c("scaf", "pos", "af_1", "af_2")], 
                          by = c("scaf", "pos"), suffix = c("_ad", "_chick"))
    gerp4_af$allele_1 <- substr(gerp4_af$af_1_ad,1,1)
    gerp4_af$allele_2 <- substr(gerp4_af$af_2_ad,1,1)
    
    gerp4_af$af_1_ad <- as.numeric(substr(gerp4_af$af_1_ad, 3, nchar(gerp4_af$af_1_ad)))
    gerp4_af$af_2_ad <- as.numeric(substr(gerp4_af$af_2_ad, 3, nchar(gerp4_af$af_2_ad)))
    gerp4_af$af_1_chick <- as.numeric(substr(gerp4_af$af_1_chick, 3, nchar(gerp4_af$af_1_chick)))
    gerp4_af$af_2_chick <- as.numeric(substr(gerp4_af$af_2_chick, 3, nchar(gerp4_af$af_2_chick)))
    
    # remove those with NAs
    gerp4_af <- subset(gerp4_af, !is.na(af_1_chick))
    
    # merge with anc and derived alleles
    gerp4_af <- left_join(gerp4_af, gerp[,c("chr", "pos", "ancestral", "derived")], by = c("scaf" = "chr", "pos"))
    
    gerp4_af <- gerp4_af %>% mutate(
      af_der_ad = case_when(allele_1 == derived ~ af_1_ad,
                            allele_2 == derived ~ af_2_ad),
      af_der_chick = case_when(allele_1 == derived ~ af_1_chick,
                               allele_2 == derived ~ af_2_chick))
    
    # make a category whether it is a unique mutation or not
    gerp4_af <- gerp4_af %>% mutate(
      unique = as.factor(case_when(af_der_ad > 0 & af_der_chick == 0 ~ "unique_in_ad",
                                   af_der_chick > 0 & af_der_ad == 0 ~ "unique_in_chick")))
      
    ### SnpEff ####
      
    # snpeff adults
    system(paste0("vcftools --vcf ", vcf, " --positions output/4_load/snpeff/chr_pos_high.txt --keep output/3_annotated_genome/unique_mut/adults_ids_subset_",
                  i,".txt --freq --out output/3_annotated_genome/unique_mut/allele_freq_high_adults_subset_", i))
    
    # snpeff chicks
    system(paste0("vcftools --vcf ", vcf, " --positions output/4_load/snpeff/chr_pos_high.txt --keep output/3_annotated_genome/unique_mut/chicks_ids_subset_",
                  i,".txt --freq --out output/3_annotated_genome/unique_mut/allele_freq_high_chicks_subset_", i))
  
    
    ## load in data, merge and clean up
    high_ad <- fread(paste0("output/3_annotated_genome/unique_mut/allele_freq_high_adults_subset_", i, ".frq"), 
                      skip = 1, col.names = c("scaf",  "pos",  "n_alleles",  "n_chr",  "af_1", "af_2"))
    
    high_ch <- fread(paste0("output/3_annotated_genome/unique_mut/allele_freq_high_chicks_subset_", i, ".frq"), 
                      skip = 1, col.names = c("scaf",  "pos",  "n_alleles",  "n_chr",  "af_1", "af_2"))
    
    # merge adult and chick af
    high_af <- left_join(high_ad[,c("scaf", "pos", "af_1", "af_2")], high_ch[,c("scaf", "pos", "af_1", "af_2")], 
                          by = c("scaf", "pos"), suffix = c("_ad", "_chick"))
    high_af$allele_1 <- substr(high_af$af_1_ad,1,1)
    high_af$allele_2 <- substr(high_af$af_2_ad,1,1)
    
    high_af$af_1_ad <- as.numeric(substr(high_af$af_1_ad, 3, nchar(high_af$af_1_ad)))
    high_af$af_2_ad <- as.numeric(substr(high_af$af_2_ad, 3, nchar(high_af$af_2_ad)))
    high_af$af_1_chick <- as.numeric(substr(high_af$af_1_chick, 3, nchar(high_af$af_1_chick)))
    high_af$af_2_chick <- as.numeric(substr(high_af$af_2_chick, 3, nchar(high_af$af_2_chick)))
    
    # remove those with NAs
    high_af <- subset(high_af, !is.na(af_1_chick))
    
    # merge with anc and derived alleles
    high_af <- left_join(high_af, high[,c("CHROM", "POS", "REF", "ALT")], by = c("scaf" = "CHROM", "pos" = "POS"))
    high_af <- high_af %>% rename(ancestral = REF, derived = ALT)
    
    high_af <- high_af %>% mutate(
      af_der_ad = case_when(allele_1 == derived ~ af_1_ad,
                            allele_2 == derived ~ af_2_ad),
      af_der_chick = case_when(allele_1 == derived ~ af_1_chick,
                               allele_2 == derived ~ af_2_chick))
    
    # make a category whether it is a unique mutation or not
    high_af <- high_af %>% mutate(
      unique = as.factor(case_when(af_der_ad > 0 & af_der_chick == 0 ~ "unique_in_ad",
                                   af_der_chick > 0 & af_der_ad == 0 ~ "unique_in_chick")))
    
    # calculate Rxy
    gerp4_af <- gerp4_af %>% mutate(l_ac = af_der_ad*(1-af_der_chick),
                                    l_ca = af_der_chick*(1-af_der_ad))
    
    high_af <- high_af %>% mutate(l_ac = af_der_ad*(1-af_der_chick),
                                    l_ca = af_der_chick*(1-af_der_ad))
    
    
  # summarise result for this iteration
  summary_gerp_i <- data.frame(sample = i,
                               method = "gerp",
                               n_unique_in_ad = nrow(subset(gerp4_af, unique == "unique_in_ad")),
                               n_unique_in_chick = nrow(subset(gerp4_af, unique == "unique_in_chick")),
                               rxy = sum(gerp4_af$l_ac) / sum(gerp4_af$l_ca))  
  
  summary_high_i <- data.frame(sample = i,
                               method = "high",
                               n_unique_in_ad = nrow(subset(high_af, unique == "unique_in_ad")),
                               n_unique_in_chick = nrow(subset(high_af, unique == "unique_in_chick")),
                               rxy = sum(high_af$l_ac) / sum(high_af$l_ca))  
  
  summary_gerp <- rbind(summary_gerp, summary_gerp_i)
  summary_high <- rbind(summary_high, summary_high_i)
  
  }
  
  
  
  summary <- rbind(summary_gerp, summary_high)
  save(summary, file = out)
  return(summary)
}

sum_unique_mutations <- get_unique_mutations(n_random_draws = 500,
                                             n_ids = 40,
                                             vcf = "data/genomic/genomes/processed/ltet_ann_chick_snp_output.vcf",
                                             out = "output/3_annotated_genome/summary_unique_mutations_40ids.RData")



