
# packages
pacman::p_load(data.table, tidyverse)

#### Kinship Structure ####
## Based on https://github.com/elhumble/Agaz_85K_workflow_2018/blob/master/3.3_IBD.R

#### Additional filtering ####

### Here we define the VCF variables
out = paste0(getwd(), "/output/1_relatedness/")
VCF= paste0(getwd(), "/data/genomic/genomes/processed/", "ltet_chicks_adults.vcf")
VCFPRUNED = paste0(out, "ltet_snps_filtered_chicks_adults_prune")

### Prune with plink
system(paste0("plink --vcf ", VCF, " --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out ", VCFPRUNED))

## Prepare file with plink for ngsrelate
### Additional filtering for LD, high MAF, hwe

# format for ngsRelate
system(paste0("plink --vcf ", VCF, " --double-id --allow-extra-chr --set-missing-var-ids @:# ",
              "--extract ", VCFPRUNED, ".prune.in ",
              "--geno 0.1 --maf 0.01 --hwe 0.001 --genome ",
              "--out ", repo, "ltet_snps_filtered_chicks_adults_ngsrelate ",
              "--recode vcf-iid"))

## format to assign categories
system(paste0("plink --vcf ", VCF, " --double-id --allow-extra-chr --set-missing-var-ids @:# ",
              "--extract ", VCFPRUNED, ".prune.in ",
              "--geno 0.1 --maf 0.01 --hwe 0.001 --genome ",
              "--out ", repo, "ltet_snps_filtered_chicks_adults_pruned_kinship_A ", # = causes an additive component file (0/1/2)
              "--recode A"))

### ngsrelate

system(paste0("/prj/blackgrouse/bin/ngsRelate/ngsRelate -h ", out, "ltet_snps_filtered_chicks_adults_ngsrelate.vcf -T GT ",
              "-O ", out, "ltet_snps_filtered_chicks_adults_ngsrelate.res -c 1"))

### Load in ngsRelate output
ngsrel <- fread(paste0(out, "ltet_snps_filtered_chicks_adults_ngsrelate.res"), fill = TRUE)

gen <- fread(paste0(out, "/ltet_snps_filtered_chicks_adults_pruned_kinship_A.genome"), header = T)

summary(gen$PI_HAT)

gen <- gen %>%
  mutate(kinship = (Z1 / 4) + (Z2 / 2)) %>%
  mutate(criteria = case_when(kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 < 0.1 ~ "Parent-offspring",
                              kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 > 0.1 & Z0 < 0.365 ~ "Full-sibling",
                              kinship >= 1/2^(7/2) & kinship < 1/2^(5/2) & Z0 > 0.365 & Z0 < 1-(1/(2^(3/2))) ~ "Second-degree",
                              kinship >= 1/2^(9/2) & kinship < 1/2^(7/2) & Z0 > 1-(1/2^(3/2)) & Z0 < 1 -(1/2^(5/2)) ~ "Third-degree",
                              kinship < 1/2^(9/2) & Z0 > 1-(1/2^(5/2)) ~ "Unrelated",
                              TRUE ~ "Unknown"))

summary(gen$kinship)                              

unknown <- subset(gen, criteria == "Unknown")
unknown$kinship

#C09 and D229192 are probably the same individual

#### Clean up for checking parents of chicks ####

# select relevant only
clean_gen <- gen%>%select(c(IID1, IID2, Z0, Z1, PI_HAT, kinship, criteria))

# change names
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

# clean
clean_gen$IID1 <- gsub("/vol/cluster-data/rchen/wgr/data/processed/alignments_dovetailrefgenome/sorted_bamfiles/", "", clean_gen$IID1)
clean_gen$IID2 <- gsub("/vol/cluster-data/rchen/wgr/data/processed/alignments_dovetailrefgenome/sorted_bamfiles/", "", clean_gen$IID2)
clean_gen$IID1 <- gsub("/vol/cluster-data/rchen/wgr/chicks_genomes/output/processed/", "", clean_gen$IID1)
clean_gen$IID2 <- gsub("/vol/cluster-data/rchen/wgr/chicks_genomes/output/processed/", "", clean_gen$IID2)

clean_gen <- left_join(clean_gen, ids[,c("file", "id")], by = c("IID1" = "file"))
clean_gen <- left_join(clean_gen, ids[,c("file", "id")], by = c("IID2" = "file"))

clean_gen <- clean_gen %>% mutate(ID1 = case_when(is.na(id.x) ~ IID1,
                                                  !is.na(id.x) ~ id.x),
                                  ID2 = case_when(is.na(id.y) ~ IID2,
                                                  !is.na(id.y) ~ id.y))

clean_gen <- clean_gen %>% select(ID1, ID2, Z0:criteria)
clean_gen$ID1 <- gsub(".sorted.bam", "", clean_gen$ID1)
clean_gen$ID2 <- gsub(".sorted.bam", "", clean_gen$ID2)

gen_chick <- subset(clean_gen, (grepl("C", clean_gen$ID1) | grepl("C", clean_gen$ID2)) & criteria == "Parent-offspring")

# in reverse to see if it matches
gen_adult <- subset(clean_gen, (grepl("D", clean_gen$ID1) | grepl("D", clean_gen$ID2)) & criteria == "Parent-offspring")
