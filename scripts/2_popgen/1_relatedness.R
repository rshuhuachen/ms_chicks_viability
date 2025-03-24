
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
