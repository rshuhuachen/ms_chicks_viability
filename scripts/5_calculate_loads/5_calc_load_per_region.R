#### load packages ####
pacman::p_load(genomation, data.table, GenomicFeatures, rtracklayer, fuzzyjoin, tibble, dplyr, GenomicRanges, tidyverse)

# function to calculate load
source("scripts/function_calculate_load.R")

#### GERP ####
load("output/4_load/gerp/gerps_all.RData")

##### convert to vcf format ####
gerps_vcf <- gerp %>% 
  dplyr::select(c(chr, pos, ancestral, derived, qual, info, format, C01:D229777)) #exclude gerp scores

gerps_vcf <- gerps_vcf %>% mutate(ID = ".", .after = pos) %>% mutate(FILTER = ".", .after = qual)

names(gerps_vcf)[1:9] <- c("CHROM", "POS", "ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT") 

### change scaf names
gerps_vcf$CHROM <- gsub(";", "__", gerps_vcf$CHROM)
gerps_vcf$CHROM <- gsub("=", "_", gerps_vcf$CHROM)

### remove genotypes: not necessary here
gerp_pos <- gerps_vcf[,c(1:9)]

##### load annotation data #####

annotation_dir <- "/vol/cluster-data/rchen/git/grouse-annotation/output"
#annotation_dir <- "/Users/vistor/Documents/Work/GitHub/PhD/grouse-annotation/output"

# extract exons
promoter=unique(gffToGRanges(paste0(annotation_dir, "/promoters.gff3")))
exons_gene=unique(gffToGRanges(paste0(annotation_dir, "/exons_gene.gff3")))
introns=unique(gffToGRanges(paste0(annotation_dir, "/introns_transcripts.gff3")))

##### annotate gerp regions ####
gerp_pos$end <- gerp_pos$POS
gerp_pos$start <- gerp_pos$POS
gerp_gr <- as(gerp_pos, "GRanges")

gerp_promoter <- as.data.frame(subsetByOverlaps(gerp_gr, promoter)) %>%
  add_column("region_promoter" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

gerp_exons <- as.data.frame(subsetByOverlaps(gerp_gr, exons_gene))%>%
  add_column("region_exon" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

gerp_intron <- as.data.frame(subsetByOverlaps(gerp_gr, introns))%>%
  add_column("region_intron" = "1") %>% dplyr::rename(chr = seqnames)%>% unique()

gerp_regions <- full_join(gerp_promoter[,c("chr", "start", "region_promoter")], gerp_exons[,c("chr", "start", "region_exon")], by = c("chr","start"))
gerp_regions <- full_join(gerp_regions, gerp_intron[,c("chr", "start", "region_intron")], by = c("chr","start"))

##### merge with genotypes ####

gerp$chr <- gsub(";", "__", gerp$chr)
gerp$chr <- gsub("=", "_", gerp$chr)

gerp_regions_geno <- left_join(gerp, gerp_regions, 
                               by = c(c("chr"), c("pos" = "start")))

##### save annotated gerps ####
save(gerp_regions_geno, file = "output/4_load/gerp/gerp_annotated_regions.RData")

##### calculate gerp loads per region ####
gerp_promoter <- subset(gerp_regions_geno, region_promoter ==1)
gerp_exon <- subset(gerp_regions_geno, region_exon ==1)
gerp_intron <- subset(gerp_regions_geno, region_intron==1)

gerp_promoter <- gerp_promoter[,c(1:(ncol(gerp_promoter)-3))]
gerp_exon <- gerp_exon[,c(1:(ncol(gerp_exon)-3))]
gerp_intron <- gerp_intron[,c(1:(ncol(gerp_intron)-3))]

load_exon <- calculate_load_gerp(vcf = gerp_exon, output_vcf = F, loadtype = "gerp_exons")
load_intron <- calculate_load_gerp(vcf = gerp_intron, output_vcf = F, loadtype = "gerp_introns")
load_promo <- calculate_load_gerp(vcf = gerp_promoter, output_vcf = F, loadtype = "gerp_promoters")

loads_regions_gerp <- rbind(load_exon, load_intron, load_promo)

#### SnpEff #####
high_gt <- read.table(file="output/4_load/snpeff/ltet_filtered_ann_aa_chick_high.vcf.gz")

##### rename columns #####
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

names(high_gt) <- c(c("CHROM", "POS", "ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT"), names$id)

### change scaf names
high_gt$CHROM <- gsub(";", "__", high_gt$CHROM)
high_gt$CHROM <- gsub("=", "_", high_gt$CHROM)

### remove genotypes: not necessary here
high <- high_gt[,c(1:9)]

##### annotate high impact mutations ####
high$end <- high$POS
high$start <- high$POS
high_gr <- as(high, "GRanges")

high_promoter <- as.data.frame(subsetByOverlaps(high_gr, promoter)) %>%
  add_column("region_promoter" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

high_exons <- as.data.frame(subsetByOverlaps(high_gr, exons_gene))%>%
  add_column("region_exon" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

high_intron <- as.data.frame(subsetByOverlaps(high_gr, introns))%>%
  add_column("region_intron" = "1") %>% dplyr::rename(chr = seqnames)%>% unique()

high_regions <- full_join(high_promoter[,c("chr", "start", "region_promoter")], high_exons[,c("chr", "start", "region_exon")], by = c("chr","start"))
high_regions <- full_join(high_regions, high_intron[,c("chr", "start", "region_intron")], by = c("chr","start"))

##### merge with geno ####

high_regions_geno <- left_join(high_gt, high_regions, 
                               by = c(c("CHROM" = "chr"), c("POS" = "start")))

##### save annotated high regions ####
save(high_regions_geno, file = "output/4_load/snpeff/highs_annotated_regions.RData")

##### calculate loads ####
high_regions_geno$CHROM <- gsub("__",";", high_regions_geno$CHROM)
high_regions_geno$CHROM <- gsub("HRSCAF_", "HRSCAF=", high_regions_geno$CHROM)

high_promoter <- subset(high_regions_geno, region_promoter ==1)
high_exon <- subset(high_regions_geno, region_exon ==1)
high_intron <- subset(high_regions_geno, region_intron==1)

high_promoter <- high_promoter[,c(1:(ncol(high_promoter)-3))]
high_exon <- high_exon[,c(1:(ncol(high_exon)-3))]
high_intron <- high_intron[,c(1:(ncol(high_intron)-3))]

load_exon_high <- calculate_load_snpeff(vcf = high_exon, output_vcf = F, loadtype = "high_exons")
load_intron_high <- calculate_load_snpeff(vcf = high_intron, output_vcf = F, loadtype = "high_introns")
load_promo_high <- calculate_load_snpeff(vcf = high_promoter, output_vcf = F, loadtype = "high_promoters")

loads_regions_high <- rbind(load_exon_high, load_intron_high, load_promo_high)

#### combine loads into one file ####
load_per_region <- rbind(loads_regions_gerp, loads_regions_high)

save(load_per_region, file = "output/loads_per_region.RData")
