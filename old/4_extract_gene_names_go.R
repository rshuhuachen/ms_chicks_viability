#### load packages ####
pacman::p_load(genomation, data.table, GenomicFeatures, rtracklayer, fuzzyjoin, tibble, dplyr, GenomicRanges, tidyverse)

### GERP annotated ####
load(file = "output/4_load/gerp/gerp_annotated_regions.RData")

gerp_exons <- subset(gerp_regions_geno, region_exon ==1)

### remove genotypes: not necessary here
gerp_exons_snps <- gerp_exons[,c(1:9)]

##### load annotation data #####

annotation_dir <- "/vol/cluster-data/rchen/git/grouse-annotation/output"

# extract promo
exons_gene=unique(gffToGRanges(paste0(annotation_dir, "/exons_gene.gff3")))

##### annotate gerp regions ####
gerp_exons_snps$end <- gerp_exons_snps$pos
gerp_exons_snps$start <- gerp_exons_snps$pos
gerp_gr <- as(gerp_exons_snps, "GRanges")

gerp_exons_snps_full <- mergeByOverlaps(gerp_gr, exons_gene)

gene_ids_exon <- data.frame(ann = unique(gerp_exons_snps_full$ID))

#### merge with 'similar to' column to get gene IDs ####
lookup <- fread("../grouse-annotation/data/lookup_ANN_gene_id.txt")
names(lookup) <- c("ann", "gene_id")

gene_ids_exon <- left_join(gene_ids_exon, lookup, by = "ann")
gene_ids_exon$gene_id <- toupper(gene_ids_exon$gene_id)

write.csv(unique(gene_ids_exon$gene_id), file = "output/gene_ids_gerp_exon_mutations.csv", quote=F, row.names = F)

### background list ####
# mutations in all exons

all <- read.table("data/genomic/genomes/processed/ltet_snps_filtered_chicks_adults.vcf.gz") %>% dplyr::select(c(V1:V9))

names(all)[1:9] <- c("CHROM", "POS", "ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT") 

### change scaf names
all$CHROM <- gsub(";", "__", all$CHROM)
all$CHROM <- gsub("=", "_", all$CHROM)

##### load annotation data #####

annotation_dir <- "/vol/cluster-data/rchen/git/grouse-annotation/output"

# extract exons
exons_gene=unique(gffToGRanges(paste0(annotation_dir, "/exons_gene.gff3")))

##### annotate gerp regions ####
all$end <- all$POS
all$start <- all$POS
all_gr <- as(all, "GRanges")

all_exons <- mergeByOverlaps(all_gr, exons_gene)

data.frame(ann = unique(all_exons$exon_name))
all_exons_ids <- data.frame(ann = c(unique(all_exons$exon_name)))
all_exons_ids$ann <- sub("(-RA).*", "", all_exons_ids$ann)

all_exons_ids <- left_join(all_exons_ids, lookup, by = "ann")
all_exons_ids$gene_id <- toupper(all_exons_ids$gene_id)

write.csv(unique(all_exons_ids$gene_id), file = "output/gene_ids_all_exon_mutations.csv", quote=F, row.names = F)

