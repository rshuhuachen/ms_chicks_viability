pacman::p_load(data.table, dplyr)
snpef <- "/vol/cluster-data/rchen/git/genetic_load_ltet/results/snpeff/annotation/ltet_ann_chick_snp_output.vcf"

# get header
names <- as.data.frame(t(fread("output/headers_vcf.txt", header = F)))
names <- rownames_to_column(names)
#select only the id names
names <- names[c(10:nrow(names)),]
ids <- fread("/vol/cluster-data/rchen/git/genetic_load_ltet/data/metadata/file_list_all_bgi_clean.csv")
ids$file_id <- gsub(".sorted.bam", "", ids$file)
merged <- full_join(ids, names, by = c("file_id" = "V1"))
merged$order <- gsub("V", "", merged$rowname)
merged$order <- as.numeric(merged$order)
merged <- merged %>% arrange(order)
merged$name <- c(merged$file_id[1:46], merged$id[47:nrow(merged)])
write.table(merged$name, "output/headers_vcf_correct.csv", quote=F, sep=",", row.names=F, col.names=F)

system(paste0("bcftools reheader --samples /vol/cluster-data/rchen/wgr/chicks_genomes/output/headers_vcf_correct.csv -o /vol/cluster-data/rchen/wgr/chicks_genomes/output/grouse_chick_adult.vcf ", snpef))

### add gerp score

all_gerps <- data.frame()
for (i in c(1:30)){
    gerp <- fread(paste0("output/beds/gerp_overlapSNP_scaf_", i, ".tsv.gz"))
    gerp_select <- gerp[,c(1,3,5)]
    all_gerps <- rbind(all_gerps, gerp_select)
}

write_tsv(all_gerps, "output/combined_gerps.tsv", quote=F, row.names=F, col.names=F)

system("ls /vol/cluster-data/rchen/wgr/chicks_genomes/output/beds/*.tsv.gz | xargs zcat | cut -f1,3,3,5 > output/combined_gerps.tsv")

### combine in vcf
zcat {input.vcf} | \
      vcf-annotate -a {input.bed} \
        -d key=INFO,ID=AA,Number=1,Type=String,Description='Ancestral Allele' \
        -c CHROM,FROM,TO,INFO/AA

#index
system("bgzip output/grouse_chick_adult.vcf")
system("tabix -p vcf -h output/grouse_chick_adult.vcf.gz")
system("output/grouse_chick_adult.vcf.gz | vcf-annotate -a output/combined_gerps.bed -d key=INFO,ID=RS,Number=1,Type=Float,Description='GERP score' -c CHROM,FROM,TO,INFO/RS")