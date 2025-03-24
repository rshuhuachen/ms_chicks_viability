cd /vol/cluster-data/rchen/git/ms_chicks_viability
VCF="/vol/cluster-data/rchen/git/ms_chicks_viability/data/genomic/genomes/processed/ltet_chicks_adults.vcf"

## filtering indels
minallele=2
maxallele=2

#total raw snps =  9,025,598 -> grep -v "#" data/genomic/genomes/processed/ltet_chicks_adults.vcf | wc -l
## one liner filtering step

# make sure to change the paths where you have the unfiltered vcf file as well as the filtered one

vcftools --vcf $VCF --min-alleles 2 --max-alleles 2 --minDP 20 --max-missing 0.7 --max-meanDP 60 --minQ 30 --recode --stdout  > data/genomic/genomes/ltet_snps_filtered_chicks_adults.vcf
bgzip data/processed/ltet_snps_filtered_chicks_adults.vcf
#7,245,059
