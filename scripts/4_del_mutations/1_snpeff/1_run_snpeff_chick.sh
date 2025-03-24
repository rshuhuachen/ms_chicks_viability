#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=8000
#SBATCH --job-name=snpef_chick
#SBATCH --nodes=10
#SBATCH --error=log/snpef_chick.err
#SBATCH --output=log/snpef_chick.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch scripts/4_del_mutations/1_snpeff/1_run_snpeff_chick.sh

java -Xmx8g -jar /vol/cluster-data/rchen/git/genetic_load_ltet/src/snpEff.jar ann -v \
  lyrurus_tetrix data/genomic/genomes/processed/ltet_snps_filtered_chicks_adults.vcf.gz \
  > data/genomic/genomes/processed/ltet_ann_chick_snp_output.vcf
  


