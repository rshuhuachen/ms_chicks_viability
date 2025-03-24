#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=12000
#SBATCH --job-name=snpcall
#SBATCH --nodes=1
#SBATCH --error=log/snpcall.err
#SBATCH --output=log/snpcall.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch scripts/1_genotyping/1_snpcall.sh

bcftools mpileup -f data/genomic/refgenome/PO2979_Lyrurus_tetrix_black_grouse.RepeatMasked.fasta -a DP -C 50 -q 20 -b data/genomic/genomes/list_bam_all.txt -I --threads 12 | bcftools call -Ov -v -m -o data/genomic/genomes/processed/ltet_chicks_adults.vcf
