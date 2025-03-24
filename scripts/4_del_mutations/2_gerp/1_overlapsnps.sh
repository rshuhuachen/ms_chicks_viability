#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=8000
#SBATCH --job-name=overlap
#SBATCH --nodes=10
#SBATCH --error=log/overlap.err
#SBATCH --output=log/overlap.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch scripts/4_del_mutations/2_gerp/1_overlapsnps.sh

convert2bed -i vcf < output/3_annotated_genome/ltet_filtered_ann_aa_chick_correct.vcf > output/3_annotated_genome/ltet_filtered_ann_aa_chick_correct.bed

for i in {1..30};
    do bedtools intersect -a output/4_load/gerp/beds/gerp_scaf_$i.bed -b output/3_annotated_genome/ltet_filtered_ann_aa_chick_correct.bed -wa -wb | cut -f 6-10 --complement > output/4_load/gerp/overlap/gerp_overlapSNP_scaf_$i.tsv
    gzip output/4_load/gerp/overlap/gerp_overlapSNP_scaf_$i.tsv;
    done
