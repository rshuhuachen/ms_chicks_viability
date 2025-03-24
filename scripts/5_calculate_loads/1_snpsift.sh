#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=8000
#SBATCH --job-name=snpsift
#SBATCH --nodes=10
#SBATCH --error=/vol/cluster-data/rchen/git/ms_purging_grouse/logs/snpsift.err
#SBATCH --output=/vol/cluster-data/rchen/git/ms_purging_grouse/logs/snpsift.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch /vol/cluster-data/rchen/git/ms_purging_grouse/scripts/5_del_mutations/1_snpsift.sh

### Filter for LOF mutations and high impact mutations ####


## High impact
zcat output/3_annotated_genome/ltet_filtered_ann_aa_chick_correct.vcf.gz| java -jar /vol/cluster-data/rchen/git/genetic_load_ltet/src/SnpSift.jar filter " ( ANN[*].IMPACT = 'HIGH' )" > output/4_load/snpeff/ltet_filtered_ann_aa_chick_high.vcf
bgzip output/4_load/snpeff/ltet_filtered_ann_aa_chick_high.vcf

## medium
zcat output/3_annotated_genome/ltet_filtered_ann_aa_chick_correct.vcf.gz | java -jar /vol/cluster-data/rchen/git/genetic_load_ltet/src/SnpSift.jar filter " ( ANN[*].IMPACT = 'MODERATE')" > output/4_load/snpeff/ltet_filtered_ann_aa_chick_moderate.vcf
bgzip output/4_load/snpeff/ltet_filtered_ann_aa_chick_moderate.vcf




