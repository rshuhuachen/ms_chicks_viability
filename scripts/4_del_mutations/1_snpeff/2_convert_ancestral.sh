#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=8000
#SBATCH --job-name=conver_anc
#SBATCH --nodes=10
#SBATCH --error=logs/convert_aa.err
#SBATCH --output=logs/convert_aa.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch scripts/4_del_mutations/1_snpeff/2_convert_ancestral.sh

# extract snp_pos_and_alleles
grep -v "^##" data/genomic/genomes/processed/ltet_ann_chick_snp_output.vcf | 
    cut -f 1,2,4,5 | gzip > data/genomic/genomes/processed/ltet_snps_filtered_chicks.vcf.tsv.gz

# anc alleles

Rscript --vanilla scripts/4_del_mutations/1_snpeff/ancient_alleles_assignment.R 2> logs/r_ancestral_ref.log 1> logs/r_ancestral_ref.log

# remove string ids
cat data/genomic/genomes/processed/ltet_ann_chick_snp_output.vcf | sed 's/\/vol\/cluster-data\/rchen\/wgr\/data\/processed\/alignments_dovetailrefgenome\/sorted_bamfiles\///g' | sed 's/\/vol\/cluster-data\/rchen\/wgr\/chicks_genomes\/output\/processed\///g' | sed 's/.sorted.bam//g' | sed 's/;HRSCAF=/__HRSCAF_/' | gzip > data/genomic/genomes/processed/ltet_snps_filtered_chicks_adults_cleannames.vcf.gz

# clean scaf names in bed
cat data/genomic/genomes/processed/anc_allele_assignment_chick.bed | sed 's/;HRSCAF=/__HRSCAF_/' > data/genomic/genomes/processed/anc_allele_assignment_cleanname_chick.bed

# pack_aa_bed
bgzip data/genomic/genomes/processed/anc_allele_assignment_cleanname_chick.bed
tabix -s 1 -b 2 -e 3 data/genomic/genomes/processed/anc_allele_assignment_cleanname_chick.bed.gz

# annotate_vcf
zcat data/genomic/genomes/processed/ltet_snps_filtered_chicks_adults_cleannames.vcf.gz | \
      vcf-annotate -a data/genomic/genomes/processed/anc_allele_assignment_cleanname_chick.bed.gz \
        -d key=INFO,ID=AA,Number=1,Type=String,Description='Ancestral Allele' \
        -c CHROM,FROM,TO,INFO/AA > data/genomic/genomes/processed/ltet_filtered_ann_aa_chick.vcf

# convert_vcf_alleles -> do this in container
#/prj/blackgrouse/bin/mambaforge/bin/singularity shell --bind /vol/cluster-data/rchen/git/ms_purging_grouse /vol/cluster-data/rchen/git/genetic_load_ltet/scripts/snpeff/1_ancestral_as_reference/workflow/containers/jvarkit_1b2aedf24.sif
#apptainer shell -u --overlay /vol/cluster-data/rchen/git/ms_purging_grouse:ro --fakeroot --bind /vol/cluster-data/rchen/git/ms_purging_grouse /vol/cluster-data/rchen/git/genetic_load_ltet/scripts/snpeff/1_ancestral_as_reference/workflow/containers/cactus_v2.5.1.sif

#java -jar /opt/jvarkit/dist/jvarkit.jar \
#      vcffilterjdk \
#      -f /vol/cluster-data/rchen/git/genetic_load_ltet/scripts/snpeff/1_ancestral_as_reference/workflow/containers/script.js data/genomic/genomes/processeddata/processed/annotated/ltet_filtered_ann_aa_chick.vcf | \
#      bgzip > data/genomic/genomes/processed/ltet_filtered_ann_aa_chick.vcf_wrongscaf.gz

#singularity pull docker://lindenb/jvarkit:1e09f06d4c05e5a148
singularity run --bind /home/nioo/rebeccash/PhD_grouse/ms_purging_grouse jvarkit_1e09f06d4c05e5a148.sif \
    java -jar /opt/jvarkit/dist/jvarkit.jar \
    vcffilterjdk \
    -f scripts/2_snpeff/script.js data/genomic/genomes/processed/ltet_filtered_ann_aa_chick.vcf | bgzip > data/genomic/genomes/processed/ltet_filtered_ann_aa_chick_wrongscaf.vcf.gz

# revert_scaf_name
zcat data/genomic/genomes/processed/ltet_filtered_ann_aa_chick_wrongscaf.vcf.gz | sed 's/__HRSCAF_/;HRSCAF=/' | bgzip > output/3_annotated_genome/ltet_filtered_ann_aa_chick_correct.vcf.gz

# index_aa_vcf
tabix -p vcf output/3_annotated_genome/ltet_filtered_ann_aa_chick_correct.vcf.gz