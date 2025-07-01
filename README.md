# Viability selection acts against exonic conserved mutations in early life

In this repository, you will find all the scripts, metadata, processed data and plots for the manuscript Viability selection acts against exonic conserved mutations in early-life (in preparation) by Chen et al.

Below you will find a quick guide how to navigate through this repository

## Data

The data that you can find in this repository are the phenotypic and metadata. All genomic data can be downloaded from NCBI under BioProject [PRJNA1085187](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1085187) with SRA BioAccession Numbers XXXX (chicks) and [40722954–40723143](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=1085187) (yearlings and adults). Please see [Chen et al. (2025)](https://github.com/rshuhuachen/ms_load_grouse) for details on the reference genome ([GCA_043882375.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_043882375.1/)) and [genome annotation](https://github.com/rshuhuachen/ms_load_grouse/tree/main/data/genomic/annotation).

The `./data/metadata` subdirectory contains information about the 30 largest scaffolds ("30_largest.scaf"") such as their size and ID, raw and processed files about the chicks (such as the sampling data, chick ID, etc. in files "Chick samples to send"", "Chick_sample_list_extractions" and "alldata_chicks""), the script used to process the raw metadata files into the clean, processed ones ("combine_metadata"), data on the file names of the raw fasta files ("file_list_all_bgi_clean" for adults), a list of chick and yearling + adult IDs ("ids_chicks" and "ids_adults", respectively) and the clean, processed metadata file ("metadata_adult_chick" which is an .RData object that can be imported in R using the load() command).

The `./data/phenotypic` subdirectory contains an .RData object that contains phenotypic information of yearlings and adults, including longevity which is the variable used for this manuscript (longevity == 1 year; yearling versus longevity \> 1; adult).

## Scripts

The `./scripts` subdirectory contains all scripts used for pre-processing, analyses and plotting figures. The subdirectories and scripts are numbered in the order they should and have been executed in.

1_genotyping: The scripts in this folder first perform a quality check on the raw .fq.gz files (`./scripts/1_genotyping/0_fastqc.sh`). If you want to reproduce these steps, you will need to download the original sequencing files and adjust the paths used in the scripts according to where you place the downloaded files. After the quality check, we use snakemake (`./scripts/1_genotyping/1_snakefile`) to align our genomes to the reference genome and prepare the files for SNPs calling (`./scripts/1_genotyping/2_snpcall.sh`). We get a depth report in bash (`./scripts/1_genotyping/3_get_depthreport.sh`) which we next analyze in R (`./scripts/1_genotyping/4_get_depthreport.R`). Next, we use the coverage statistics to filter our SNPs accordingly (`./scripts/1_genotyping/5_snpfilter.sh`). Lastly, we count the number of raw reads for reporting purposes (`./scripts/1_genotyping/6_count_reads.R`)

The data are output in the directory `data/genomic/genomes/processed/` but due to the size of the files, these are not uploaded on github. This includes the .vcf file, for example.

2_popgen: Here, we calculate relatedness in order to identify recaptured individuals.

3_inbreeding: In this directory we use BCFtools to identify ROHs and calculate FROH

4_del_mutations: Here, we compute the predicted deleterious mutations using SnpEff and GERP. After executing SnpEff, we assign the ancestral alleles based on the reconstructed genome of the last common ancestor of the black grouse and the white ptarmigan. This file is also used to calculate the GERP scores.

5_calculate_loads: We next use the predicted deleterious mutations to calculate the total mutation loads. First, we filter for the deleterious mutations only (high impact for SnpEff, GERP scores \>= 4 for GERP) and use a function specified in `scripts/5_calculate_loads/function_calculate_load.R` to calculate the loads. The loads are only calculated based on the 29 largest autosomal scaffolds and exclude mutations that contain warning messages (from SnpEff). To calculate the total load stratified by genomic region, we first annotate each mutation according to its location, subset the mutations accordingly, and calculate the total mutation load based on each subset.

6_models: In this directory, we use the package brms to compute Bayesian linear mixed effect models. We first run the models for inbreeding, then the total load (genome wide), followed by region-specific load. We region-specific load models are executed through snakemake which can be called with the script in `./scripts/6_models/snakefile_models_region`

7_plotting: Lastly, we run model diagnostics on all models, extract the results, and plot them accordingly. This is done separately for figure 1 (`./scripts/7_plotting/f1_plot_froh_load_age.R`, although the grouse illustrations have been added manually with Figma), figure 2 (`./scripts/7_plotting/f2_plot_load_regions.R` and the supplementary figures (`./scripts/7_plotting/plot_supps.R`

## Output

The majority of the output (e.g. model output) is too large to be stored on github. However, what can be found on github are the model diagnosis plots (`./output/model_diagnosis`), FROH estimates (`./output/froh_chick_adult.RData`), the processed individual load estimates (`./output/loads.RData` and `./output/loads_per_region.RData`) and the intervals of the Bayesian model outputs plus R2 values (`./output/intervals*`).

## Package versions

The genotyping analyses should be executed on a HPU, and due to the resource demanding process of the brms models, I'd recommend to also run these on a HPU. All analyses were done with bash and/or R(Studio). The package versions can be found within the manuscript.
