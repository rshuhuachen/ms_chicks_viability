#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=8000
#SBATCH --job-name=fasqc
#SBATCH --nodes=1
#SBATCH --error=/vol/cluster-data/rchen/wgr/chicks_genomes/log/fastqc.err
#SBATCH --output=/vol/cluster-data/rchen/wgr/chicks_genomes/log/fastqc.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch /vol/cluster-data/rchen/wgr/chicks_genomes/workflow/fastqc.sh

for file in /vol/cluster-data/rchen/wgr/chicks_genomes/data/BMK231125-BT148-ZX01-0201/BMK_DATA_20231229114919_1/Data/Unknown_*; 
do
    mv "$file" "${file/Unknown_BT148-00100/C}"; 
done

for file in /vol/cluster-data/rchen/wgr/chicks_genomes/data/BMK231125-BT148-ZX01-0201/BMK_DATA_20231229114919_1/Data/*fq.gz; 
do
    /vol/biotools/bin/fastqc $file --threads 8 -o /vol/cluster-data/rchen/wgr/chicks_genomes/output/report/fastqc;
done
