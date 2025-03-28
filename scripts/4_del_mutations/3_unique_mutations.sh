#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=8000
#SBATCH --job-name=unique_mut
#SBATCH --nodes=10
#SBATCH --error=log/unique_mut.err
#SBATCH --output=log/unique_mut.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch /vol/cluster-data/rchen/git/ms_chicks_viability/scripts/4_del_mutations/3_unique_mutations.sh

Rscript /vol/cluster-data/rchen/git/ms_chicks_viability/scripts/4_del_mutations/3_unique_mutations.R &> log/unique_mut.out