#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=8000
#SBATCH --job-name=snpsift
#SBATCH --nodes=10
#SBATCH --error=log/snpsift.err
#SBATCH --output=log/snpsift.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch scripts/5_calculate_loads/2_get_gerp.sh

Rscript scripts/5_calculate_loads/2_get_gerp.R &> log/snpsift.out