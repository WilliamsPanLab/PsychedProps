#!/bin/bash
#
#SBATCH --job-name=ttest
#SBATCH --time=16:30:00
#SBATCH -n 1
#SBATCH --mem=20G
#SBATCH -p leanew1,normal  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=no@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
module load R
Rscript facewise_ttest_mdma.R

