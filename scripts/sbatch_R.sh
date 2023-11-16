#!/bin/bash
#
#SBATCH --job-name=facewise_stats
#SBATCH --time=48:00:00
#SBATCH -n 1
#SBATCH --mem=35G
#SBATCH -p leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
module load R/4.1
# iteration is only argument
Rscript facewise_Task_L.R $1
