#!/bin/bash
#
#SBATCH --job-name=pixelwise_stats
#SBATCH --time=48:00:00
#SBATCH -n 1
#SBATCH --mem=35G
#SBATCH -p leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
module load R/4.1
# run group-level pixelwise
Rscript pixelwise_mixedEfMods_LSD.R
