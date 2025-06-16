#!/bin/bash
#
#SBATCH --job-name=Rscript
#SBATCH --time=12:00:00
#SBATCH -n 1
#SBATCH --mem=25G
#SBATCH -p leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
# will need R
module load R/4.1
Rscript GenTDistr_Spins_Psil.R 
