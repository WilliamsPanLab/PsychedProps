#!/bin/bash
#
#SBATCH --job-name=Rscript
#SBATCH --time=5:00:00
#SBATCH -n 1
#SBATCH --mem=25G
#SBATCH -p leanew1,normal  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
# will need R
module load R/4.1
Rscript GenTDistr_Spins_VIS.R
Rscript GenTDistr_Spins_MOT.R
Rscript GenTDistr_Spins_FPN.R
