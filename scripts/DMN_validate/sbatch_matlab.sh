#!/bin/bash
#
#SBATCH --job-name=FPNSpins
#SBATCH --time=48:00:00
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH -p leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
# will need matlab
module load matlab
#matlab -nodisplay -r "MOTSpins2mat_fslr"
#matlab -nodisplay -r "VISSpins2mat_fslr"
matlab -nodisplay -r "FPNSpins2mat_fslr"
