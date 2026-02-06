#!/bin/bash
#
#SBATCH --job-name=magAgg
#SBATCH --time=48:00:00
#SBATCH -n 1
#SBATCH --mem=25G
#SBATCH -p normal,leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
# will need matlab
module load matlab
matlab -nodisplay -r "loop_Spun_Extracts_FPN"
