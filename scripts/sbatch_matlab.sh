#!/bin/bash
#
#SBATCH --job-name=magAgg
#SBATCH --time=3:00:00
#SBATCH -n 1
#SBATCH --mem=25G
#SBATCH -p leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
subj=$1
# will need matlab
module load matlab
matlab -nodisplay -r "Calc_AvgMagnitude_psil('$subj')"
