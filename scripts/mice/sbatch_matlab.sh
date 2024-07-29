#!/bin/bash
#
#SBATCH --job-name=simStreamComb
#SBATCH --time=3:00:00
#SBATCH -n 1
#SBATCH --mem=25G
#SBATCH -p leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

# dumb cd
#cd /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/mice/
# will need matlab
module load matlab
matlab -nodisplay -r "Extract_All_angles_Diaz"
