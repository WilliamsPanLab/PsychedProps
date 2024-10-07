#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=1:00:00
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH -p leanew1,normal  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=no@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

# will need matlab
module load matlab
# and freesurfer
module load biology
module load freesurfer/7.3.2

# subject name is input argument
subj=$1

# cd to scripts directory
cd /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts

# Center and bandpass the data
matlab -nodisplay -r "partialC_and_BP_and_TDif_ISPOT('$subj')"

