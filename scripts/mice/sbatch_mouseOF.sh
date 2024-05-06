#!/bin/bash
#
#SBATCH --job-name=simStreamComb
#SBATCH --time=48:00:00
#SBATCH -n 1
#SBATCH --mem=25G
#SBATCH -p leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

# will need matlab
module load matlab

# subject name is input argument
subj=$1

matlab -nodisplay -r "Mouse_OF('$subj')"
