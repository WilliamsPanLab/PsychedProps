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
# and freesurfer
module load biology
module load freesurfer/7.3.2
# and workbench
module load workbench

# subject name is input argument
subj=$1
# sesh is input 2
sesh=$2

matlab -nodisplay -r "ExampleFrames_MDMA('$subj','$sesh','rs1')"
matlab -nodisplay -r "ExampleFrames_MDMA('$subj','$sesh','rs2')"

# run matlab
#matlab -nodisplay -r "Lagged_CircFC('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Lagged_CircFC('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Lagged_CircFC('$subj','$sesh','wm')"
#matlab -nodisplay -r "Lagged_CircFC('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Lagged_CircFC('$subj','$sesh','emotion')"
