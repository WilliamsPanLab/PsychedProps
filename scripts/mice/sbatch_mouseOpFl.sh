#!/bin/bash
#
#SBATCH --job-name=MouseOF
#SBATCH --time=48:00:00
#SBATCH -n 1
#SBATCH --mem=30G
#SBATCH -p leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

# will need matlab
module load matlab
# for gif maker
ml system
ml ffmpeg/5.0

# subject name is input argument
subj=$1

# optical flow
matlab -nodisplay -r "Mouse_OpFl_Diaz('$subj',1)"
matlab -nodisplay -r "Mouse_OpFl_Diaz('$subj',2)"
matlab -nodisplay -r "Mouse_OpFl_Diaz('$subj',3)"
matlab -nodisplay -r "Mouse_OpFl_Diaz('$subj',4)"
matlab -nodisplay -r "Mouse_OpFl_Diaz('$subj',5)"
#matlab -nodisplay -r "Mouse_OpFl_Dex('$subj',6)"
