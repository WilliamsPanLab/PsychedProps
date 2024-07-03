#!/bin/bash
#
#SBATCH --job-name=MouseOF
#SBATCH --time=100:00:00
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
matlab -nodisplay -r "Mouse_OpFl_ket_mixed('$subj',1)"
matlab -nodisplay -r "Mouse_OpFl_ket_mixed('$subj',2)"
matlab -nodisplay -r "Mouse_OpFl_ket_mixed('$subj',3)"
matlab -nodisplay -r "Mouse_OpFl_ket_mixed('$subj',4)"
matlab -nodisplay -r "Mouse_OpFl_ket_mixed('$subj',5)"
