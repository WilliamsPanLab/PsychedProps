#!/bin/bash
#
#SBATCH --job-name=MouseOF
#SBATCH --time=3:00:00
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

# extract DMN segregation (LSD)
#matlab -nodisplay -r "Extract_DMNSeg_mice('$subj',1)"
#matlab -nodisplay -r "Extract_DMNSeg_mice('$subj',2)"
#matlab -nodisplay -r "Extract_DMNSeg_mice('$subj',3)"
#matlab -nodisplay -r "Extract_DMNSeg_mice('$subj',4)"
#matlab -nodisplay -r "Extract_DMNSeg_mice('$subj',5)"
#matlab -nodisplay -r "Extract_DMNSeg_mice('$subj',6)"

# DMN autocorrelation
matlab -nodisplay -r "Extract_TAutoCor_mice('$subj',1)"
matlab -nodisplay -r "Extract_TAutoCor_mice('$subj',2)"
matlab -nodisplay -r "Extract_TAutoCor_mice('$subj',3)"
matlab -nodisplay -r "Extract_TAutoCor_mice('$subj',4)"
matlab -nodisplay -r "Extract_TAutoCor_mice('$subj',5)"
matlab -nodisplay -r "Extract_TAutoCor_mice('$subj',6)"

# optical flow
#matlab -nodisplay -r "Extract_DMNMag_Diaz_mice('$subj',1)"
#matlab -nodisplay -r "Extract_DMNMag_Diaz_mice('$subj',2)"
#matlab -nodisplay -r "Extract_DMNMag_Diaz_mice('$subj',3)"
#matlab -nodisplay -r "Extract_DMNMag_Diaz_mice('$subj',4)"
#matlab -nodisplay -r "Extract_DMNMag_Diaz_mice('$subj',5)"
#matlab -nodisplay -r "Extract_DMNMag_Diaz_mice('$subj',6)"
