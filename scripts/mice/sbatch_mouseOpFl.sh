#!/bin/bash
#
#SBATCH --job-name=MouseOF
#SBATCH --time=20:00:00
#SBATCH -n 1
#SBATCH --mem=27G
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

# note that dex and diaz have different mouse IDs (for the most part): this is not intended to run all commented-out scripts at once on single mice
# instead, remove #'s from only the scripts you need: you'll want to process one drug cohort at a time
# also note that not all drugs have 6 acquisitions normatively (in addition to some limited within-drug-cohort variability in acquisitions)

# run optical flow (LSD)
#matlab -nodisplay -r "Mouse_OpFl('$subj',1)"
#matlab -nodisplay -r "Mouse_OpFl('$subj',2)"
#matlab -nodisplay -r "Mouse_OpFl('$subj',3)"
#matlab -nodisplay -r "Mouse_OpFl('$subj',4)"
#matlab -nodisplay -r "Mouse_OpFl('$subj',5)"
#matlab -nodisplay -r "Mouse_OpFl('$subj',6)"

# run OF for dex
#matlab -nodisplay -r "Mouse_OpFl_Dex('$subj',1)"
#matlab -nodisplay -r "Mouse_OpFl_Dex('$subj',2)"
#matlab -nodisplay -r "Mouse_OpFl_Dex('$subj',3)"
#matlab -nodisplay -r "Mouse_OpFl_Dex('$subj',4)"
#matlab -nodisplay -r "Mouse_OpFl_Dex('$subj',5)"
#matlab -nodisplay -r "Mouse_OpFl_Dex('$subj',6)"

# run of for diaz
#matlab -nodisplay -r "Mouse_OpFl_Diaz('$subj',1)"
#matlab -nodisplay -r "Mouse_OpFl_Diaz('$subj',2)"
#matlab -nodisplay -r "Mouse_OpFl_Diaz('$subj',3)"
#matlab -nodisplay -r "Mouse_OpFl_Diaz('$subj',4)"
#matlab -nodisplay -r "Mouse_OpFl_Diaz('$subj',5)"
#matlab -nodisplay -r "Mouse_OpFl_Diaz('$subj',6)"

# extract DMN magnitudes (LSD)
#matlab -nodisplay -r "Extract_DMNMag_mice('$subj',1)"
#matlab -nodisplay -r "Extract_DMNMag_mice('$subj',2)"
#matlab -nodisplay -r "Extract_DMNMag_mice('$subj',3)"
#matlab -nodisplay -r "Extract_DMNMag_mice('$subj',4)"
#matlab -nodisplay -r "Extract_DMNMag_mice('$subj',5)"
#matlab -nodisplay -r "Extract_DMNMag_mice('$subj',6)"

# extract DMN agnles (LSD)
matlab -nodisplay -r "Extract_RelativeAngles_mice('$subj',1)"
matlab -nodisplay -r "Extract_RelativeAngles_mice('$subj',2)"
matlab -nodisplay -r "Extract_RelativeAngles_mice('$subj',3)"
matlab -nodisplay -r "Extract_RelativeAngles_mice('$subj',4)"
matlab -nodisplay -r "Extract_RelativeAngles_mice('$subj',5)"
matlab -nodisplay -r "Extract_RelativeAngles_mice('$subj',6)"

# extract DMN segregation (LSD)
#matlab -nodisplay -r "Extract_DMNSeg_mice('$subj',1)"
#matlab -nodisplay -r "Extract_DMNSeg_mice('$subj',2)"
#matlab -nodisplay -r "Extract_DMNSeg_mice('$subj',3)"
#matlab -nodisplay -r "Extract_DMNSeg_mice('$subj',4)"
#matlab -nodisplay -r "Extract_DMNSeg_mice('$subj',5)"
#matlab -nodisplay -r "Extract_DMNSeg_mice('$subj',6)"

#dex mag
#matlab -nodisplay -r "Extract_DMNMag_Dex_mice('$subj',1)"
#matlab -nodisplay -r "Extract_DMNMag_Dex_mice('$subj',2)"
#matlab -nodisplay -r "Extract_DMNMag_Dex_mice('$subj',3)"
#matlab -nodisplay -r "Extract_DMNMag_Dex_mice('$subj',4)"
#matlab -nodisplay -r "Extract_DMNMag_Dex_mice('$subj',5)"
#matlab -nodisplay -r "Extract_DMNMag_Dex_mice('$subj',6)"

# DMN autocorrelation (LSD)
#matlab -nodisplay -r "Extract_TAutoCor_mice('$subj',1)"
#matlab -nodisplay -r "Extract_TAutoCor_mice('$subj',2)"
#matlab -nodisplay -r "Extract_TAutoCor_mice('$subj',3)"
#matlab -nodisplay -r "Extract_TAutoCor_mice('$subj',4)"
#matlab -nodisplay -r "Extract_TAutoCor_mice('$subj',5)"
#matlab -nodisplay -r "Extract_TAutoCor_mice('$subj',6)"

# diaz mag
#matlab -nodisplay -r "Extract_DMNMag_Diaz_mice('$subj',1)"
#matlab -nodisplay -r "Extract_DMNMag_Diaz_mice('$subj',2)"
#matlab -nodisplay -r "Extract_DMNMag_Diaz_mice('$subj',3)"
#matlab -nodisplay -r "Extract_DMNMag_Diaz_mice('$subj',4)"
#matlab -nodisplay -r "Extract_DMNMag_Diaz_mice('$subj',5)"
#matlab -nodisplay -r "Extract_DMNMag_Diaz_mice('$subj',6)"


echo spun metrics
#matlab -nodisplay -r "Extract_DMNMag_Spun_mice('$subj',1)"
#matlab -nodisplay -r "Extract_DMNMag_Spun_mice('$subj',2)"
#matlab -nodisplay -r "Extract_DMNMag_Spun_mice('$subj',3)"
#matlab -nodisplay -r "Extract_DMNMag_Spun_mice('$subj',4)"
#matlab -nodisplay -r "Extract_DMNMag_Spun_mice('$subj',5)"
#matlab -nodisplay -r "Extract_DMNMag_Spun_mice('$subj',6)"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_mice('$subj',1)"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_mice('$subj',2)"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_mice('$subj',3)"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_mice('$subj',4)"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_mice('$subj',5)"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_mice('$subj',6)"
#matlab -nodisplay -r "Extract_DMNSeg_mice_Spun('$subj',1)"
#matlab -nodisplay -r "Extract_DMNSeg_mice_Spun('$subj',2)"
#matlab -nodisplay -r "Extract_DMNSeg_mice_Spun('$subj',3)"
#matlab -nodisplay -r "Extract_DMNSeg_mice_Spun('$subj',4)"
#matlab -nodisplay -r "Extract_DMNSeg_mice_Spun('$subj',5)"
#matlab -nodisplay -r "Extract_DMNSeg_mice_Spun('$subj',6)"
#matlab -nodisplay -r "Extract_TAutoCor_mice_Spun('$subj',1)"
#matlab -nodisplay -r "Extract_TAutoCor_mice_Spun('$subj',2)"
#matlab -nodisplay -r "Extract_TAutoCor_mice_Spun('$subj',3)"
#matlab -nodisplay -r "Extract_TAutoCor_mice_Spun('$subj',4)"
#matlab -nodisplay -r "Extract_TAutoCor_mice_Spun('$subj',5)"
#matlab -nodisplay -r "Extract_TAutoCor_mice_Spun('$subj',6)"
