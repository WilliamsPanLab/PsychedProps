#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=1:00:00
#SBATCH -n 1
#SBATCH --mem=20G
#SBATCH -p leanew1,normal # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------


############    Buffaflow   ###################
##             _.-````'-,,
##   _,.,_ ,-'`           ``''-.,
## /)     (\                   '``-.
##((      ) )                ->    `\
## \)    (_/               ->       )\
##  |       /)           '    ,'    / \
##  `\    ^'            '     (    /  ))
##    |      _/\ ,     /    ,,`\   (  "`
##     \Y,   |  \  \  /____| / \_ \##
# This script is to prep the environment of a slurm node, 
##       `)_/    \  \  )    (->  (->
# Launch optical flow, caclulate angular distances, generate figures
##                \( \(     |/   |/
# and delete interim files upon completion
##    mic & dwb  /_(/_(    /_(  /_(
ml matlab
# and freesurfer
module load biology
module load fsl
module load freesurfer/7.3.2
# and workbench
module load workbench
module load contribs poldrack anaconda/5.0.0-py36
module load python/3.9
# subject name is input argument
subj=$1
# sesh is input 2
sesh=$2

# SNR Masks: create precursors with fslmaths
# /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/antimask_subj.sh $subj
# SNR Masks: calculate TSNR
# python /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/TSNR_mask_1_ExtractSNR_subjectwise.py $subj

# downsample baseline to fsaverage5 for NMF
/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/DS_surf_ts_mdma_fs5_rs.sh $subj ses-00

# interpolate fs4 time series to faces and between-timepoints
#matlab -nodisplay -r "InterpolateTS('$subj','$sesh','rs1')"
#matlab -nodisplay -r "InterpolateTS('$subj','$sesh','rs2')"
#matlab -nodisplay -r "InterpolateTS('$subj','$sesh','emotion')"
#matlab -nodisplay -r "InterpolateTS('$subj','$sesh','gambling')"
#matlab -nodisplay -r "InterpolateTS('$subj','$sesh','wm')"

# combine angular time series with magnitude time series
#matlab -nodisplay -r "Combine_FacewiseTS('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Combine_FacewiseTS('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Combine_FacewiseTS('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Combine_FacewiseTS('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Combine_FacewiseTS('$subj','$sesh','wm')"

#############################
#### module III: Calc. Angles
#############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module III: Angular distance calculation"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
# make output directory outside scratch
# mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj} 

# extract relative angles
# group
#matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','wm')"

# extract autocorr
#matlab -nodisplay -r "Extract_AutoCor('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_AutoCor('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_AutoCor('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_AutoCor('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_AutoCor('$subj','$sesh','wm')"

# extract entropy
#matlab -nodisplay -r "Extract_NGSC('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_NGSC('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_NGSC('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_NGSC('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_NGSC('$subj','$sesh','wm')"

# extract amygdalar FC (loops over tasks internally)
#matlab -nodisplay -r "Extract_AmygFC('$subj','$sesh')"

# make 2d histograms of DMN angle/magnitudes
#python3 Viz_AngMag.py $subj $sesh rs1
#python3 Viz_AngMag.py $subj $sesh rs2
#python3 Viz_AngMag.py $subj $sesh emotion
#python3 Viz_AngMag.py $subj $sesh gambling
#python3 Viz_AngMag.py $subj $sesh wm

#################
echo "OpFl complete"
