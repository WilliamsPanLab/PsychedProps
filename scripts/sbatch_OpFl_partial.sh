#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=48:00:00
#SBATCH -n 1
#SBATCH --mem=13G
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

# downsample baseline for NMF
#/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/DS_surf_ts_mdma_fs5_concat.sh $subj

# SNR Masks: create precursors with fslmaths
# /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/antimask_subj.sh $subj
# SNR Masks: calculate TSNR
# python /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/TSNR_mask_1_ExtractSNR_subjectwise.py $subj

# downsample baseline to fsaverage5 for NMF
#/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/DS_surf_ts_mdma_fs5_rs.sh $subj ses-00

# downsample NMF output to fs4 
#/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/DS_surf_Networks_fs5tofs4_ind.sh $subj

# smooth subject's fs4 nets
#/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/Smooth_DownSampled_Nets.sh $subj

# motion masking
#matlab -nodisplay -r "MotMask('$subj','$sesh','rs1')"
#matlab -nodisplay -r "MotMask('$subj','$sesh','rs2')"
#matlab -nodisplay -r "MotMask('$subj','$sesh','emotion')"
#matlab -nodisplay -r "MotMask('$subj','$sesh','gambling')"
#matlab -nodisplay -r "MotMask('$subj','$sesh','wm')"


#############################
#### module III: Calc. Angles
#############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module III: Angular distance calculation"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
# make output directory outside scratch
# mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj} 

# extract relative angles
matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','rs1')"
matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','rs2')"
matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','wm')"
matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','gambling')"

# extract relative angles
#matlab -nodisplay -r "Extract_RelativeAngles_verts('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles_verts('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles_verts('$subj','$sesh','wm')"
#matlab -nodisplay -r "Extract_RelativeAngles_verts('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_RelativeAngles_verts('$subj','$sesh','emotion')"

# Extract DMN magnitudes
#matlab -nodisplay -r "Extract_DMNMag('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNMag('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNMag('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_DMNMag('$subj','$sesh','wm')"
#matlab -nodisplay -r "Extract_DMNMag('$subj','$sesh','emotion')"

# calculate distance maps
#matlab -nodisplay -r "Calc_AvgBup('$subj')"

# extract autocorr
#matlab -nodisplay -r "Extract_TAutoCor('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_TAutoCor('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_TAutoCor('$subj','$sesh','wm')"

# extract entropy
#matlab -nodisplay -r "Extract_NGSC('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_NGSC('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_NGSC('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_NGSC('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_NGSC('$subj','$sesh','wm')"

# extract DMNSeg
#matlab -nodisplay -r "Extract_DMNSeg('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNSeg('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNSeg('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_DMNSeg('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_DMNSeg('$subj','$sesh','wm')"

# extract DMN win
matlab -nodisplay -r "Extract_DMNWin('$subj','$sesh','rs1')"
matlab -nodisplay -r "Extract_DMNWin('$subj','$sesh','rs2')"
matlab -nodisplay -r "Extract_DMNWin('$subj','$sesh','emotion')"
matlab -nodisplay -r "Extract_DMNWin('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_DMNWin('$subj','$sesh','wm')"

# Streamlines L
#matlab -nodisplay -r "OpFlStreamlines_Left('$subj','$sesh','rs1')"
#matlab -nodisplay -r "OpFlStreamlines_Left('$subj','$sesh','rs2')"
#matlab -nodisplay -r "OpFlStreamlines_Left('$subj','$sesh','emotion')"
#matlab -nodisplay -r "OpFlStreamlines_Left('$subj','$sesh','gambling')"
#matlab -nodisplay -r "OpFlStreamlines_Left('$subj','$sesh','wm')"

# Streamlines R
#matlab -nodisplay -r "OpFlStreamlines_Right('$subj','$sesh','rs1')"
#matlab -nodisplay -r "OpFlStreamlines_Right('$subj','$sesh','rs2')"
#matlab -nodisplay -r "OpFlStreamlines_Right('$subj','$sesh','emotion')"
#matlab -nodisplay -r "OpFlStreamlines_Right('$subj','$sesh','gambling')"
#matlab -nodisplay -r "OpFlStreamlines_Right('$subj','$sesh','wm')"

# DS Streams
#matlab -nodisplay -r "DS_streams('$subj','$sesh','rs1')"
#matlab -nodisplay -r "DS_streams('$subj','$sesh','rs2')"
#matlab -nodisplay -r "DS_streams('$subj','$sesh','emotion')"
#matlab -nodisplay -r "DS_streams('$subj','$sesh','gambling')"
#matlab -nodisplay -r "DS_streams('$subj','$sesh','wm')"

# make 2d histograms of DMN angle/magnitudes
#python3 Viz_AngMag.py $subj $sesh rs1
#python3 Viz_AngMag.py $subj $sesh rs2
#python3 Viz_AngMag.py $subj $sesh emotion
#python3 Viz_AngMag.py $subj $sesh gambling
#python3 Viz_AngMag.py $subj $sesh wm

#################
echo "OpFl complete"
