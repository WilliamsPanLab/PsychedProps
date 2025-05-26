#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=1:00:00
#SBATCH -n 1
#SBATCH --mem=25G
#SBATCH -p normal,leanew1  # Queue names you can submit to
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


############################
#### module I: Preprocessing
############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module I: preproc"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
# will need matlab
module load matlab
# and freesurfer
module load biology
module load freesurfer/7.3.2
# and workbench
module load workbench
module load contribs poldrack anaconda/5.0.0-py36
# subject name is input argument
subj=$1
# sesh is input 2
sesh=$2

# mask resting-state out from aggregate cifti
matlab -nodisplay -r "RS_mask_psil('$subj','$sesh')"

# Downsample the data 
/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/DS_surf_ts_psil.sh $1 $2
#sleep 20

# cd to workaround addpath in matlab shell call
cd /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts

# mask images: 6+ continuous frames only
matlab -nodisplay -r "MotMask_psil('$subj','$sesh','rs1')"
matlab -nodisplay -r "MotMask_psil('$subj','$sesh','rs2')"
matlab -nodisplay -r "MotMask_psil('$subj','$sesh','rs3')"
matlab -nodisplay -r "MotMask_psil('$subj','$sesh','rs4')"
matlab -nodisplay -r "MotMask_psil('$subj','$sesh','rs5')"
matlab -nodisplay -r "MotMask_psil('$subj','$sesh','rs6')"

############################
#### module II: Optical Flow
############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module II: Optical Flow"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
# Calculate Optical Flow
#matlab -nodisplay -r "OpFl_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "OpFl_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "OpFl_psil('$subj','$sesh','rs3')"
#matlab -nodisplay -r "OpFl_psil('$subj','$sesh','rs4')"
#matlab -nodisplay -r "OpFl_psil('$subj','$sesh','rs5')"
#matlab -nodisplay -r "OpFl_psil('$subj','$sesh','rs6')"

echo "Subcortical runs"
#matlab -nodisplay -r "-r SubCort_OpFl_caud_R_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "-r SubCort_OpFl_caud_R_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "-r SubCort_OpFl_caud_R_psil('$subj','$sesh','emotion')"
#matlab -nodisplay -r "-r SubCort_OpFl_caud_R_psil('$subj','$sesh','gambling')"
#matlab -nodisplay -r "-r SubCort_OpFl_caud_R_psil('$subj','$sesh','wm')"

#matlab -nodisplay -r "-r SubCort_OpFl_caud_L_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "-r SubCort_OpFl_caud_L_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "-r SubCort_OpFl_caud_L_psil('$subj','$sesh','emotion')"
#matlab -nodisplay -r "-r SubCort_OpFl_caud_L_psil('$subj','$sesh','gambling')"
#matlab -nodisplay -r "-r SubCort_OpFl_caud_L_psil('$subj','$sesh','wm')"

#matlab -nodisplay -r "-r SubCort_OpFl_hippo_R_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "-r SubCort_OpFl_hippo_R_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "-r SubCort_OpFl_hippo_R_psil('$subj','$sesh','emotion')"
#matlab -nodisplay -r "-r SubCort_OpFl_hippo_R_psil('$subj','$sesh','gambling')"
#matlab -nodisplay -r "-r SubCort_OpFl_hippo_R_psil('$subj','$sesh','wm')"

#matlab -nodisplay -r "-r SubCort_OpFl_hippo_L_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "-r SubCort_OpFl_hippo_L_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "-r SubCort_OpFl_hippo_L_psil('$subj','$sesh','emotion')"
#matlab -nodisplay -r "-r SubCort_OpFl_hippo_L_psil('$subj','$sesh','gambling')"
#matlab -nodisplay -r "-r SubCort_OpFl_hippo_L_psil('$subj','$sesh','wm')"

# RS filepaths
childfp=/scratch/users/apines/data/mdma/${subj}/${sesh}

#############################
#### module III: Calc. Angles
#############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module III: Angular distance calculation"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
# make output directory outside scratch
#mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj} 

# extract relative angles
#matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs3')"
#matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs4')"
#matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs5')"
#matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs6')"

# extract DMN FC
#matlab -nodisplay -r "Extract_DMNSeg_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNSeg_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNSeg_psil('$subj','$sesh','rs3')"
#matlab -nodisplay -r "Extract_DMNSeg_psil('$subj','$sesh','rs4')"
#matlab -nodisplay -r "Extract_DMNSeg_psil('$subj','$sesh','rs5')"
#matlab -nodisplay -r "Extract_DMNSeg_psil('$subj','$sesh','rs6')"

# extract DMN temporal autocor
#matlab -nodisplay -r "Extract_TAutoCor_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_psil('$subj','$sesh','rs3')"
#matlab -nodisplay -r "Extract_TAutoCor_psil('$subj','$sesh','rs4')"
#matlab -nodisplay -r "Extract_TAutoCor_psil('$subj','$sesh','rs5')"
#matlab -nodisplay -r "Extract_TAutoCor_psil('$subj','$sesh','rs6')"

# extract vertexwise angles
#matlab -nodisplay -r "Extract_RelativeAngles_verts_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles_verts_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles_verts_psil('$subj','$sesh','rs3')"
#matlab -nodisplay -r "Extract_RelativeAngles_verts_psil('$subj','$sesh','rs4')"
#matlab -nodisplay -r "Extract_RelativeAngles_verts_psil('$subj','$sesh','rs5')"
#matlab -nodisplay -r "Extract_RelativeAngles_verts_psil('$subj','$sesh','rs6')"

# Extract DMN mag
#matlab -nodisplay -r "Extract_DMNMag_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNMag_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNMag_psil('$subj','$sesh','rs3')"
#matlab -nodisplay -r "Extract_DMNMag_psil('$subj','$sesh','rs4')"
#matlab -nodisplay -r "Extract_DMNMag_psil('$subj','$sesh','rs5')"
#matlab -nodisplay -r "Extract_DMNMag_psil('$subj','$sesh','rs6')"


#################
echo "OpFl complete"
