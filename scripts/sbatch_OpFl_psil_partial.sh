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
# matlab -nodisplay -r "RS_mask_psil('$subj','$sesh')"

# Downsample the data 
# /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/DS_surf_ts_psil.sh $1 $2
#sleep 20

# cd to workaround addpath in matlab shell call
cd /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts

# mask images: 6+ continuous frames only
#matlab -nodisplay -r "MotMask_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "MotMask_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "MotMask_psil('$subj','$sesh','rs3')"
#matlab -nodisplay -r "MotMask_psil('$subj','$sesh','rs4')"
#matlab -nodisplay -r "MotMask_psil('$subj','$sesh','rs5')"
#matlab -nodisplay -r "MotMask_psil('$subj','$sesh','rs6')"

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

# RS filepaths
childfp=/scratch/users/apines/data/mdma/${subj}/${sesh}

# interpolate fs4 time series to faces and between-timepoints
#matlab -nodisplay -r "InterpolateTS_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "InterpolateTS_psil('$subj','$sesh','rs2')"

#############################
#### module III: Calc. Angles
#############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module III: Angular distance calculation"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
# make output directory outside scratch
#mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj} 

# extract relative angles
matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs1')"
matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs2')"
matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs3')"
matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs4')"
matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs5')"
matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs6')"

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

# combine angular time series with magnitude time series
#matlab -nodisplay -r "Combine_FacewiseTS('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Combine_FacewiseTS('$subj','$sesh','rs2')"

# OPFL TO VERTS
#matlab -nodisplay -r "OpFl_toVerts_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "OpFl_toVerts_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "OpFl_toVerts_psil('$subj','$sesh','rs3')"
#matlab -nodisplay -r "OpFl_toVerts_psil('$subj','$sesh','rs4')"
#matlab -nodisplay -r "OpFl_toVerts_psil('$subj','$sesh','rs5')"
#matlab -nodisplay -r "OpFl_toVerts_psil('$subj','$sesh','rs6')"

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


### OUTSIDE OF SUBJECT LEVEL

# PSIL MATS TO CSV FOR R

# R SCRIPT FOR VERTEXWISE

#############################
#### module IV: Create figures
#############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module IV: Creating figures"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"

### filepaths for angular carpetplots
#ATSl=[childFP '/' subj '_' sesh '_Prop_TS_dmn_L.csv'];
#ATSr=[childFP '/' subj '_' sesh '_Prop_TS_dmn_L.csv'];

### make figure directory outside of scratch
#mkdir /oak/stanford/groups/leanew1/users/apines/data/p50/${subj}/${sesh}/figs

### create interpolated carpetplot of bold
#python3 Viz_ITS.py $subj $sesh

### create angular carpetplot 
#python3 Viz_ATS.py $subj $sesh

### create positivity plots
#python3 Viz_AngMag.py $subj $sesh

### copy xcpd motion plots into figs dir
#cp /scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/${subj}/figures/${subj}_${sesh}_task-rs_acq-mb_dir-pe?_run-0_space-fsLR_desc-censoring_motion.svg /oak/stanford/groups/leanew1/users/apines/data/p50/${subj}/${sesh}/figs

#################
echo "OpFl complete"
