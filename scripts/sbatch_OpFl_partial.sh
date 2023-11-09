#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=1:00:00
#SBATCH -n 1
#SBATCH --mem=20G
#SBATCH -p leanew1  # Queue names you can submit to
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
ml biology
ml workbench

# subject name is input argument
subj=$1
# sesh is input 2
sesh=$2
# scripts directory
scripts=/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts

#############################
#### module III: Calc. Angles
#############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module III: Angular distance calculation"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"

# make output directory outside scratch
mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj} 

# list preprocessing out directory
xcpdir=/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/${subj}/${sesh}/func

# concatenate resting-state ciftis
${scripts}/Concat_rs_ciftis.sh $xcpdir

# calculate dorsomedial thalamus fc across the cortex
matlab -nodisplay -r "DMThalFC('$subj','$sesh')"

# downsample the dorsomedial thalamus FC maps
${scripts}/DS_surf_DMThalFC.sh ${subj} ${sesh}

# extract relative angles
#matlab -nodisplay -r "Extract_RelativeAngles_fs4('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles_fs4('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles_fs4('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_RelativeAngles_fs4('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_RelativeAngles_fs4('$subj','$sesh','wm')"

# Calculate restting-state FC propagation change correspondence
matlab -nodisplay -r "WithinSubj_DeltafMR('$subj','$sesh','wm')"

#################
echo "OpFl complete"
