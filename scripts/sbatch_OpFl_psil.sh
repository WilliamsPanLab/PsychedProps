#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=8:00:00 # need >1 (try 10-12 to be safe) if running opfl
#SBATCH -n 4 # try 4 if running opfl
#SBATCH --mem=16G # up to 25 if running opfl
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
#matlab -nodisplay -r "RS_mask_psil('$subj','$sesh')"

# Downsample the data 
#/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/DS_surf_ts_psil.sh $1 $2
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


#############################
#### module III: Calc. Angles
#############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module III: Angular distance calculation"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"

# extract relative angles
#matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs3')"
#matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs4')"
#matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs5')"
#matlab -nodisplay -r "Extract_RelativeAngles_psil('$subj','$sesh','rs6')"

# and magnitudes
#matlab -nodisplay -r "Extract_DMNMag_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNMag_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNMag_psil('$subj','$sesh','rs3')"
#matlab -nodisplay -r "Extract_DMNMag_psil('$subj','$sesh','rs4')"
#matlab -nodisplay -r "Extract_DMNMag_psil('$subj','$sesh','rs5')"
#matlab -nodisplay -r "Extract_DMNMag_psil('$subj','$sesh','rs6')"

# DMN seg
#matlab -nodisplay -r "Extract_DMNSeg_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNSeg_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNSeg_psil('$subj','$sesh','rs3')"
#matlab -nodisplay -r "Extract_DMNSeg_psil('$subj','$sesh','rs4')"
#matlab -nodisplay -r "Extract_DMNSeg_psil('$subj','$sesh','rs5')"
#matlab -nodisplay -r "Extract_DMNSeg_psil('$subj','$sesh','rs6')"


# autocor
#matlab -nodisplay -r "Extract_TAutoCor_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_psil('$subj','$sesh','rs3')"
#matlab -nodisplay -r "Extract_TAutoCor_psil('$subj','$sesh','rs4')"
#matlab -nodisplay -r "Extract_TAutoCor_psil('$subj','$sesh','rs5')"
#matlab -nodisplay -r "Extract_TAutoCor_psil('$subj','$sesh','rs6')"

#############################
#### module IV: Calc Spun Metrics
#############################

# extract magnitudes
#matlab -nodisplay -r "Extract_DMNMag_Spun_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNMag_Spun_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNMag_Spun_psil('$subj','$sesh','rs3')"
#matlab -nodisplay -r "Extract_DMNMag_Spun_psil('$subj','$sesh','rs4')"
#matlab -nodisplay -r "Extract_DMNMag_Spun_psil('$subj','$sesh','rs5')"
#matlab -nodisplay -r "Extract_DMNMag_Spun_psil('$subj','$sesh','rs6')"

echo spunMags done

# extract relative angles
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_psil('$subj','$sesh','rs3')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_psil('$subj','$sesh','rs4')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_psil('$subj','$sesh','rs5')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_psil('$subj','$sesh','rs6')"

echo spunAngles done

# extract DMN FC
#matlab -nodisplay -r "Extract_DMNSeg_Spun_psil('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNSeg_Spun_psil('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNSeg_Spun_psil('$subj','$sesh','rs3')"
#matlab -nodisplay -r "Extract_DMNSeg_Spun_psil('$subj','$sesh','rs4')"
#matlab -nodisplay -r "Extract_DMNSeg_Spun_psil('$subj','$sesh','rs5')"
#matlab -nodisplay -r "Extract_DMNSeg_Spun_psil('$subj','$sesh','rs6')"

echo spunFC done

# extract AutoCor
matlab -nodisplay -r "Extract_TAutoCor_Spun_psil('$subj','$sesh','rs1')"
matlab -nodisplay -r "Extract_TAutoCor_Spun_psil('$subj','$sesh','rs2')"
matlab -nodisplay -r "Extract_TAutoCor_Spun_psil('$subj','$sesh','rs3')"
matlab -nodisplay -r "Extract_TAutoCor_Spun_psil('$subj','$sesh','rs4')"
matlab -nodisplay -r "Extract_TAutoCor_Spun_psil('$subj','$sesh','rs5')"
matlab -nodisplay -r "Extract_TAutoCor_Spun_psil('$subj','$sesh','rs6')"

echo spunAutocor done
#################
echo "OpFl complete"
