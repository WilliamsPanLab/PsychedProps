#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=10:00:00 # need >1 (try 10-12 to be safe) if running opfl
#SBATCH -c 4 # try 4 if running opfl
#SBATCH --mem=16G # up to 25 if running opfl
#SBATCH -p normal,owners,anishm,leanew1  # Queue names you can submit to
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
#matlab -nodisplay -r "RS_mask_lsd('$subj','$sesh')"

# BP data
#matlab -nodisplay -r "BP_cifti('$subj','$sesh','rs1')"
#matlab -nodisplay -r "BP_cifti('$subj','$sesh','rs2')"
#matlab -nodisplay -r "BP_cifti('$subj','$sesh','mus')"

# Downsample the data 
#/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/DS_surf_ts_lsd.sh $1 $2

#sleep 20

# cd to workaround addpath in matlab shell call
cd /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts

# mask images: 6+ continuous frames only
#matlab -nodisplay -r "MotMask_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "MotMask_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "MotMask_lsd('$subj','$sesh','mus')"

############################
#### module II: Optical Flow
############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module II: Optical Flow"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
# Calculate Optical Flow
#matlab -nodisplay -r "OpFl_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "OpFl_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "OpFl_lsd('$subj','$sesh','mus')"


#############################
#### module III: Calc. Angles
#############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module III: Angular distance calculation"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"

# extract relative angles
matlab -nodisplay -r "Extract_RelativeAngles_lsd('$subj','$sesh','rs1')"
matlab -nodisplay -r "Extract_RelativeAngles_lsd('$subj','$sesh','rs2')"
matlab -nodisplay -r "Extract_RelativeAngles_lsd('$subj','$sesh','mus')"

# and magnitudes
#matlab -nodisplay -r "Extract_DMNMag_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNMag_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNMag_lsd('$subj','$sesh','mus')"

# new for revision
#matlab -nodisplay -r "Extract_VISMag_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_VISMag_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_VISMag_lsd('$subj','$sesh','mus')"

#matlab -nodisplay -r "Extract_MOTMag_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_MOTMag_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_MOTMag_lsd('$subj','$sesh','mus')"

#matlab -nodisplay -r "Extract_FPNMag_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_FPNMag_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_FPNMag_lsd('$subj','$sesh','mus')"

# DMN seg
#matlab -nodisplay -r "Extract_DMNSeg_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNSeg_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNSeg_lsd('$subj','$sesh','mus')"
#
# new for revision
#matlab -nodisplay -r "Extract_VISSeg_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_VISSeg_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_VISSeg_lsd('$subj','$sesh','mus')"

#matlab -nodisplay -r "Extract_MOTSeg_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_MOTSeg_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_MOTSeg_lsd('$subj','$sesh','mus')"

#matlab -nodisplay -r "Extract_FPNSeg_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_FPNSeg_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_FPNSeg_lsd('$subj','$sesh','mus')"

# autocor
#matlab -nodisplay -r "Extract_TAutoCor_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_lsd('$subj','$sesh','mus')"

# new for revision

#matlab -nodisplay -r "Extract_TAutoCor_lsd_VIS('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_lsd_VIS('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_lsd_VIS('$subj','$sesh','mus')"

#matlab -nodisplay -r "Extract_TAutoCor_lsd_MOT('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_lsd_MOT('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_lsd_MOT('$subj','$sesh','mus')"

#matlab -nodisplay -r "Extract_TAutoCor_lsd_FPN('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_lsd_FPN('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_lsd_FPN('$subj','$sesh','mus')"

#############################
#### module IV: Calc Spun Metrics
#############################

# extract magnitudes
#matlab -nodisplay -r "Extract_DMNMag_Spun_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNMag_Spun_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNMag_Spun_lsd('$subj','$sesh','mus')"

echo spunMags done

# extract relative angles
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_lsd('$subj','$sesh','mus')"

echo spunAngles done

# extract DMN FC
#matlab -nodisplay -r "Extract_DMNSeg_Spun_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNSeg_Spun_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNSeg_Spun_lsd('$subj','$sesh','mus')"

echo spunFC done

# extract AutoCor
#matlab -nodisplay -r "Extract_TAutoCor_Spun_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_Spun_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_Spun_lsd('$subj','$sesh','mus')"

echo spunAutocor done






# REVISIONS CALLS FOR UPDATE TO THIS TO TEST OTHER NETS


##### MOT

# extract magnitudes
#matlab -nodisplay -r "Extract_MOTMag_Spun_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_MOTMag_Spun_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_MOTMag_Spun_lsd('$subj','$sesh','mus')"

echo spunMags done

# extract relative angles
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_lsd_MOT('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_lsd_MOT('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_lsd_MOT('$subj','$sesh','mus')"

echo spunAngles done

# extract DMN FC
#matlab -nodisplay -r "Extract_MOTSeg_Spun_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_MOTSeg_Spun_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_MOTSeg_Spun_lsd('$subj','$sesh','mus')"

echo spunFC done

# extract AutoCor
#matlab -nodisplay -r "Extract_TAutoCor_MOT_Spun_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_MOT_Spun_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_MOT_Spun_lsd('$subj','$sesh','mus')"

echo spunAutocor done


############## VIS

# extract magnitudes
#matlab -nodisplay -r "Extract_VISMag_Spun_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_VISMag_Spun_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_VISMag_Spun_lsd('$subj','$sesh','mus')"

echo spunMags done

# extract relative angles
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_lsd_VIS('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_lsd_VIS('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_lsd_VIS('$subj','$sesh','mus')"

echo spunAngles done

# extract DMN FC
#matlab -nodisplay -r "Extract_VISSeg_Spun_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_VISSeg_Spun_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_VISSeg_Spun_lsd('$subj','$sesh','mus')"

echo spunFC done

# extract AutoCor
#matlab -nodisplay -r "Extract_TAutoCor_VIS_Spun_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_VIS_Spun_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_VIS_Spun_lsd('$subj','$sesh','mus')"

echo spunAutocor done



################# FPN


# extract magnitudes
#matlab -nodisplay -r "Extract_FPNMag_Spun_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_FPNMag_Spun_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_FPNMag_Spun_lsd('$subj','$sesh','mus')"

echo spunMags done

# extract relative angles
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_lsd_FPN('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_lsd_FPN('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_lsd_FPN('$subj','$sesh','mus')"

echo spunAngles done

# extract DMN FC
#matlab -nodisplay -r "Extract_FPNSeg_Spun_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_FPNSeg_Spun_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_FPNSeg_Spun_lsd('$subj','$sesh','mus')"

echo spunFC done

# extract AutoCor
#matlab -nodisplay -r "Extract_TAutoCor_FPN_Spun_lsd('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_FPN_Spun_lsd('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_FPN_Spun_lsd('$subj','$sesh','mus')"

echo spunAutocor done











#################
echo "OpFl complete"
