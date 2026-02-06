#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=90:00:00
#SBATCH -c 4
#SBATCH --mem=16G
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
# Downsample the data 
#/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/DS_surf_ts_mdma.sh $1 $2
#sleep 20

# cd to workaround addpath in matlab shell call
cd /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts

# mask images: 8+ continuous frames only
#matlab -nodisplay -r "MotMask('$subj','$sesh','rs1')"
#matlab -nodisplay -r "MotMask('$subj','$sesh','rs2')"
#matlab -nodisplay -r "MotMask('$subj','$sesh','emotion')"
#matlab -nodisplay -r "MotMask('$subj','$sesh','gambling')"
#matlab -nodisplay -r "MotMask('$subj','$sesh','wm')"

############################
#### module II: Optical Flow
############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module II: Optical Flow"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
# Calculate Optical Flow
#matlab -nodisplay -r "OpFl_mdma('$subj','$sesh','rs1')"
#matlab -nodisplay -r "OpFl_mdma('$subj','$sesh','rs2')"
#matlab -nodisplay -r "OpFl_mdma('$subj','$sesh','emotion')"
#matlab -nodisplay -r "OpFl_mdma('$subj','$sesh','gambling')"
#matlab -nodisplay -r "OpFl_mdma('$subj','$sesh','wm')"

#############################
#### module III: Calc. DMN Metrics
#############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module III: Angular distance calculation"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
# make output directory outside scratch
#mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj} 

# extract magnitudes
#matlab -nodisplay -r "Extract_DMNMag('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNMag('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNMag('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_DMNMag('$subj','$sesh','wm')"
#matlab -nodisplay -r "Extract_DMNMag('$subj','$sesh','emotion')"

# extract relative angles
#matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','wm')"

# extract DMN FC
#matlab -nodisplay -r "Extract_DMNSeg('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNSeg('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNSeg('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_DMNSeg('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_DMNSeg('$subj','$sesh','wm')"

# extract AutoCor
#matlab -nodisplay -r "Extract_TAutoCor('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_TAutoCor('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_TAutoCor('$subj','$sesh','wm')"

# for revisions
# mag
#matlab -nodisplay -r "Extract_VISMag('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_VISMag('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_VISMag('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_VISMag('$subj','$sesh','wm')"
#matlab -nodisplay -r "Extract_VISMag('$subj','$sesh','emotion')"

#matlab -nodisplay -r "Extract_MOTMag('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_MOTMag('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_MOTMag('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_MOTMag('$subj','$sesh','wm')"
#matlab -nodisplay -r "Extract_MOTMag('$subj','$sesh','emotion')"

#matlab -nodisplay -r "Extract_FPNMag('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_FPNMag('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_FPNMag('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_FPNMag('$subj','$sesh','wm')"
#matlab -nodisplay -r "Extract_FPNMag('$subj','$sesh','emotion')"


# seg
#matlab -nodisplay -r "Extract_VISSeg('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_VISSeg('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_VISSeg('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_VISSeg('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_VISSeg('$subj','$sesh','wm')"

#matlab -nodisplay -r "Extract_MOTSeg('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_MOTSeg('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_MOTSeg('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_MOTSeg('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_MOTSeg('$subj','$sesh','wm')"

#matlab -nodisplay -r "Extract_FPNSeg('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_FPNSeg('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_FPNSeg('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_FPNSeg('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_FPNSeg('$subj','$sesh','wm')"

# ta
#matlab -nodisplay -r "Extract_TAutoCor_VIS('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_VIS('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_VIS('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_TAutoCor_VIS('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_TAutoCor_VIS('$subj','$sesh','wm')"

#matlab -nodisplay -r "Extract_TAutoCor_MOT('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_MOT('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_MOT('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_TAutoCor_MOT('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_TAutoCor_MOT('$subj','$sesh','wm')"

#matlab -nodisplay -r "Extract_TAutoCor_FPN('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_FPN('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_FPN('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_TAutoCor_FPN('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_TAutoCor_FPN('$subj','$sesh','wm')"

echo #############################
echo #### module IV: Calc Spun Metrics
echo #############################

# extract magnitudes
#matlab -nodisplay -r "Extract_DMNMag_Spun('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNMag_Spun('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNMag_Spun('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_DMNMag_Spun('$subj','$sesh','wm')"
#matlab -nodisplay -r "Extract_DMNMag_Spun('$subj','$sesh','emotion')"

echo spunMags done

# extract relative angles
#matlab -nodisplay -r "Extract_RelativeAngles_Spun('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun('$subj','$sesh','wm')"

echo spunAngles done

# extract DMN FC
#matlab -nodisplay -r "Extract_DMNSeg_Spun('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_DMNSeg_Spun('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_DMNSeg_Spun('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_DMNSeg_Spun('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_DMNSeg_Spun('$subj','$sesh','wm')"

echo spunFC done

# extract AutoCor
#matlab -nodisplay -r "Extract_TAutoCor_Spun('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_Spun('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_Spun('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_TAutoCor_Spun('$subj','$sesh','gambling')"
#matlab -nodisplay -r "Extract_TAutoCor_Spun('$subj','$sesh','wm')"

echo spunAutocor done

##################################################
##################################################
##################################################
### ADDING OTHER NETS FOR REVISIONS
##################################################
##################################################
##################################################

echo #############################
echo #### module IV: Calc Spun Metrics - MOT
echo #############################

# extract magnitudes
#matlab -nodisplay -r "Extract_MOTMag_Spun('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_MOTMag_Spun('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_MOTMag_Spun('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_MOTMag_Spun('$subj','$sesh','wm')"
#matlab -nodisplay -r "Extract_MOTMag_Spun('$subj','$sesh','emotion')"

echo motor spunMags done

# extract relative angles
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_MOT('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_MOT('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_MOT('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_MOT('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_RelativeAngles_Spun_MOT('$subj','$sesh','wm')"

echo motor spunAngles done

# extract DMN FC
#matlab -nodisplay -r "Extract_MOTSeg_Spun('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_MOTSeg_Spun('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_MOTSeg_Spun('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_MOTSeg_Spun('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_MOTSeg_Spun('$subj','$sesh','wm')"

echo spunFC done

# extract AutoCor
#matlab -nodisplay -r "Extract_TAutoCor_MOT_Spun('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_MOT_Spun('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_MOT_Spun('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_TAutoCor_MOT_Spun('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_TAutoCor_MOT_Spun('$subj','$sesh','wm')"

echo spunAutocor done


echo #############################
echo #### module IV: Calc Spun Metrics - VIS
echo #############################

# extract magnitudes
#matlab -nodisplay -r "Extract_VISMag_Spun('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_VISMag_Spun('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_VISMag_Spun('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_VISMag_Spun('$subj','$sesh','wm')"
#matlab -nodisplay -r "Extract_VISMag_Spun('$subj','$sesh','emotion')"

echo VIS spunMags done

# extract relative angles
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_VIS('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_VIS('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_VIS('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_VIS('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_RelativeAngles_Spun_VIS('$subj','$sesh','wm')"

echo visual spunAngles done

# extract DMN FC
#matlab -nodisplay -r "Extract_VISSeg_Spun('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_VISSeg_Spun('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_VISSeg_Spun('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_VISSeg_Spun('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_VISSeg_Spun('$subj','$sesh','wm')"

echo spunFC done

# extract AutoCor
#matlab -nodisplay -r "Extract_TAutoCor_VIS_Spun('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_VIS_Spun('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_VIS_Spun('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_TAutoCor_VIS_Spun('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_TAutoCor_VIS_Spun('$subj','$sesh','wm')"

echo spunAutocor done




echo #############################
echo #### module IV: Calc Spun Metrics - FPN
echo #############################

# extract magnitudes
#matlab -nodisplay -r "Extract_FPNMag_Spun('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_FPNMag_Spun('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_FPNMag_Spun('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_FPNMag_Spun('$subj','$sesh','wm')"
#matlab -nodisplay -r "Extract_FPNMag_Spun('$subj','$sesh','emotion')"

echo FPN spunMags done

# extract relative angles
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_FPN('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_FPN('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_FPN('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_RelativeAngles_Spun_FPN('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_RelativeAngles_Spun_FPN('$subj','$sesh','wm')"

echo FPN spunAngles done

# extract DMN FC
#matlab -nodisplay -r "Extract_FPNSeg_Spun('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_FPNSeg_Spun('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_FPNSeg_Spun('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_FPNSeg_Spun('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_FPNSeg_Spun('$subj','$sesh','wm')"

echo spunFC done

# extract AutoCor
#matlab -nodisplay -r "Extract_TAutoCor_FPN_Spun('$subj','$sesh','rs1')"
#matlab -nodisplay -r "Extract_TAutoCor_FPN_Spun('$subj','$sesh','rs2')"
#matlab -nodisplay -r "Extract_TAutoCor_FPN_Spun('$subj','$sesh','emotion')"
#matlab -nodisplay -r "Extract_TAutoCor_FPN_Spun('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_TAutoCor_FPN_Spun('$subj','$sesh','wm')"

echo spunAutocor done























#################
echo "OpFl complete"
