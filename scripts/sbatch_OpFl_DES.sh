#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=2:30:00
#SBATCH -n 1
#SBATCH --mem=18G
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
# sesh is input 2: only baseline for connectome
sesh=ses-BL

# Downsample the data 
/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/DS_surf_ts_mdma.sh $subj $sesh
sleep 20

# mask images: 8+ continuous frames only
matlab -nodisplay -r "MotMask_DES('$subj','$sesh','rs1')"
matlab -nodisplay -r "MotMask_DES('$subj','$sesh','rs2')"
matlab -nodisplay -r "MotMask_DES('$subj','$sesh','emotion')"
matlab -nodisplay -r "MotMask_DES('$subj','$sesh','gambling')"
matlab -nodisplay -r "MotMask_DES('$subj','$sesh','wm')"

############################
#### module II: Optical Flow
############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module II: Optical Flow"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
# Calculate Optical Flow
matlab -nodisplay -r "OpFl_mdma('$subj','$sesh','rs1')"
matlab -nodisplay -r "OpFl_mdma('$subj','$sesh','rs2')"
matlab -nodisplay -r "OpFl_mdma('$subj','$sesh','emotion')"
matlab -nodisplay -r "OpFl_mdma('$subj','$sesh','gambling')"
matlab -nodisplay -r "OpFl_mdma('$subj','$sesh','wm')"

#############################
#### module III: Calc. Angles
#############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module III: Angular distance calculation"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
# make output directory outside scratch
mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj}

# extract relative angles
matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','rs1')"
matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','rs2')"
matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','emotion')"
matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','wm')"

# extract DMN mag
matlab -nodisplay -r "Extract_DMNMag('$subj','$sesh','rs1')"
matlab -nodisplay -r "Extract_DMNMag('$subj','$sesh','rs2')"
matlab -nodisplay -r "Extract_DMNMag('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_DMNMag('$subj','$sesh','wm')"
matlab -nodisplay -r "Extract_DMNMag('$subj','$sesh','emotion')"

# extract DMN FC
matlab -nodisplay -r "Extract_DMNSeg('$subj','$sesh','rs1')"
matlab -nodisplay -r "Extract_DMNSeg('$subj','$sesh','rs2')"
matlab -nodisplay -r "Extract_DMNSeg('$subj','$sesh','emotion')"
matlab -nodisplay -r "Extract_DMNSeg('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_DMNSeg('$subj','$sesh','wm')"

# extract AutoCor
matlab -nodisplay -r "Extract_TAutoCor('$subj','$sesh','rs1')"
matlab -nodisplay -r "Extract_TAutoCor('$subj','$sesh','rs2')"
matlab -nodisplay -r "Extract_TAutoCor('$subj','$sesh','emotion')"
matlab -nodisplay -r "Extract_TAutoCor('$subj','$sesh','gambling')"
matlab -nodisplay -r "Extract_TAutoCor('$subj','$sesh','wm')"


#################
echo "OpFl complete"

