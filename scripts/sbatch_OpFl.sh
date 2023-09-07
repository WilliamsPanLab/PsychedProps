#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=10:20:00
#SBATCH -n 1
#SBATCH --mem=35G
#SBATCH -p leanew1,normal  # Queue names you can submit to
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
/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/fs_5/DS_surf_ts_mdma_fs5.sh $1 $2

# cd to workaround addpath in matlab shell call
cd /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts

# mask images: 8+ continuous frames only
matlab -nodisplay -r "MotMask('$subj','$sesh')"

############################
#### module II: Optical Flow
############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module II: Optical Flow"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
# Calculate Optical Flow
matlab -nodisplay -r "OpFl_mdma_fs5('$subj','$sesh')"

# RS filepaths
childfp=/scratch/users/apines/data/mdma/${subj}/${sesh}
rsIn=${childfp}/${subj}_${sesh}_OpFl_rs_fs5.mat

# interpolate fs5 time series to faces and between-timepoints
matlab -nodisplay -r "InterpolateTS('$subj','$sesh')"

#############################
#### module III: Calc. Angles
#############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module III: Angular distance calculation"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
# make output directory outside scratch
mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj} 

# extract relative angles
matlab -nodisplay -r "Extract_RelativeAngles('$subj','$sesh','$rsIn')"

# combine angular time series with magnitude time series
matlab -nodisplay -r "Combine_FacewiseTS('$subj','$sesh')"

#############################
#### module IV: Create figures
#############################
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"
echo "Starting module IV: Creating figures"
echo "ΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔΔ"

### filepaths for angular carpetplots
ATSl=[childFP '/' subj '_' sesh '_Prop_TS_dmn_L.csv'];
ATSr=[childFP '/' subj '_' sesh '_Prop_TS_dmn_L.csv'];

### make figure directory outside of scratch
mkdir /oak/stanford/groups/leanew1/users/apines/data/p50/${subj}/${sesh}/figs

### create interpolated carpetplot of bold
python3 Viz_ITS.py $subj $sesh

### create angular carpetplot 
python3 Viz_ATS.py $subj $sesh

### create positivity plots
# some type of radar plot?

### copy xcpd motion plots into figs dir

echo "OpFl complete"
