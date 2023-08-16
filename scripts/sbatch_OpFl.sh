#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=0:20:00
#SBATCH -n 1
#SBATCH --mem=35G
#SBATCH -p leanew1,normal  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=no@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------


########################################
##             _.-````'-,_
##
##   _,.,_ ,-'`           `'-.,_
##
## /)     (\                   '``-.
##
##((      ) )                      `\
##
## \)    (_/                        )\
##
##  |       /)           '    ,'    / \
##
##  `\    ^'            '     (    /  ))
##
##    |      _/\ ,     /    ,,`\   (  "`
##
##     \Y,   |  \  \  | ````| / \_ \##
# This script is to prep the environment of a slurm node, 
##       `)_/    \  \  )    ( >  ( >
# launch OpFl (Orientation and Positivity of Functional Images),
##                \( \(     |/   |/
# and delete interim files upon completion
##    mic & dwb  /_(/_(    /_(  /_(


############################
#### module I: Preprocessing
############################

echo "Starting module I: preproc"

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

# cd to scripts directory
cd /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/fs_5


# mask images: 10+ continuous frames only
matlab -nodisplay -r "MotMask('$subj','$sesh')"
### add feed-in for FD time series

############################
#### module II: Optical Flow
############################

echo "Starting module II: Optical Flow"

# Calculate Optical Flow
matlab -nodisplay -r "OpFl_mdma_fs5('$subj','$sesh')"

# RS filepaths
childfp=/scratch/users/apines/data/mdma/${subj}/${sesh}
rsIn=${childfp}/${subj}_${sesh}_OpFl_rs_fs5.mat
rsOut=${childfp}/${subj}_${sesh}_PGGDist_rs_fs5.mat

# interpolate fs5 time series to faces and between-timepoints

#############################
#### module III: Calc. Angles
#############################

echo "Starting module III: Angular distance calculation"

# make output directory outside scratch
mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj} 
mkdir -p /oak/stanford/groups/leanew1/users/apines/data/p50/${subj}/${sesh}/
./run_Extract_BUTD_ResultantVecs_Gran_fs5.sh /share/software/user/restricted/matlab/R2018a/ $rsIn $childfp/OpFl_timeseries_L_fs5.mat $childfp/OpFl_timeseries_R_fs5.mat

# make a simple opfl txt file with whole-cortex amplitude and SD time series
# matlab -nodisplay -r "OpFl_AmpSD('$subj','$sesh')"

# this is dumb, but switch back
#cd /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts


#############################
#### module IV: Plot report
#############################

echo "Starting module IV: Plotting report"

# extract time course of precuneus
# extract ROIs into text file (50 and 100 are precun)
wb_command -cifti-convert -to-text ${childfp}/${subj}_${sesh}_rs_concat_Parcellated.ptseries.nii $ROITS

# extract peaks, delays, and magnitudes from ptseries and global signal
#ml python/3
#python3 /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/derive_pulses.py ${subj} ${sesh}

# extract GS from tsv
rs1ConfTsv=${xcpd_outdir}${subj}_${sesh}_task-rs_acq-mb_dir-pe0_run-0_design.tsv
rs2ConfTsv=${xcpd_outdir}${subj}_${sesh}_task-rs_acq-mb_dir-pe1_run-0_design.tsv

####### extract FC of DM thalamus
matlab -nodisplay -r "ROIfc('$subj','$sesh')"

### transform this into grayplot viz: interpolated BOLD grayplot, GS, Orientation grayplot, half radar for positvity by delta dmn and negativity by delta dmn
python3 Viz_grayplots.py $subj $sesh $childfp

# extract amygdalar TS


# make grayplots of 3k mgh file (thalamus, whole-brain, OpFl Ampltiude) along with time series (GS, thalamus, precuneus, OpFl Amplitude)
# bring it all together: generate grayplots and power spectral density within
#python Viz_grayplots.py $subj $sesh $childfp

# quick cleanup of interim files?


echo "OpFl complete"
