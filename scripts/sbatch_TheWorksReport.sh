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

# will need matlab
module load matlab
# and freesurfer
module load biology
module load freesurfer/7.3.2
# and workbench
module load workbench

# subject name is input argument
subj=$1
# sesh is input 2
sesh=$2

# Downsample the data 
# /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/fs_5/DS_surf_ts_mdma_fs5.sh $1 $2

# cd to scripts directory
# cd /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/fs_5

# Calculate Optical Flow
# matlab -nodisplay -r "OpFl_mdma_fs5('$subj','$sesh')"

# RS filepaths
# childfp=/scratch/users/apines/data/mdma/${subj}/${sesh}
# rsIn=${childfp}/${subj}_${sesh}_OpFl_rs_fs5.mat
# rsOut=${childfp}/${subj}_${sesh}_PGGDist_rs_fs5.mat

# make output directory outside scratch
# mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj} 
mkdir -p /oak/stanford/groups/leanew1/users/apines/data/p50/${subj}/${sesh}/
# ./run_Extract_BUTD_ResultantVecs_Gran_fs5.sh /share/software/user/restricted/matlab/R2018a/ $rsIn $childfp/OpFl_timeseries_L_fs5.mat $childfp/OpFl_timeseries_R_fs5.mat

# make a simple opfl txt file with whole-cortex amplitude and SD time series
# matlab -nodisplay -r "OpFl_AmpSD('$subj','$sesh')"

# this is dumb, but switch back
#cd /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts

# derive FC matrix
#matlab -nodisplay -r "fs4FC('$subj','$sesh')"

# derive PG
#module load python/3.9
#source /oak/stanford/groups/leanew1/users/apines/dmap/bin/activate
#python fcmat_to_threshCosSim.py $subj $sesh

# visualize PG
#PGpng=${childfp}_${subj}_${sesh}_PG.png
#matlab -nodisplay -r "vis_fs4PG('$subj','$sesh','$PGpng')"

# extract time course of precuneus
# extract ROIs into text file (50 and 100 are precun)
#wb_command -cifti-convert -to-text ${childfp}/${subj}_${sesh}_rs_concat_Parcellated.ptseries.nii $ROITS

#### set xcpd_outdir
# MDMA
xcpd_outdir=/scratch/groups/leanew1/xcpd_outMDMA_36p_despike_bp/xcp_d/${subj}/${sesh}/func/
# PTRD
#xcpd_outdir=/scratch/groups/leanew1/xcpd_outPTRD_36p_despike_bp/xcp_d/${subj}/${sesh}/func/
# P50 (KET)
# xcpd_outdir=/scratch/groups/leanew1/xcpd_outP50_36p_despike_bp/xcp_d/${subj}/${sesh}/func/
# set script outdir
# childfp=/oak/stanford/groups/leanew1/users/apines/data/p50/${subj}/${sesh}/

# extract subcortical time courses
rs1xcpSubcort_fp=${xcpd_outdir}${subj}_${sesh}_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_atlas-subcortical_den-91k_timeseries.ptseries.nii
rs2xcpSubcort_fp=${xcpd_outdir}${subj}_${sesh}_task-rs_acq-mb_dir-pe1_run-0_space-fsLR_atlas-subcortical_den-91k_timeseries.ptseries.nii
#wmxcpSubcort_fp=${xcpd_outdir}${subj}_${sesh}_task-wm_acq-mb_dir-pe0_run-0_space-fsLR_atlas-subcortical_den-91k_timeseries.ptseries.nii
#midxcpSubcort_fp=${xcpd_outdir}${subj}_${sesh}_task-mid_acq-mb_dir-pe0_run-0_space-fsLR_atlas-subcortical_den-91k_timeseries.ptseries.nii
#gambxcpSubcort_fp=${xcpd_outdir}${subj}_${sesh}_task-gambling_acq-mb_dir-pe0_run-0_space-fsLR_atlas-subcortical_den-91k_timeseries.ptseries.nii
#emoxcpSubcort_fp=${xcpd_outdir}${subj}_${sesh}_task-emotion_acq-mb_dir-pe0_run-0_space-fsLR_atlas-subcortical_den-91k_timeseries.ptseries.nii
SubcortTS1=${childfp}${subj}_${sesh}_rs_SubCortROIS.txt
SubcortTS2=${childfp}${subj}_${sesh}_rs2_SubCortROIS.txt
#SubcortTS3=${childfp}${subj}_${sesh}_wm_SubCortROIS.txt
# note mdma doesnt have multiband mid
#SubcortTS4=${childfp}${subj}_${sesh}_mid_SubCortROIS.txt
#SubcortTS5=${childfp}${subj}_${sesh}_gambling_SubCortROIS.txt
#SubcortTS6=${childfp}${subj}_${sesh}_emotion_SubCortROIS.txt

# subcortical
######### extract dorsomedial thalamus time series: from 9 and 34 from Tian_Subcortex_S3_3T
wb_command -cifti-convert -to-text $rs1xcpSubcort_fp $SubcortTS1
wb_command -cifti-convert -to-text $rs2xcpSubcort_fp $SubcortTS2
# scale 3 tian parc location
scale3fp=/oak/stanford/groups/leanew1/users/apines/maps/Tian_Subcortex_S3_3T_32k.dlabel.nii
# should match labels here https://github.com/yetianmed/subcortex/blob/master/Group-Parcellation/3T/Subcortex-Only/Tian_Subcortex_S3_3T_label.txt

# extract peaks, delays, and magnitudes from ptseries and global signal
#ml python/3
#python3 /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/derive_pulses.py ${subj} ${sesh}

# extract GS from tsv
rs1ConfTsv=${xcpd_outdir}${subj}_${sesh}_task-rs_acq-mb_dir-pe0_run-0_design.tsv
rs2ConfTsv=${xcpd_outdir}${subj}_${sesh}_task-rs_acq-mb_dir-pe1_run-0_design.tsv

####### extract FC of DM thalamus
matlab -nodisplay -r "ROIfc('$subj','$sesh')"

########3# make grayplot of grayordinates org. by FC to those regions (sep L and R)

#python Viz_grayplots.py $subj $sesh $childfp



# extract amygdalar TS

# plot both together

# Optical flow streamlines
#cd /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/fs_5
#matlab -nodisplay -r "OpFlStreamlines('$subj',$sesh','Streamlines')"

# optical flow CFC

# make grayplots of 3k mgh file (thalamus, whole-brain, OpFl Ampltiude) along with time series (GS, thalamus, precuneus, OpFl Amplitude)
# bring it all together: generate grayplots and power spectral density within
#python Viz_grayplots.py $subj $sesh $childfp


#### THIS WILL COME AFTER SESH-LEVEL RUNS ARE COMPLETED INDIVIDUALLY
# t-test FC mat
# make cfc
# t-test CFC
# Optical flow t-test (amplitude)
# Optical flow t-test (TD%)
# pull it together into one visual report

