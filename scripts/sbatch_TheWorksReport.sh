#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=7:00:00
#SBATCH -n 1
#SBATCH --mem=85G
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
/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/fs_5/DS_surf_ts_mdma_fs5.sh $1 $2

# cd to scripts directory
cd /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/fs_5

# Calculate Optical Flow
matlab -nodisplay -r "OpFl_mdma_fs5('$subj','$sesh')"

# RS filepaths
childfp=/scratch/users/apines/data/mdma/${subj}/${sesh}
rsIn=${childfp}/${subj}_${sesh}_OpFl_rs_fs5.mat
rsOut=${childfp}/${subj}_${sesh}_PGGDist_rs_fs5.mat

# make output directory outside scratch
mkdir /oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/${subj} 

./run_Extract_BUTD_ResultantVecs_Gran_fs5.sh /share/software/user/restricted/matlab/R2018a/ $rsIn $childfp/OpFl_timeseries_L_fs5.mat $childfp/OpFl_timeseries_R_fs5.mat

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
#AgTS=${childfp}/${subj}_${sesh}_rs_concat.dtseries.nii
#ROITS=${childfp}_${subj}_${sesh}_CortROIS.txt
# parcellate into ROIs
#wb_command -cifti-parcellate $AgTS /oak/stanford/groups/leanew1/users/apines/maps/Schaefer2018_100Parcels_7Networks_order.dlabel.nii COLUMN ${childfp}/${subj}_${sesh}_rs_concat_Parcellated.ptseries.nii
# extract ROIs into text file (50 and 100 are precun)
#wb_command -cifti-convert -to-text ${childfp}/${subj}_${sesh}_rs_concat_Parcellated.ptseries.nii $ROITS

# extract time course of thalamus
#rs1xcpSubcort_fp=/scratch/groups/leanew1/xcpd_outMDMA/xcp_d/${subj}/${sesh}/func/${subj}_${sesh}_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_atlas-subcortical_den-91k_timeseries.ptseries.nii
#rs2xcpSubcort_fp=/scratch/groups/leanew1/xcpd_outMDMA/xcp_d/${subj}/${sesh}/func/${subj}_${sesh}_task-rs_acq-mb_dir-pe1_run-0_space-fsLR_atlas-subcortical_den-91k_timeseries.ptseries.nii
#SubcortTS1=${childfp}_${subj}_${sesh}_SubCortROIS1.txt
#SubcortTS2=${childfp}_${subj}_${sesh}_SubCortROIS2.txt
#wb_command -cifti-convert -to-text $rs1xcpSubcort_fp $SubcortTS1
#wb_command -cifti-convert -to-text $rs2xcpSubcort_fp $SubcortTS2
# should match labels here https://github.com/yetianmed/subcortex/blob/master/Group-Parcellation/3T/Subcortex-Only/Tian_Subcortex_S3_3T_label.txt


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

