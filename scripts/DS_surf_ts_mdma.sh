subj=$1
sesh=$2

# set freesurfer dir
export FREESURFER_HOME=/share/software/user/open/freesurfer/6.0.0
# set subjs dir
export SUBJECTS_DIR=/share/software/user/open/freesurfer/6.0.0/subjects
# set freesurfer license
export FS_LICENSE=/oak/stanford/groups/leanew1/users/apines/license.txt
# mdma fp
parentfp=/scratch/groups/leanew1/xcpd_outMDMA_36p_despike_bp/xcp_d/${subj}/${sesh}/func
childfp=/scratch/users/apines/data/mdma/${subj}/${sesh}

# make output dirs
mkdir ${childfp} -p

# concat phase encoding directions
wb_command -cifti-merge ${childfp}/${subj}_${sesh}_rs_concat.dtseries.nii -cifti ${parentfp}/${subj}_${sesh}_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii -cifti ${parentfp}/${subj}_${sesh}_task-rs_acq-mb_dir-pe1_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii

# subject's aggregated time series
AgTS=${childfp}/${subj}_${sesh}_rs_concat.dtseries.nii

# separate hemispheres - left
wb_command -cifti-separate $AgTS COLUMN -metric CORTEX_LEFT ${childfp}/${subj}_${sesh}_L_AggTS.func.gii 

# right hemi
wb_command -cifti-separate $AgTS COLUMN -metric CORTEX_RIGHT ${childfp}/${subj}_${sesh}_R_AggTS.func.gii

### resample both hemis to 2.5k vertices
# left hemisphere
wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_AggTS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_AggTS_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii

# right hemisphere
wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_AggTS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_AggTS_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii

# convert to mgh for reading individual hemisphere time series into matlab
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_L_AggTS_3k.func.gii ${childfp}/${subj}_${sesh}_L_AggTS_3k.mgh
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_R_AggTS_3k.func.gii ${childfp}/${subj}_${sesh}_R_AggTS_3k.mgh


