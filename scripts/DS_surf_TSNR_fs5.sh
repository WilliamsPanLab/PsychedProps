# set freesurfer dir
export FREESURFER_HOME=/share/software/user/open/freesurfer/6.0.0
# set subjs dir
export SUBJECTS_DIR=/share/software/user/open/freesurfer/6.0.0/subjects
# set freesurfer license
export FS_LICENSE=/oak/stanford/groups/leanew1/users/apines/license.txt

# parent dir
parentfp=/oak/stanford/groups/leanew1/users/apines/maps

# child dir is parent dir, all files in this script are intermediate files
childfp=${parentfp}

# aggregated tsnr maps (dscalar)
AgNet=${parentfp}/average_tsnr_rs_mdma.dscalar.nii

# separate hemispheres - left
wb_command -cifti-separate $AgNet COLUMN -metric CORTEX_LEFT ${childfp}/mdma_L_AggTSNR.func.gii 
# right hemi
wb_command -cifti-separate $AgNet COLUMN -metric CORTEX_RIGHT ${childfp}/mdma_R_AggTSNR.func.gii

### resample both hemis to 10k vertices
# left hemisphere
wb_command -metric-resample ${childfp}/mdma_L_AggTSNR.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/mdma_L_AggTSNR_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii
	
# right hemisphere
wb_command -metric-resample ${childfp}/mdma_R_AggTSNR.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/mdma_R_AggTSNR_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii