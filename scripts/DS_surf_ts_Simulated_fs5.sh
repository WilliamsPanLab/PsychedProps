# set freesurfer dir
export FREESURFER_HOME=/share/software/user/open/freesurfer/6.0.0
# set subjs dir
export SUBJECTS_DIR=/share/software/user/open/freesurfer/6.0.0/subjects
# set freesurfer license
export FS_LICENSE=/oak/stanford/groups/leanew1/users/apines/license.txt

# subject's aggregated time series
AgTS=~/ciftiout_Sym_Smooth.dtseries.nii

# separate hemispheres - left
wb_command -cifti-separate $AgTS COLUMN -metric CORTEX_LEFT ~/Sim_L_AggTS.func.gii

# right hemi
wb_command -cifti-separate $AgTS COLUMN -metric CORTEX_RIGHT ~/Sim_R_AggTS.func.gii

### resample both hemis to 10k vertices
# left hemisphere
wb_command -metric-resample ~/Sim_L_AggTS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ADAP_BARY_AREA ~/Sim_L_AggTS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii

# right hemisphere
wb_command -metric-resample ~/Sim_R_AggTS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii ADAP_BARY_AREA ~/Sim_R_AggTS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii

# convert to mgh for reading individual hemisphere time series into matlab
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ~/Sim_L_AggTS_10k.func.gii ~/Sim_L_AggTS_10k.mgh
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ~/Sim_R_AggTS_10k.func.gii ~/Sim_R_AggTS_10k.mgh
