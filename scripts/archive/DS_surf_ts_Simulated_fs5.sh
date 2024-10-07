seed=$1
# set freesurfer dir
export FREESURFER_HOME=/share/software/user/open/freesurfer/6.0.0
# set subjs dir
export SUBJECTS_DIR=/share/software/user/open/freesurfer/6.0.0/subjects
# set freesurfer license
export FS_LICENSE=/oak/stanford/groups/leanew1/users/apines/license.txt

# subject's aggregated time series and other filenames
AgTS=/scratch/users/apines/ciftiout_Sym_${seed}.dtseries.nii
LeftHemi=/scratch/users/apines/Sim_L_AggTS_${seed}.func.gii
RightHemi=/scratch/users/apines/Sim_R_AggTS_${seed}.func.gii
LeftHemi_10=/scratch/users/apines/Sim_L_AggTS_${seed}_10k.func.gii
RightHemi_10=/scratch/users/apines/Sim_R_AggTS_${seed}_10k.func.gii
LeftHemi_10_mgh=/scratch/users/apines/Sim_L_AggTS_${seed}_10k.mgh
RightHemi_10_mgh=/scratch/users/apines/Sim_R_AggTS_${seed}_10k.mgh


# separate hemispheres - left
wb_command -cifti-separate $AgTS COLUMN -metric CORTEX_LEFT ${LeftHemi}

# right hemi
wb_command -cifti-separate $AgTS COLUMN -metric CORTEX_RIGHT ${RightHemi}

### resample both hemis to 10k vertices
# left hemisphere
wb_command -metric-resample ${LeftHemi} /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ADAP_BARY_AREA ${LeftHemi_10} -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii

# right hemisphere
wb_command -metric-resample ${RightHemi} /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii ADAP_BARY_AREA ${RightHemi_10} -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii

# convert to mgh for reading individual hemisphere time series into matlab
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${LeftHemi_10} ${LeftHemi_10_mgh}
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${RightHemi_10} ${RightHemi_10_mgh}
