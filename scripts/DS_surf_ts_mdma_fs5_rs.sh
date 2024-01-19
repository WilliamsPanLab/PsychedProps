subj=$1
sesh=$2

# set freesurfer dir
export FREESURFER_HOME=/share/software/user/open/freesurfer/6.0.0
# set subjs dir
export SUBJECTS_DIR=/share/software/user/open/freesurfer/6.0.0/subjects
# set freesurfer license
export FS_LICENSE=/oak/stanford/groups/leanew1/users/apines/license.txt
# mdma fp
parentfp=/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/${subj}/${sesh}/func
childfp=/scratch/users/apines/data/mdma/${subj}/${sesh}

# make output dirs
mkdir ${childfp} -p

# separate hemispheres - left - ap resting state
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii COLUMN -metric CORTEX_LEFT ${childfp}/${subj}_${sesh}_L_rs1_TS.func.gii
# separate hemispheres - right - ap resting state
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii COLUMN -metric CORTEX_RIGHT ${childfp}/${subj}_${sesh}_R_rs1_TS.func.gii
# separate hemispheres - left - pa resting state
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_task-rs_acq-mb_dir-pe1_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii COLUMN -metric CORTEX_LEFT ${childfp}/${subj}_${sesh}_L_rs2_TS.func.gii
# separate hemispheres - right - pa resting state
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_task-rs_acq-mb_dir-pe1_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii COLUMN -metric CORTEX_RIGHT ${childfp}/${subj}_${sesh}_R_rs2_TS.func.gii

### resample both hemis to 10k vertices
# rs1 - left
wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_rs1_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_rs1_TS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii
# rs1 - right
wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_rs1_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_rs1_TS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii
# rs2 - left
wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_rs2_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_rs2_TS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii
# rs2 - right
wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_rs2_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_rs2_TS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii

# convert to mgh for reading individual hemisphere time series into matlab
# rs 1 - left
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_L_rs1_TS_10k.func.gii ${childfp}/${subj}_${sesh}_L_rs1_TS_10k.mgh
# rs 1 - right
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_R_rs1_TS_10k.func.gii ${childfp}/${subj}_${sesh}_R_rs1_TS_10k.mgh
# rs 2 - left
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_L_rs2_TS_10k.func.gii ${childfp}/${subj}_${sesh}_L_rs2_TS_10k.mgh
# rs 2 - right
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_R_rs2_TS_10k.func.gii ${childfp}/${subj}_${sesh}_R_rs2_TS_10k.mgh

