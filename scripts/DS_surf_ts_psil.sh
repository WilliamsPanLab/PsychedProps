subj=$1
sesh=$2

# set freesurfer dir
export FREESURFER_HOME=/share/software/user/open/freesurfer/6.0.0
# set subjs dir
export SUBJECTS_DIR=/share/software/user/open/freesurfer/6.0.0/subjects
# set freesurfer license
export FS_LICENSE=/oak/stanford/groups/leanew1/users/apines/license.txt
# psil fp
parentfp=/oak/stanford/groups/leanew1/SHARED_DATASETS/private/WashU_psilocybin/${subj}/${subj}_${sesh}/func
childfp=/scratch/users/apines/data/psil/${subj}/${sesh}

# make output dirs
mkdir ${childfp} -p

# separate hemispheres - resting state 1
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_rs1.dtseries.nii COLUMN -metric CORTEX_LEFT ${childfp}/${subj}_${sesh}_L_rs1_TS.func.gii
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_rs1.dtseries.nii COLUMN -metric CORTEX_RIGHT ${childfp}/${subj}_${sesh}_R_rs1_TS.func.gii

# separate hemispheres - resting state 2
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_rs2.dtseries.nii COLUMN -metric CORTEX_LEFT ${childfp}/${subj}_${sesh}_L_rs2_TS.func.gii
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_rs2.dtseries.nii COLUMN -metric CORTEX_RIGHT ${childfp}/${subj}_${sesh}_R_rs2_TS.func.gii

# resting state 3
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_rs3.dtseries.nii COLUMN -metric CORTEX_LEFT ${childfp}/${subj}_${sesh}_L_rs3_TS.func.gii
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_rs3.dtseries.nii COLUMN -metric CORTEX_RIGHT ${childfp}/${subj}_${sesh}_R_rs3_TS.func.gii

# resting state 4
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_rs4.dtseries.nii COLUMN -metric CORTEX_LEFT ${childfp}/${subj}_${sesh}_L_rs4_TS.func.gii
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_rs4.dtseries.nii COLUMN -metric CORTEX_RIGHT ${childfp}/${subj}_${sesh}_R_rs4_TS.func.gii

# resting state 5
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_rs5.dtseries.nii COLUMN -metric CORTEX_LEFT ${childfp}/${subj}_${sesh}_L_rs5_TS.func.gii
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_rs5.dtseries.nii COLUMN -metric CORTEX_RIGHT ${childfp}/${subj}_${sesh}_R_rs5_TS.func.gii

# resting state 6
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_rs6.dtseries.nii COLUMN -metric CORTEX_LEFT ${childfp}/${subj}_${sesh}_L_rs6_TS.func.gii
wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_rs6.dtseries.nii COLUMN -metric CORTEX_RIGHT ${childfp}/${subj}_${sesh}_R_rs6_TS.func.gii

### resample both hemis to 2.5k vertices
# rs1 - left
wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_rs1_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_rs1_TS_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii
# rs1 - right
wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_rs1_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_rs1_TS_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii
# rs2 - left
wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_rs2_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_rs2_TS_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii
# rs2 - right
wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_rs2_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_rs2_TS_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii
# rs3 - left
wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_rs3_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_rs3_TS_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii
# rs3 - right
wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_rs3_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_rs3_TS_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii
# rs4 - left
wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_rs4_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_rs4_TS_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii
# rs4 - right
wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_rs4_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_rs4_TS_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii
# rs5 - left
wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_rs5_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_rs5_TS_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii
# rs5 - right
wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_rs5_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_rs5_TS_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii
# rs6 - left
wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_rs6_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_rs6_TS_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii
# rs6 - right
wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_rs6_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_rs6_TS_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii

# convert to mgh for reading individual hemisphere time series into matlab
# rs 1 - left
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_L_rs1_TS_3k.func.gii ${childfp}/${subj}_${sesh}_L_rs1_TS_3k.mgh
# rs 1 - right
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_R_rs1_TS_3k.func.gii ${childfp}/${subj}_${sesh}_R_rs1_TS_3k.mgh
# rs 2 - left
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_L_rs2_TS_3k.func.gii ${childfp}/${subj}_${sesh}_L_rs2_TS_3k.mgh
# rs 2 - right
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_R_rs2_TS_3k.func.gii ${childfp}/${subj}_${sesh}_R_rs2_TS_3k.mgh
# rs 3 - left
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_L_rs3_TS_3k.func.gii ${childfp}/${subj}_${sesh}_L_rs3_TS_3k.mgh
# rs 3 - right
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_R_rs3_TS_3k.func.gii ${childfp}/${subj}_${sesh}_R_rs3_TS_3k.mgh
# rs 4 - left
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_L_rs4_TS_3k.func.gii ${childfp}/${subj}_${sesh}_L_rs4_TS_3k.mgh
# rs 4 - right
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_R_rs4_TS_3k.func.gii ${childfp}/${subj}_${sesh}_R_rs4_TS_3k.mgh
# rs 5 - left
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_L_rs5_TS_3k.func.gii ${childfp}/${subj}_${sesh}_L_rs5_TS_3k.mgh
# rs 5 - right
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_R_rs5_TS_3k.func.gii ${childfp}/${subj}_${sesh}_R_rs5_TS_3k.mgh
# rs 6 - left
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_L_rs6_TS_3k.func.gii ${childfp}/${subj}_${sesh}_L_rs6_TS_3k.mgh
# rs 6 - right
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_R_rs6_TS_3k.func.gii ${childfp}/${subj}_${sesh}_R_rs6_TS_3k.mgh
