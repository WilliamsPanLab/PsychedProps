subj=$1
sesh=$2

# set freesurfer dir
export FREESURFER_HOME=/share/software/user/open/freesurfer/6.0.0
# set subjs dir
export SUBJECTS_DIR=/share/software/user/open/freesurfer/6.0.0/subjects
# set freesurfer license
export FS_LICENSE=/oak/stanford/groups/leanew1/users/apines/license.txt

# parent dir
parentfp=/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/${subj}/${sesh}/func
# child dir is parent dir, all files in this script are intermediate files
childfp=${parentfp}

# subject's aggregated DMThalFC maps (dscalar)
DMThalFC=${parentfp}/${subj}_${sesh}_dmtFC.dscalar.nii

# separate hemispheres - left
wb_command -cifti-separate $DMThalFC COLUMN -metric CORTEX_LEFT ${childfp}/${subj}_${sesh}_L_dmtFC.func.gii 
# right hemi
wb_command -cifti-separate $DMThalFC COLUMN -metric CORTEX_RIGHT ${childfp}/${subj}_${sesh}_R_dmtFC.func.gii

### resample both hemis to 3k vertices
# left hemisphere
wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_dmtFC.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_dmtFC_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii

# right
wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_dmtFC.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_dmtFC_3k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii

