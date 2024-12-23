module load workbench
# set freesurfer dir
export FREESURFER_HOME=/share/software/user/open/freesurfer/6.0.0
# set subjs dir
export SUBJECTS_DIR=/share/software/user/open/freesurfer/6.0.0/subjects
# set freesurfer license
export FS_LICENSE=/oak/stanford/groups/leanew1/users/apines/license.txt

# Loop through all subjects from sub-MDMA001 to sub-MDMA017
for i in $(seq -f "%03g" 1 17); do
    subj=sub-MDMA${i}

    # Loop through all sessions from ses-00 to ses-03
    for j in $(seq -f "%02g" 0 3); do
        sesh=ses-${j}

        # mdma fp
        parentfp=/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/${subj}/${sesh}/func
        childfp=/scratch/users/xuezhang/data/mdma/${subj}/${sesh}

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
        # separate hemispheres - left - ap emotion
        wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_task-emotion_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii COLUMN -metric CORTEX_LEFT ${childfp}/${subj}_${sesh}_L_emotion_TS.func.gii
        # separate hemispheres - right - ap emotion
        wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_task-emotion_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii COLUMN -metric CORTEX_RIGHT ${childfp}/${subj}_${sesh}_R_emotion_TS.func.gii
        # separate hemispheres - left - ap gambling
        wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_task-gambling_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii COLUMN -metric CORTEX_LEFT ${childfp}/${subj}_${sesh}_L_gambling_TS.func.gii
        # separate hemispheres - right - ap gambling
        wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_task-gambling_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii COLUMN -metric CORTEX_RIGHT ${childfp}/${subj}_${sesh}_R_gambling_TS.func.gii
        # separate hemispheres - left - ap working memory
        wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_task-wm_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii COLUMN -metric CORTEX_LEFT ${childfp}/${subj}_${sesh}_L_wm_TS.func.gii
        # separate hemispheres - right - ap working memory
        wb_command -cifti-separate ${parentfp}/${subj}_${sesh}_task-wm_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii COLUMN -metric CORTEX_RIGHT ${childfp}/${subj}_${sesh}_R_wm_TS.func.gii

        ### resample both hemis to 10k vertices
        # rs1 - left
        wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_rs1_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_rs1_TS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii
        # rs1 - right
        wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_rs1_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_rs1_TS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii
        # rs2 - left
        wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_rs2_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_rs2_TS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii
        # rs2 - right
        wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_rs2_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_rs2_TS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii
        # emotion - left
        wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_emotion_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_emotion_TS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii
        # emotion - right
        wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_emotion_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_emotion_TS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii
        # gambling - left
        wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_gambling_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_gambling_TS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii
        # gambling - right
        wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_gambling_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_gambling_TS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii
        # wm - left
        wb_command -metric-resample ${childfp}/${subj}_${sesh}_L_wm_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_L_wm_TS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii
        # wm - right
        wb_command -metric-resample ${childfp}/${subj}_${sesh}_R_wm_TS.func.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii ADAP_BARY_AREA ${childfp}/${subj}_${sesh}_R_wm_TS_10k.func.gii -area-metrics /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii /oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii

        # convert to mgh for reading individual hemisphere time series into matlab
        # rs 1 - left
        /share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_L_rs1_TS_10k.func.gii ${childfp}/${subj}_${sesh}_L_rs1_TS_10k.mgh
        # rs 1 - right
        /share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_R_rs1_TS_10k.func.gii ${childfp}/${subj}_${sesh}_R_rs1_TS_10k.mgh
        # rs 2 - left
        /share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_L_rs2_TS_10k.func.gii ${childfp}/${subj}_${sesh}_L_rs2_TS_10k.mgh
        # rs 2 - right
        /share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_R_rs2_TS_10k.func.gii ${childfp}/${subj}_${sesh}_R_rs2_TS_10k.mgh
        # emotion - left
        /share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_L_emotion_TS_10k.func.gii ${childfp}/${subj}_${sesh}_L_emotion_TS_10k.mgh
        # emotion - right
        /share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_R_emotion_TS_10k.func.gii ${childfp}/${subj}_${sesh}_R_emotion_TS_10k.mgh
        # gambling - left
        /share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_L_gambling_TS_10k.func.gii ${childfp}/${subj}_${sesh}_L_gambling_TS_10k.mgh
        # gambling - right
        /share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_R_gambling_TS_10k.func.gii ${childfp}/${subj}_${sesh}_R_gambling_TS_10k.mgh
        # wm - left
        /share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_L_wm_TS_10k.func.gii ${childfp}/${subj}_${sesh}_L_wm_TS_10k.mgh
        # wm - right
        /share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${subj}_${sesh}_R_wm_TS_10k.func.gii ${childfp}/${subj}_${sesh}_R_wm_TS_10k.mgh

        echo "Finished processing $subj $sesh."
    done

    echo "Finished all sessions for $subj."
done


