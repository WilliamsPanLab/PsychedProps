
# Script Review and Comments

## Summary of Required Dependencies and Recommendations

### General Requirements for the Repository
This repository contains a mix of **Bash**, **Python**, and **MATLAB** scripts. To run these scripts, you need access to the following software and data directories:
- **Software**:
  - Freesurfer
  - Workbench
  - Python 3 (with `nibabel` installed)
  - `xcp_d` and `fmriprep` singularity files
  - Matlab R2022b (recommended; R2020a may cause errors)
- **Data Directories** (consider including these into the current repo):
  - `/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder`
  - `/oak/stanford/groups/leanew1/users/apines/fs5surf/lh.Medial_wall.label`
  - `/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases`
  - `/home/users/apines/sbatchTemplate.sh`
  - `/oak/stanford/groups/leanew1/users/apines/fs5surf`
  - `/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/`
  - `/oak/stanford/groups/leanew1/users/apines/libs`
  - `/oak/stanford/groups/leanew1/users/apines/fs4surf`
  - `/oak/stanford/groups/leanew1/users/apines/libs/lukaslang-ofd-614a2ffc50d6`

### General Recommendations
1. **Script Naming**: Rename scripts to a consistent format like `STEPxx_xx_Description` to ensure clarity in workflow order.
2. **Code Cleanup**: Remove any unused scripts to reduce confusion.
3. **Script Order**: 
   - Run `xcp_d` immediately after `fmriprep` to streamline time-consuming processes.
   - Consider integrating all preprocessing steps into `sbatch_OpFI.sh`?
8. **Downsampling Clarification**: Explain why `NMF` is run on `fs5` (10k vertices) while `Optical Flow` uses `fs4` (3k vertices), as this could cause confusion.
9. **Log Directory**: Ensure `/log` exists before running scripts like `fmriprep` that generate logs.
6. **Usage Guide**: 
   - Include instructions on running the scripts, specifying required modules and software.
   - Add a reference link to the Sherlock HPC documentation for job submission.
### Example Usage
For ease of understanding, include example commands to run each script. For instance:
```bash
ml biology
ml workbench
matlab -nodesktop < TSNR_mask_3_avSNR-mw_to_label.m
```

---

## Specific Script Comments

### 1. `full_fmriprep.sh`
- **Input Data**: Update to use `p50` dataset instead of `connectome`.

### 2. `TSNR_mask_0_antimask.sh`
- **Requirements**: Needs FSL and Workbench.
- **Usage Scope**: Is this script intended only for **MDMA**? Please confirm.
- **Changes Made**: modified the script to loop for subjects and sessions and added a helper script to submit the job: `sbatch_TSNR_mask_0.sh`

### 3. `TSNR_mask_1_ExtractSNR_subjectwise.py`
- **Requirements**: 
  - Load the Python module (Python 3) before running.
  - Ensure `nibabel` is installed.
- **Execution**: `python3 -u TSNR_mask_1_ExtractSNR_subjectwise.py`
- **File Access**: Requires access to `hcp.gradients.dscalar.nii`.
- **Changes Made**: modified the script to loop for subjects and sessions and added a helper script to submit the job: `sbatch_TSNR_mask_1.sh`

### 4. `TSNR_mask_2_Combine_SNRmaps.py`
- **File Permissions**: Modify permissions for `sub-MDMA0*meanTSNR.dscalar.nii`.
- **Output Folder**: Ensure the output directory is created.

### 5. `/DS_surf_TSNR_fs5.sh`
- **Execution**: Run on `reward/attention`. Requires `wb_command` and FreeSurfer.

### 6. `TSNR_mask_3_avSNR-mw_to_label.m`
- **Tool Requirements**: 
  - Tool folder required at `/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder`.
  - Verify if **xcp** can proceed before the `tSNR mask` is created.
- **Matlab Execution**: `matlab -nodesktop < TSNR_mask_3_avSNR-mw_to_label.m`
- **File Access**: Needs access to `/oak/stanford/groups/leanew1/users/apines/fs5surf/lh.Medial_wall.label`.

### 7. `DS_surf_ts_mdma_fs5.sh`
- **Changes Made**: modified the script to loop for subjects and sessions and added a helper script to submit the job: `sbatch_DS_surf_ts_mdma_fs5.sh`
- **Requirements**: Needs `Freesurfer` and `wb_command`; access to `/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases`.
- **Execution on Sherlock**: Load `biology` and `workbench` modules.

### 8. `MotMask.m`
- **Execution Order**: This script is likely to be run within `sbatch_OpFI.sh` after `DS_surf_ts_mdma.sh`, as it checks for downsampled data **3k**.

### 9. `Step_1_CreatePrepData.m`
- **Access Requirements**: Needs access to the tool folder and `/oak/stanford/groups/leanew1/users/apines/fs5surf/`.
- **Execution**: Load Matlab and run `matlab -nodesktop < Step_1_CreatePrepData.m`.
- **Output**: Specify the expected output in the README for clarity.

### 10. `Step_2nd_ParcellationInitialize.m`
- **README Clarification**: README suggests determining the optimal solution, but `k=4` appears pre-determined.
- **Access Requirements**: Access to `cp /home/users/apines/sbatchTemplate.sh` and `lh.Mask_SNR.label`.

### 11. `Step_3rd_SelRobustInit.m`
- **Access Requirements**: Needs access to the tool folder.
- **MATLAB Version**: Use R2022b.
- **Execution**: `matlab -nodisplay < Step_3rd_SelRobustInit.m`

### 12. `Step_4_Extract_Group_Atlas.m`
- **Modules**: Load `biology` and `workbench`.
- **MATLAB Compatibility**: R2022b is required due to compatibility issues with `wb_command`.

### 13. `DS_surf_Networks_fs5tofs4.sh`
- **Access Requirements**: Needs access to `/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii`.

### 14. `NMF/Netgiis_2_mat.m`
- **Duplication**: Consider removing the script outside the `NMF` folder to avoid confusion.

### 15. `OpFL_mdma.m`
- **Dependencies**: 
  - Access to `/oak/stanford/groups/leanew1/users/apines/libs/`.
  - Access to `/oak/stanford/groups/leanew1/users/apines/fs4surf/`.
