#!/bin/bash

ml biology
ml workbench

in_dir=/scratch/users/apines/VISspins_fs4_tmp
out_dir=/scratch/users/apines/VISspins_fslr_gii
mkdir -p $out_dir

# Sphere and area files (same as before)
sphere5_L=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii
sphereLR_L=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii
area5_L=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii
areaLR_L=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii

sphere5_R=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii
sphereLR_R=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii
area5_R=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii
areaLR_R=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii

# Loop through all 10k spins
for i in $(seq -w 00001 10000); do
    echo "Upsampling permutation $i to fs_LR 32k"

    # Left hemisphere
    wb_command -metric-resample \
        ${in_dir}/perm${i}.L.func.gii \
        $sphere5_L \
        $sphereLR_L \
        ADAP_BARY_AREA \
        ${out_dir}/perm${i}.L.32k.func.gii \
        -area-metrics $area5_L $areaLR_L

    # Right hemisphere
    wb_command -metric-resample \
        ${in_dir}/perm${i}.R.func.gii \
        $sphere5_R \
        $sphereLR_R \
        ADAP_BARY_AREA \
        ${out_dir}/perm${i}.R.32k.func.gii \
        -area-metrics $area5_R $areaLR_R
done
