#!/bin/bash

ml biology
ml workbench

in_dir=/scratch/users/apines/VISspins_fs4_tmp
out_dir=/scratch/users/apines/VISspins_fs4_gii
mkdir -p $out_dir

# Sphere and area files (same as before)
sphere5_L=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii
sphere4_L=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.L.3k_fsavg_L.surf.gii
area5_L=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii
area4_L=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.L.midthickness_va_avg.3k_fsavg_L.shape.gii

sphere5_R=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii
sphere4_R=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4_std_sphere.R.3k_fsavg_R.surf.gii
area5_R=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii
area4_R=/oak/stanford/groups/leanew1/users/apines/standard_mesh_atlases/resample_fsaverage/fsaverage4.R.midthickness_va_avg.3k_fsavg_R.shape.gii

for i in $(seq -w 00001 10000); do
    echo $i
    wb_command -metric-resample ${in_dir}/perm${i}.L.func.gii $sphere5_L $sphere4_L ADAP_BARY_AREA ${out_dir}/perm${i}.L.3k.func.gii -area-metrics $area5_L $area4_L
    wb_command -metric-resample ${in_dir}/perm${i}.R.func.gii $sphere5_R $sphere4_R ADAP_BARY_AREA ${out_dir}/perm${i}.R.3k.func.gii -area-metrics $area5_R $area4_R
done
