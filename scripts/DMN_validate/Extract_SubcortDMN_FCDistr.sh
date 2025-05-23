#!/bin/bash
ml biology
ml workbench

# input files
dscalar_file=~/GroupAvg_DMNFC_map.dscalar.nii
dlabel_file=/oak/stanford/groups/leanew1/users/apines/maps/Tian_Subcortex_S1_3T_32k.dscalar.nii

# temporary outputs
temp_dscalar_txt=~/tmp_dscalar_values.txt
temp_label_indices_txt=~/tmp_label_indices.txt
# permanent
output_csv=~/DMNFC_voxelwise_values.csv

# extract label index per voxel
wb_command -cifti-convert -to-text "$dscalar_file" "$temp_dscalar_txt"

# Extract label indices per voxel (as integers from dlabel)
wb_command -cifti-convert -to-text "$dlabel_file" "$temp_label_indices_txt"

# Combine into CSV: label,value
paste -d, "$temp_label_indices_txt" "$temp_dscalar_txt" > "$output_csv"
echo "Saved all subcortical label indices and DMNFC values to:"
echo "$output_csv"

