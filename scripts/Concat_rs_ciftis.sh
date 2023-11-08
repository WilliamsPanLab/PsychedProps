#!/bin/bash

input_directory=$1

# Collect the file paths into an array using find

files=()
while IFS= read -r -d '' file; do
  files+=("-cifti" "$file")
done < <(find "$input_directory" -type f -name "*task-rs*Smoothed_bold.dtseries.nii" -print0)

# Concatenate files
wb_command -cifti-merge "${input_directory}/concatenated_rs.dtseries.nii" "${files[@]}"

echo "Concatenation completed. The output file is 'concatenated_rs.dtseries.nii'."

