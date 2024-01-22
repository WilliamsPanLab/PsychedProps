#!/bin/bash

input_directory=$1
output_directory="/scratch/groups/leanew1/xcpd_outP50_36p_bp/bv_rs"

# Loop over subdirectories in the input directory
for subdirectory in "$input_directory"/*; do
  # Check if the subdirectory is a directory
  if [ -d "$subdirectory" ]; then
    # Extract subject name from the subdirectory
    subj=$(basename "$subdirectory")

    # Set the correct subdirectory path based on the subject name
    if [[ "$subj" == sub-MDMA* ]]; then
      subj_directory="$subdirectory/ses-00/func"
    elif [[ "$subj" == sub-CONN* ]]; then
      subj_directory="$subdirectory/ses-BL/func"
    else
      echo "Unknown subject format: $subj. Skipping."
      continue
    fi

    echo "Searching in directory: $subj_directory"

    # Collect the file paths into an array using find command
    files=()
    while IFS= read -r -d '' file; do
      files+=("-cifti" "$file")
    done < <(find "$subj_directory" -type f -name "*denoisedSmoothed_bold.dtseries.nii" -print0)

    # Check if any files were found
    if [ ${#files[@]} -eq 0 ]; then
      echo "No files found for subject $subj. Skipping."
      continue
    fi

    # Concatenate the files using wb_command
    output_file="$output_directory/${subj}_concatenated.dtseries.nii"
    wb_command -cifti-merge "$output_file" "${files[@]}"

    echo "Concatenation completed for subject $subj. The output file is '$output_file'."
  fi
done

