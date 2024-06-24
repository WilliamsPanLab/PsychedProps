#!/bin/bash

# Set the destination folder
destination_folder="/scratch/users/apines/p50_mice/forTransfer"

# Create an associative array to keep track of pre file counts for each mouse
declare -A pre_counts

# Iterate over the source directories and find the .h5 files
for filepath in /scratch/users/apines/p50_mice/proc/20200228/*postLSD*/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5; do
    # Extract the mouse ID from the filepath
    mouse_id=$(echo $filepath | grep -oP 'm[0-9]+')
    echo $mouse_id
    # Increment the pre count for this mouse
    if [[ -z "${pre_counts[$mouse_id]}" ]]; then
        pre_counts[$mouse_id]=1
    else
        pre_counts[$mouse_id]=$((pre_counts[$mouse_id] + 1))
    fi

    # Determine the new filename
    new_filename="${mouse_id}_post${pre_counts[$mouse_id]}.h5"

    # Move the file to the destination folder with the new name
    cp "$filepath" "$destination_folder/$new_filename"
done

