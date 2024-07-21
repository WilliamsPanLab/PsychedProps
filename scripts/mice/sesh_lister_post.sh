#!/bin/bash

# Change to the directory
cd /scratch/users/apines/p50_mice/proc

# Use find to search for 'masked_dff.h5' files in directories matching the pattern '*pre*_1'
find . -type f -name "masked_dff.h5" -path "*/20200228/*postLSD*/masked_dff.h5" | while read filepath; do
# Extract the date directory and the mouse directory from the path
    date_dir=$(echo $filepath | cut -d'/' -f2)
    mouse_dir=$(echo $filepath | cut -d'/' -f3)

    # Check if the mouse directory matches the 'preXXX_1' pattern
    if [[ "$mouse_dir" =~ post* ]]; then
        # Print the date and mouse directory
        echo "${date_dir} ${mouse_dir}"
    fi
done > /scratch/users/apines/p50_mice/demo/ipynb/sess_postDrug.txt
