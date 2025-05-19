#!/bin/bash

data_dir="/scratch/groups/leanew1/xcpd_outP50_36p_bp/bv_rs"
subjects=($(ls "$data_dir"/sub-CONN*_concatenated.dtseries.nii | sed -E 's|.*/(sub-CONN[0-9]+)_.*|\1|' | sort -u))
# Loop through subjects
for subject in "${subjects[@]}"; do
	# Submit sbatch job
       	sbatch sbatch_DMNFC_map.sh $subject
	echo $subject
done
