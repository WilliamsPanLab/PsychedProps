#!/bin/bash

subjects=("sub-P50001" "sub-P50002" "sub-P50003" "sub-P50004" "sub-P50005" "sub-P50006" "sub-P50007" "sub-P50008"
          "sub-P50009" "sub-P50010" "sub-P50011" "sub-P50012" "sub-P50013" "sub-P50014" "sub-P50014")
sessions=("ses-00" "ses-01" "ses-02" "ses-03")

# Loop through subjects
for subject in "${subjects[@]}"; do
    # Loop through sessions
    for session in "${sessions[@]}"; do
	# Define the file path to check
	file_to_check="/scratch/users/apines/SimStreams/${subject}_${session}_rs2_streamConnectivity_L_cantseeme.mat"
        if [ -e "$file_to_check" ]; then
            echo "File exists for $subject $session. Skipping inner loop."
        else
		# Submit sbatch job
        	sbatch sbatch_OpFl_ket.sh $subject $session
	fi
	done
done
