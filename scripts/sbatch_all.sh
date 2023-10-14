#!/bin/bash

subjects=("sub-MDMA002" "sub-MDMA003" "sub-MDMA005"
          "sub-MDMA007" "sub-MDMA008" "sub-MDMA009" "sub-MDMA011"
          "sub-MDMA012" "sub-MDMA013" "sub-MDMA014" "sub-MDMA015"
          "sub-MDMA016" "sub-MDMA017")
sessions=("ses-00" "ses-01" "ses-02" "ses-03")

# Loop through subjects
for subject in "${subjects[@]}"; do
    # Loop through sessions
    for session in "${sessions[@]}"; do
        # Submit sbatch job
        sbatch sbatch_OpFl.sh $subject $session
	done
done
