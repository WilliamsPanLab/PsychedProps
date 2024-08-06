#!/bin/bash
base_dir="/scratch/users/apines/PsiloData"

# Loop through subject directories
for subject_dir in "$base_dir"/*/; do
    # Extract subject name from the directory path
    subject=$(basename "$subject_dir")

    # Loop through session directories for the current subject
    for session_dir in "$subject_dir"/*/; do
        # Extract session name from the directory path and remove the subject prefix
        session=$(basename "$session_dir" | sed "s/^$subject//")

        # Remove leading underscores or other separators
        session=${session##*_}

        # Check if the session name starts with "Baseline"
        if [[ $session == Baseline* ]]; then
            # Define the file path to check
            file_to_check="/scratch/users/apines/data/psil/${subject}/${session}/${subject}_${session}_rs2_k4_Prop_Feats_groCantSeeMe.csv"

            if [ -e "$file_to_check" ]; then
                echo "File exists for $subject $session. Skipping inner loop."
            else
                # Submit sbatch job
                echo $subject
                echo $session
                sbatch sbatch_DMNFC_map_psil.sh "$subject" "$session"
            fi
        fi
    done
done

