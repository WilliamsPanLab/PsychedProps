#!/bin/bash
#
#SBATCH --job-name=sim_sub
#SBATCH --time=12:00:00
#SBATCH -n 1
#SBATCH --mem=5G
#SBATCH -p leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
# Number of times to submit the job

total_submissions=100
output_directory="/scratch/users/apines/SimStreams"
# Loop to submit the job
for ((i = 1; i <= total_submissions; i++)); do
    # test if output file exists
    output_file="$output_directory/${i}_streamConnectivity_L.mat"
    if [ -e "$output_file" ]; then
    	echo "Output file $output_file exists, skipping job submission for iteration $i"
    else
	echo "Submitting job with argument: $i"
        sbatch sbatch_OpFl_Simulated.sh $i
    
        if ((i < total_submissions)); then
            echo "Sleeping for 3 minutes before the next submission..."
            sleep 2  # Sleep for a bit
         else
            echo "All jobs submitted."
         fi
    fi
done
