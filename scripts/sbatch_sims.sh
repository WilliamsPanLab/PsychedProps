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
sleep 10001
# Loop to submit the job
for ((i = 1; i <= total_submissions; i++)); do
    echo "Submitting job with argument: $i"
    sbatch sbatch_OpFl_Simulated.sh $i
    
    if ((i < total_submissions)); then
        echo "Sleeping for 3 minutes before the next submission..."
        sleep 280  # Sleep for a bit
    else
        echo "All jobs submitted."
    fi
done
