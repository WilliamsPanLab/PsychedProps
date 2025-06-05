#!/bin/bash
#
#SBATCH --job-name=sbatchEm
#SBATCH --time=6:50:00
#SBATCH -n 1
#SBATCH --mem=12G
#SBATCH -p leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

subjects=($(seq -f "S%02g" 1 20))
sessions=("LSD" "PCB")

# Loop through subjects
for subject in "${subjects[@]}"; do
    # Loop through sessions
    for session in "${sessions[@]}"; do
	# Define the file path to check
  #	file_to_check="/scratch/users/apines/SimStreams/${subject}_${session}_rs2_streamConnectivity_L_cantseeme.mat"
  #      if [ -e "$file_to_check" ]; then
  #          echo "File exists for $subject $session. Skipping inner loop."
  #      else
		# Submit sbatch job
		echo $subject
		echo $session
   	    	sbatch sbatch_OpFl_lsd.sh $subject $session
  #	fi
  	done
done
