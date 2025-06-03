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

subjects=("sub-MDMA001" "sub-MDMA002" "sub-MDMA003" "sub-MDMA005"
          "sub-MDMA007" "sub-MDMA008" "sub-MDMA009" "sub-MDMA011"
          "sub-MDMA012" "sub-MDMA013" "sub-MDMA014" "sub-MDMA015"
          "sub-MDMA016" "sub-MDMA017")
sessions=("ses-00" "ses-01" "ses-02" "ses-03")

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
   	    	sbatch sbatch_OpFl.sh $subject $session
  		#sbatch sbatch_OpFl_partial.sh $subject $session
  #	fi
  	done
done
