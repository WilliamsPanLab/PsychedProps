#!/bin/bash
#
#SBATCH --job-name=simStreamComb
#SBATCH --time=1:50:00
#SBATCH -n 1
#SBATCH --mem=25G
#SBATCH -p leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

# will need matlab
module load matlab
# and freesurfer
module load biology
module load freesurfer/7.3.2
# and workbench
module load workbench

# combine streamlines
matlab -nodisplay -r "SigTestStreams('$1','$2')"
