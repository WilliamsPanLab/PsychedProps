#!/bin/bash
#
#SBATCH --job-name=simStreamComb
#SBATCH --time=20:50:00
#SBATCH -n 4
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

# aggregate info
matlab -nodisplay -r "SigTestSubjStreams('$1','$2')"
