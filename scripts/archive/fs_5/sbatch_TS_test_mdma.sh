#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=9:00:00
#SBATCH -n 1
#SBATCH --mem=125G
#SBATCH -p leanew1,normal  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=no@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

# will need matlab
module load matlab
# and freesurfer
module load biology
module load freesurfer/7.3.2
# and workbench
module load workbench
matlab -nodisplay -r "WithinSubj_ts_Conjuction"

