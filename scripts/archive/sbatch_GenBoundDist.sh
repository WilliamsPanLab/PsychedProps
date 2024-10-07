#!/bin/bash
#
#SBATCH --job-name=GenBoundDist
#SBATCH --time=1:00:00
#SBATCH -n 1
#SBATCH --mem=16G
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

subj=$1
SolutionNumber=$2
DMNnumber=$3
merge=$4
MergeVec1=$5
MergeVec2=$6

# combine streamlines
matlab -nodisplay -r "GenBoundDistMap('$subj','$SolutionNumber','$DMNnumber','$merge',['$MergeVec1' ' ' '$MergeVec2'])"
