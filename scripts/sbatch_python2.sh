#!/bin/bash
#
#SBATCH --job-name=Python
#SBATCH --time=3:00:00
#SBATCH -n 1
#SBATCH --mem=24G
#SBATCH -p leanew1 # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=no@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
# export not working on non-interactive jobs?
#export PATH="/scratch/users/apines/miniconda3/bin:$PATH"
#conda activate testpines
ml python/3.9
# probably cant patch up enough libraries with just sherlock's modules
#ml py-numpy/1.26.3_py312
# above just breaks the conda env pathing
#python /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/mice/Group_Mask_and_DS_oneHalf.py
python3 /oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/mice/DS_Smooth_oneSixth_Drug_ExampleFrames.py
