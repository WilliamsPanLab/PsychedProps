#!/bin/bash
#
#SBATCH --job-name=DS_des
#SBATCH --time=2:00:00
#SBATCH -n 1
#SBATCH --mem=20G
#SBATCH -p normal,leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
subj=$1

ml biology
ml workbench

/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/DS_surf_ts_mdma_fs5_rs.sh $subj ses-BL
