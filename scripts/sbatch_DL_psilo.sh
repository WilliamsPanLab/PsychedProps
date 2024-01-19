#!/bin/bash
#
#SBATCH --job-name=rclone
#SBATCH --time=4:00:00
#SBATCH -n 4
#SBATCH --mem=30G
#SBATCH -p leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

subj=$1

/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts/DL_psilo.sh $subj
