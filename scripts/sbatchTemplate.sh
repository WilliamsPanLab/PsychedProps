#!/bin/bash
#
#SBATCH --job-name=Matlab
#SBATCH --time=4:00:00
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH -p normal,leanew1 # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=no@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
module load matlab
