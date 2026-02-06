#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=30:00:00
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH -p leanew1,normal # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

./DS_MOTSpins.sh
