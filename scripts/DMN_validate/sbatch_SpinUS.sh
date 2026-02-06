#!/bin/bash
#
#SBATCH --job-name=OpFl
#SBATCH --time=40:00:00
#SBATCH -n 1
#SBATCH --mem=20G
#SBATCH -p leanew1,normal # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

./US_FPNSpins.sh
echo done with FPN
./US_VISSpins.sh
echo done with VIS
./US_MOTSpins.sh
echo done with MOT
