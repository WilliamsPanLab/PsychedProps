#!/bin/bash
#
#SBATCH -J DS_surf
#SBATCH --time=8:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH -p leanew1,normal  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH -o logs/DS_surf.%N.%j.out
#SBATCH -e logs/DS_surf.%N.%j.err
#SBATCH --mail-user=xuezhang@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------


cd /scratch/users/xuezhang
./DS_surf_ts_mdma_fs5.sh


exitcode=$?
exit $exitcode