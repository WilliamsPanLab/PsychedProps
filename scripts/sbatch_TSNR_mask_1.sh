#!/bin/bash
#
#SBATCH -J tSNR_1
#SBATCH --time=8:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH -p leanew1,normal  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH -o logs/tSNR_1.%N.%j.out
#SBATCH -e logs/tSNR_1.%N.%j.err
#SBATCH --mail-user=xuezhang@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

ml python/3.12.1
python3 /scratch/users/xuezhang/TSNR_mask_1_ExtractSNR_subjectwise.py


exitcode=$?
exit $exitcode

