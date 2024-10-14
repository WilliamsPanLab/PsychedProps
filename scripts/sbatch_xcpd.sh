#!/bin/sh
#SBATCH -J xcp_d_p50
#SBATCH --array=1-17
#SBATCH --time=9:00:00
#SBATCH -n 2
#SBATCH --cpus-per-task=2
#SBATCH --mem=115G
#SBATCH -p leanew1,normal # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH -o logs/xcpd.%N.%j.out
#SBATCH -e logs/xcpd.%N.%j.err
#SBATCH --mail-user=xuezhang@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------
# Calculate the subject ID based on the job array index
subj=sub-MDMA$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
echo "Processing $subj..."
singularity run --cleanenv \
	-B /oak/stanford/groups/leanew1/users/apines:/home/xcp_d \
	-B /oak/stanford/groups/leanew1/users/apines/license.txt:/opt/freesurfer/license.txt \
	/oak/stanford/groups/leanew1/users/apines/xcp_d-0.3.0.simg -v /scratch/users/xuezhang/p50_copy/bids/data/derivatives/fmriprep-20.2.3/fmriprep -v /scratch/groups/leanew1/xcpd_outP50_36p_bp -w /scratch/groups/leanew1/xcpd_workP50_36p participant --participant_label $subj -p 36P --smoothing 6 --nthreads 4 --cifti
