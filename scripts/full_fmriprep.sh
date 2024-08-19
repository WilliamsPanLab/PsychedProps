#!/bin/bash
#
#SBATCH -J fprep
#SBATCH --time=48:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=110G
#SBATCH -p normal,leanew1  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH -o logs/fprep.%N.%j.out
#SBATCH -e logs/fprep.%N.%j.err
#SBATCH --mail-user=apines@stanford.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

#set -x

# NOTE: logs folder must exist!
#STUDY="/scratch/users/apines/data/engage/data"
STUDY="/oak/stanford/groups/leanew1/SHARED_DATASETS/private/connectome/bids/data"
BIDS_DIR="${STUDY}/bids-input"
DERIVS_DIR="${STUDY}/derivatives/fmriprep-20.2.3"
WF_DIR="${SCRATCH}/tmp"

containers_home="/oak/stanford/groups/leanew1/containers/poldrack"
fmriprep_container="fmriprep-20.2.3.simg"

# Prepare some writeable bind-mount points.
TEMPLATEFLOW_HOST_HOME=$HOME/.cache/templateflow
FMRIPREP_HOST_CACHE=$HOME/.cache/fmriprep
mkdir -p ${TEMPLATEFLOW_HOST_HOME}
mkdir -p ${FMRIPREP_HOST_CACHE}

# Prepare derivatives folder
mkdir -p ${BIDS_DIR}
mkdir -p ${DERIVS_DIR}
mkdir -p ${WF_DIR}

# Make sure FS_LICENSE is defined in the container.
export SINGULARITYENV_FS_LICENSE=/oak/stanford/groups/leanew1/users/apines/license.txt

# Designate a templateflow bind-mount point
export SINGULARITYENV_TEMPLATEFLOW_HOME="/templateflow"
SINGULARITY_CMD="singularity run --cleanenv \
-B ${BIDS_DIR}:/data \
-B ${DERIVS_DIR}:/derivatives \
-B ${TEMPLATEFLOW_HOST_HOME}:${SINGULARITYENV_TEMPLATEFLOW_HOME} \
-B ${WF_DIR}:/tmp \
${containers_home}/${fmriprep_container}"

opts="-w /tmp/tmp_workflow \
-vvv \
--nthreads 16 \
--mem_mb 30000 \
--output-spaces MNI152NLin2009cAsym:res-2 T1w fsnative \
--use-aroma \
--cifti-output \
--skip_bids_validation \
--resource-monitor \
--notrack"

# Parse the participants.tsv file and extract one subject ID from the line corresponding to this SLURM task.
# subject=$(head -n 1 ${BIDS_DIR}/participants.tsv)
# subject=$( sed -n -E "$((${SLURM_ARRAY_TASK_ID}))s/sub-(\S*)\>.*/\1/gp" ${BIDS_DIR}/participants.tsv )
# AP 9-20: just input subj as arg since we are running one at a time anyways
subject=$1
echo $subject

# Remove IsRunning files from FreeSurfer
#find ${LOCAL_FREESURFER_DIR}/sub-$subject/ -name "*IsRunning*" -type f -delete

# Compose the command line
cmd="${SINGULARITY_CMD} ${BIDS_DIR} ${DERIVS_DIR} participant --participant-label ${subject} ${opts}"

# Setup done, run the command
echo Running task ${SLURM_ARRAY_TASK_ID}
echo Commandline: $cmd
eval $cmd
exitcode=$?

# Output results to a table
echo "sub-$subject   ${SLURM_ARRAY_TASK_ID}    $exitcode" \
      >> logs/${SLURM_JOB_NAME}.${SLURM_ARRAY_JOB_ID}.tsv
echo Finished tasks ${SLURM_ARRAY_TASK_ID} with exit code $exitcode
exit $exitcode

