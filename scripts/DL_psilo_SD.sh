# only argument is the name of the folder to download
subjName=$1

module load system
module load rclone 

mkdir /scratch/users/apines/PsiloData/${subjName}
mkdir /oak/stanford/groups/leanew1/SHARED_DATASETS/private/WashU_psilocybin/${subjName}

rclone copy box:/psilo_pfm/${subjName} --include "*/*rsfMRI_sd.nii.gz" /oak/stanford/groups/leanew1/SHARED_DATASETS/private/WashU_psilocybin/${subjName}

# you can use below to copy it from SHARED to scratch for proc.
rsync -av --relative /oak/stanford/groups/leanew1/SHARED_DATASETS/private/WashU_psilocybin/${subjName}/./*/func/*rsfMRI_sd.nii.gz /scratch/users/apines/PsiloData/${subjName}/
