# only argument is the name of the folder to download
subjName=$1

module load system
module load rclone 

mkdir /oak/stanford/groups/leanew1/SHARED_DATASETS/private/WashU_psilocybin/${subjName}
mkdir /scratch/users/apines/PsiloData/${subjName}

rclone copy box:/psilo_pfm/${subjName} --include "*/*noGSR_sm4.dtseries.nii" /oak/stanford/groups/leanew1/SHARED_DATASETS/private/WashU_psilocybin/${subjName}

# you can use below to copy it from SHARED to scratch for proc.
rsync -av --relative /oak/stanford/groups/leanew1/SHARED_DATASETS/private/WashU_psilocybin/${subjName}/./*/func/*noGSR_sm4.dtseries.nii /scratch/users/apines/PsiloData/${subjName}/

