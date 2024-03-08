# only argument is the name of the folder to download
subjName=$1

module load system
module load rclone 

mkdir /scratch/users/apines/PsiloData/${subjName}

rclone copy box:/psilo_pfm/${subjName} --include "*/*rsfMRI_sd.nii.gz" /scratch/users/apines/PsiloData/${subjName}
