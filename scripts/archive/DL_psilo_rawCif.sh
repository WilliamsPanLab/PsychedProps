# only argument is the name of the folder to download
subjName=$1

module load system
module load rclone 

mkdir /scratch/users/apines/PsiloData/${subjName}

rclone copy box:/psilo_pfm/${subjName} --include "*/*upck_faln_dbnd_xr3d_dc_atl.ctx.dtseries.nii" /scratch/users/apines/PsiloData/${subjName}
