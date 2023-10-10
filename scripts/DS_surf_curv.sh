# set freesurfer dir
export FREESURFER_HOME=/share/software/user/open/freesurfer/6.0.0
# set subjs dir
export SUBJECTS_DIR=/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/freesurfer
# set freesurfer license
export FS_LICENSE=/oak/stanford/groups/leanew1/users/apines/license.txt

childfp=/scratch/users/apines/data/mdma/${1}

# left
mri_surf2surf --srcsubject $1 --trgsubject fsaverage4 --hemi lh --sval /oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/freesurfer/${1}/surf/lh.curv.pial --tval ${childfp}/${1}_fs4space_curv.L.func.gii
# right
mri_surf2surf --srcsubject $1 --trgsubject fsaverage4 --hemi rh --sval /oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/freesurfer/${1}/surf/rh.curv.pial --tval ${childfp}/${1}_fs4space_curv.R.func.gii
