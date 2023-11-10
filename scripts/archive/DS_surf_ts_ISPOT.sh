subj=$1
# set freesurfer dir
export FREESURFER_HOME=/share/software/user/open/freesurfer/6.0.0
# set subjs dir
export SUBJECTS_DIR=/share/software/user/open/freesurfer/6.0.0/subjects
# set freesurfer license
export FS_LICENSE=/oak/stanford/groups/leanew1/users/apines/license.txt
# ispot fp
parentfp=/oak/stanford/groups/leanew1/SHARED_DATASETS/private/ispot/bids/data/derivatives/fmriprep-20.2.1/fmriprep/
childfp=/scratch/users/apines/data/ispot

# make output dirs
mkdir ${childfp}/${1}

# left gng
/share/software/user/open/freesurfer/6.0.0/bin/mri_vol2surf --src ${parentfp}/${1}/ses-00/func/${1}_ses-00_task-gonogo_acq-sb_dir-pe1_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii.gz --out ${childfp}/${1}/gng_lh.func.gii --mni152reg --hemi lh --trgsubject fsaverage4
# right gng
/share/software/user/open/freesurfer/6.0.0/bin/mri_vol2surf --src ${parentfp}/${1}/ses-00/func/${1}_ses-00_task-gonogo_acq-sb_dir-pe1_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii.gz --out ${childfp}/${1}/gng_rh.func.gii --mni152reg --hemi rh --trgsubject fsaverage4
# convert them all to mgh
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${1}/gng_lh.func.gii ${childfp}/${1}/gng_lh.mgh
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${1}/gng_rh.func.gii ${childfp}/${1}/gng_rh.mgh

# left con faces
/share/software/user/open/freesurfer/6.0.0/bin/mri_vol2surf --src ${parentfp}/${1}/ses-00/func/${1}_ses-00_task-conscious_acq-sb_dir-pe1_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii.gz --out ${childfp}/${1}/con_lh.func.gii --mni152reg --hemi lh --trgsubject fsaverage4
# right con faces
/share/software/user/open/freesurfer/6.0.0/bin/mri_vol2surf --src ${parentfp}/${1}/ses-00/func/${1}_ses-00_task-conscious_acq-sb_dir-pe1_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii.gz --out ${childfp}/${1}/con_rh.func.gii --mni152reg --hemi rh --trgsubject fsaverage4
# convert them to mgh
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${1}/con_lh.func.gii ${childfp}/${1}/con_lh.mgh
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${1}/con_rh.func.gii ${childfp}/${1}/con_rh.mgh

# left noncon
/share/software/user/open/freesurfer/6.0.0/bin/mri_vol2surf --src ${parentfp}/${1}/ses-00/func/${1}_ses-00_task-nonconscious_acq-sb_dir-pe1_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii.gz --out ${childfp}/${1}/noncon_lh.func.gii --mni152reg --hemi lh --trgsubject fsaverage4
# right noncon
/share/software/user/open/freesurfer/6.0.0/bin/mri_vol2surf --src ${parentfp}/${1}/ses-00/func/${1}_ses-00_task-nonconscious_acq-sb_dir-pe1_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii.gz --out ${childfp}/${1}/noncon_rh.func.gii --mni152reg --hemi rh --trgsubject fsaverage4
# convert them to mgh
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${1}/noncon_lh.func.gii ${childfp}/${1}/noncon_lh.mgh
/share/software/user/open/freesurfer/6.0.0/bin/mri_convert.bin ${childfp}/${1}/noncon_rh.func.gii ${childfp}/${1}/noncon_rh.mgh

