# set freesurfer dir
export FREESURFER_HOME=/share/software/user/open/freesurfer/6.0.0
# set subjs dir
export SUBJECTS_DIR=/share/software/user/open/freesurfer/6.0.0/subjects
# set freesurfer license
export FS_LICENSE=/oak/stanford/groups/leanew1/users/apines/license.txt

# convert dlabel to giftis
wb_command -cifti-separate ~/null_lL_WG33/Gordon333_FreesurferSubcortical.32k_fs_LR.dlabel.nii COLUMN -label CORTEX_LEFT /oak/stanford/groups/leanew1/users/apines/maps/Gordon_333_left.label.gii

wb_command -cifti-separate ~/null_lL_WG33/Gordon333_FreesurferSubcortical.32k_fs_LR.dlabel.nii COLUMN -label CORTEX_RIGHT /oak/stanford/groups/leanew1/users/apines/maps/Gordon_333_right.label.gii


