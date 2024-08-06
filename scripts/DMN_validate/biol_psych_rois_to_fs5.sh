ml biology
ml freesurfer
# AG L
/share/software/user/open/freesurfer/6.0.0/bin/mri_vol2surf --src /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/rois/Defaultmode_AG_L.nii.gz --out /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/rois/Defaultmode_AG_L.func.gii --mni152reg --hemi lh --trgsubject fsaverage5
# AG R
/share/software/user/open/freesurfer/6.0.0/bin/mri_vol2surf --src /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/rois/Defaultmode_AG_R.nii.gz --out /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/rois/Defaultmode_AG_R.func.gii --mni152reg --hemi rh --trgsubject fsaverage5
# medial has to go to both hemis
# amPFC
/share/software/user/open/freesurfer/6.0.0/bin/mri_vol2surf --src /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/rois/Defaultmode_amPFC_M.nii.gz --out /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/rois/Defaultmode_amPFC_M_L.func.gii --mni152reg --hemi lh --trgsubject fsaverage5
/share/software/user/open/freesurfer/6.0.0/bin/mri_vol2surf --src /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/rois/Defaultmode_amPFC_M.nii.gz --out /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/rois/Defaultmode_amPFC_M_R.func.gii --mni152reg --hemi rh --trgsubject fsaverage5
# pcc
/share/software/user/open/freesurfer/6.0.0/bin/mri_vol2surf --src /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/rois/Defaultmode_PCC_M.nii.gz --out /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/rois/Defaultmode_PCC_M_L.func.gii --mni152reg --hemi lh --trgsubject fsaverage5
/share/software/user/open/freesurfer/6.0.0/bin/mri_vol2surf --src /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/rois/Defaultmode_PCC_M.nii.gz --out /oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/rois/Defaultmode_PCC_M_R.func.gii --mni152reg --hemi rh --trgsubject fsaverage5
