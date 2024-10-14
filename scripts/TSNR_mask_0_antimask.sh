ml biology
ml fsl
ml workbench

for i in $(seq -f "%03g" 1 17); do
  subj=sub-MDMA${i}

  # create subject output folder
  mkdir -p /scratch/users/xuezhang/data/mdma/${subj}/ses-00
  # use bet to make liberal brain mask for first resting state
  bet /oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/${subj}/ses-00/func/${subj}_ses-00_task-rs_acq-mb_dir-pe0_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe0_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_bet1.nii.gz -f .1
  # use fslmaths to binarize
  fslmaths /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe0_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_bet1.nii.gz -bin /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe0_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_bet1_bin.nii.gz
  # get inverse
  fslmaths /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe0_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_bet1_bin.nii.gz -mul -1 /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe0_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_bet1_bin-neg.nii.gz
  # add 1 so out-of-mask becomes in-mask
  fslmaths /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe0_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_bet1_bin-neg.nii.gz -add 1 /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe0_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_non-brain.nii.gz
  # make antimasked time series for SNR calculation
  fslmaths /oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/${subj}/ses-00/func/${subj}_ses-00_task-rs_acq-mb_dir-pe0_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz -mas /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe0_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_non-brain.nii.gz /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe0_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_antimasked.nii.gz

  # use bet to make liberal brain mask for second resting state
  bet /oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/${subj}/ses-00/func/${subj}_ses-00_task-rs_acq-mb_dir-pe1_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe1_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_bet1.nii.gz -f .1
  # use fslmaths to binarize
  fslmaths /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe1_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_bet1.nii.gz -bin /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe1_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_bet1_bin.nii.gz
  # get inverse
  fslmaths /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe1_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_bet1_bin.nii.gz -mul -1 /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe1_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_bet1_bin-neg.nii.gz
  # add 1 so out-of-mask becomes in-mask
  fslmaths /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe1_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_bet1_bin-neg.nii.gz -add 1 /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe1_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_non-brain.nii.gz
  # make antimasked time series for SNR calculation
  fslmaths /oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/${subj}/ses-00/func/${subj}_ses-00_task-rs_acq-mb_dir-pe1_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz -mas /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe1_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_non-brain.nii.gz /scratch/users/xuezhang/data/mdma/${subj}/ses-00/${subj}_ses-00_task-rs_acq-mb_dir-pe1_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_antimasked.nii.gz
done