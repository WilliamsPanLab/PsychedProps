import sys
import nibabel as nb
import numpy as np

subj=sys.argv[1]

# images parent folder
parentFP='/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/'

# find ses-00 raw rest image (rs1)
rs1fp=parentFP + subj + '/ses-00/func/' + subj + '_ses-00_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_bold.dtseries.nii'
# extract mean grayordinate-wise signal
rs1=nb.load(rs1fp)
rs1_data=rs1.get_fdata()
meanrs1=np.mean(rs1_data,axis=0)
# get nifti version to extract values outside of brain
rs1niffp='/scratch/users/apines/data/mdma/' + subj + '/ses-00/' + subj + '_ses-00_task-rs_acq-mb_dir-pe0_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_antimasked.nii.gz'
rs1nif=nb.load(rs1niffp).get_fdata()
# extract meanSD
sdrs1=np.std(rs1nif,axis=3)
mean_sdsrs1=np.mean(sdrs1)
# caclulate tsnr
tsnr_rs1=meanrs1/mean_sdsrs1

# find ses-00 raw rest image (rs2)
rs2fp=parentFP + subj + '/ses-00/func/' + subj + '_ses-00_task-rs_acq-mb_dir-pe1_run-0_space-fsLR_den-91k_bold.dtseries.nii'
# extract mean grayordinate-wise signal
rs2=nb.load(rs2fp)
rs2_data=rs2.get_fdata()
meanrs2=np.mean(rs2_data,axis=0)
# get nifti version to extract values outside of brain
rs2niffp='/scratch/users/apines/data/mdma/' + subj + '/ses-00/' + subj + '_ses-00_task-rs_acq-mb_dir-pe1_run-0_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold_antimasked.nii.gz'
rs2nif=nb.load(rs2niffp).get_fdata()
# extract meanSD
sdrs2=np.std(rs2nif,axis=3)
mean_sdsrs2=np.mean(sdrs2)
# caclulate tsnr
tsnr_rs2=meanrs2/mean_sdsrs2

# average tsnr across rs maps
tsnr_rs=(tsnr_rs1+tsnr_rs2)/2

# load in a dscalar to use as reference image for creating a new cifti
refCif=nb.load('/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients.dscalar.nii')
# dumb matrix combination just to match dimensions, dimensions 2:26 are extraneous
tsnr_datamat=np.zeros((26,91282))
tsnr_datamat[0,]=tsnr_rs
# create new cifti to save out
cifti_img = nb.cifti2.Cifti2Image(tsnr_datamat, header=refCif.header)

# save out to scratch
output_cifti_fp = '/scratch/users/apines/data/mdma/' + subj + '/ses-00/' + subj + 'meanTSNR.dscalar.nii'

# Save the CIFTI image to the specified file path
nb.save(cifti_img, output_cifti_fp)