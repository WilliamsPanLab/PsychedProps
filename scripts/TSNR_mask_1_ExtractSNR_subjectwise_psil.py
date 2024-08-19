import sys
import nibabel as nb
import numpy as np

subj=sys.argv[1]

# images parent folder
parentFP='/scratch/users/apines/PsiloData/' + subj + '/'+ subj + '_Baseline1/func/'

# find basline raw rest image
rs1fp=parentFP + subj + '_Baseline1_upck_faln_dbnd_xr3d_dc_atl.ctx.dtseries.nii'
# extract mean grayordinate-wise signal
rs1=nb.load(rs1fp)
rs1_data=rs1.get_fdata()
meanrs1=np.mean(rs1_data,axis=0)

# load in sd map
sdfp=parentFP + subj + '_Baseline1_rsfMRI_sd_antimasked.nii.gz'
sdrs1=nb.load(sdfp).get_fdata()
# extract mean SD
mean_sdsrs1=np.mean(sdrs1)
# caclulate tsnr
tsnr_rs=meanrs1/mean_sdsrs1

# load in a dscalar to use as reference image for creating a new cifti
refCif=nb.load('/home/users/apines/Template_Cifti.dtseries.nii')
# dumb matrix combination just to match dimensions, dimensions 2:26 are extraneous
tsnr_datamat=np.zeros((1026,91206))
tsnr_datamat[0,]=tsnr_rs
# create new cifti to save out
cifti_img = nb.cifti2.Cifti2Image(tsnr_datamat, header=refCif.header)

# save out to scratch
output_cifti_fp = parentFP + subj + '_meanTSNR.dtseries.nii'
# Save the CIFTI image to the specified file path
nb.save(cifti_img, output_cifti_fp)
