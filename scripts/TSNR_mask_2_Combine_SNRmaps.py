import numpy as np
import nibabel as nb

# List of subjects
subjects = ['sub-MDMA001', 'sub-MDMA002', 'sub-MDMA003', 'sub-MDMA005',
            'sub-MDMA007', 'sub-MDMA008', 'sub-MDMA009',
            'sub-MDMA011', 'sub-MDMA012', 'sub-MDMA013', 'sub-MDMA014', 'sub-MDMA015',
            'sub-MDMA016', 'sub-MDMA017']

# Initialize an array to accumulate tsnr_rs values
total_tsnr_rs = np.zeros((1, 91282))

# Iterate over each subject
for subj in subjects:
    # Load the saved CIFTI file
    cifti_file_path = '/scratch/users/apines/data/mdma/' + subj + '/ses-00/' + subj + 'meanTSNR.dscalar.nii'
    cifti_img = nb.load(cifti_file_path)

    # Extract tsnr_rs values and accumulate
    tsnr_rs_subject = cifti_img.get_fdata()[0, :]
    total_tsnr_rs += tsnr_rs_subject

# Calculate the average tsnr_rs across all subjects
average_tsnr_rs = total_tsnr_rs / len(subjects)

# use last-loaded image as reference cifti
refCif=cifti_img
# dumb matrix combination just to match dimensions, dimensions 2:26 are extraneous
tsnr_datamat=np.zeros((26,91282))
tsnr_datamat[0,]=average_tsnr_rs
# create new cifti to save out
a_cifti_img = nb.cifti2.Cifti2Image(tsnr_datamat, header=refCif.header)

# Save the averaged CIFTI image
output_average_cifti_fp = '/oak/stanford/groups/leanew1/users/apines/maps/average_tsnr_rs_mdma.dscalar.nii'
nb.save(a_cifti_img, output_average_cifti_fp)

