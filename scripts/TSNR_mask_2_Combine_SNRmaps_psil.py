import numpy as np
import nibabel as nb

# List of subjects: note missing raw cifti for PS19
subjects = ['PS03', 'PS16', 'PS18', 'PS21','PS24']

# Initialize an array to accumulate tsnr_rs values
total_tsnr_rs = np.zeros((1, 91206))

# Iterate over each subject
for subj in subjects:
    # Load the saved CIFTI file
    cifti_file_path = '/scratch/users/apines/PsiloData/' + subj + '/'+ subj + '_Baseline1/func/' + subj + '_meanTSNR.dtseries.nii'
    cifti_img = nb.load(cifti_file_path)

    # Extract tsnr_rs values and accumulate
    tsnr_rs_subject = cifti_img.get_fdata()[0, :]
    total_tsnr_rs += tsnr_rs_subject

# Calculate the average tsnr_rs across all subjects
average_tsnr_rs = total_tsnr_rs / len(subjects)

# use last-loaded image as reference cifti
refCif=cifti_img
# dumb matrix combination just to match dimensions, dimensions 2:26 are extraneous
tsnr_datamat=np.zeros((1026,91206))
tsnr_datamat[0,]=average_tsnr_rs
# create new cifti to save out
a_cifti_img = nb.cifti2.Cifti2Image(tsnr_datamat, header=refCif.header)

# Save the averaged CIFTI image
output_average_cifti_fp = '/oak/stanford/groups/leanew1/users/apines/maps/average_tsnr_rs_psil.dtseries.nii'
nb.save(a_cifti_img, output_average_cifti_fp)

