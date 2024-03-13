# ketamine subjects acquired before the scanner upgrade have some negative bold signal in raw images in a few select voxels (mostly inferior cerebellum)
# this script 0's them out for confound regression. These voxels are not included in final analyses

import nibabel as nib
import glob
import sys
import os
import numpy as np

def process_cifti_files(folder_path, extension='.dtseries.nii'):
    """
    Process all CIFTI files in the specified folder, replacing negative values with zeros.
    Args:
        folder_path (str): Path to the folder containing CIFTI files.
        extension (str): Extension of the CIFTI files to process.
    """

    # Construct the pattern to match files
    search_pattern = os.path.join(folder_path, f'*acq-mb*{extension}')
    # Find all files in the folder matching the pattern
    cifti_files = glob.glob(search_pattern)

    for cifti_file in cifti_files:
        print(f'Processing {cifti_file}...')
        
        # Load the CIFTI file
        img = nib.load(cifti_file)

        # Get data from the CIFTI file
        data = img.get_fdata()

        # Replace negative values with zeros
        data[data < 0] = 0

        # and nans, 0'ed out in post subj 7 mb as well
        data[np.isnan(data)]=0

        # Create a new CIFTI image with the modified data
        new_img = nib.Cifti2Image(data, header=img.header, nifti_header=img.nifti_header)

        # Generate the output file path
        output_file = cifti_file

        # Save the new CIFTI file
        nib.save(new_img, output_file)
        print(f'Saved modified file to {output_file}')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python script.py <folder_path>")
    else:
        folder_path = sys.argv[1]
        process_cifti_files(folder_path)

