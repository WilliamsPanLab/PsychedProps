import scipy.io as sio
import numpy as np
from skimage.transform import resize
import h5py

# Input file path (MAT file) and output file path
# will have to make it subject specific
input_mat_file='/oak/stanford/groups/leanew1/users/apines/p50_mice/proc2/proc/20200228/thy1gc6s_0p3mgkg_m2000_postLSD0p3mgkg_0/atlasAlignedToData.mat'
output_mat_file = '/oak/stanford/groups/leanew1/users/apines/p50_mice/proc2/proc/20200228/thy1gc6s_0p3mgkg_m2000_postLSD0p3mgkg_0/test_downsampled_atlas.mat'

# Load the .mat file
mat_data = sio.loadmat(input_mat_file)

# extract values
atlas = mat_data['atlasAlignedToData']

# Set downscale factor
downscale_factor = 6  # Change this to your desired factor

# Downsample using nearest-neighbor interpolation (order=0)
downsampled_atlas = resize(
    atlas,
    (atlas.shape[0] // downscale_factor, atlas.shape[1] // downscale_factor),
    order=0,  # Nearest-neighbor to preserve labels
    anti_aliasing=False,
    preserve_range=True
).astype(np.int32)

# Save the downsampled atlas to a new .mat file
sio.savemat(output_mat_file, {'downsampled_atlas': downsampled_atlas})

print(f"Downsampled atlas saved to {output_mat_file}")

