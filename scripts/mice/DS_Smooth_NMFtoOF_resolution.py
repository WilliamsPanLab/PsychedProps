import os
import numpy as np
import scipy as sp
import h5py
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from skimage.transform import downscale_local_mean, resize
from skimage.measure import block_reduce as downscale_local_mean
from scipy.ndimage import gaussian_filter
from skimage.measure import block_reduce
from scipy.io import savemat

data_dir = '/scratch/users/apines/p50_mice/proc'
subjlist = pd.read_csv('/scratch/users/apines/p50_mice/demo/ipynb/sess_postLSD.txt', delimiter=' ', names=['date', 'session'])
fname = 'masked_dff.h5'
data_key = 'vid'

# note it should already be downsampled by factor of 2. So downscale factor is 3 rather than the 6 in OG -> OF downsample.
downscale = 3
ds_dims = (70, 67)

# load in hard-coded mask filepath from predrug session. One subject (mouse) because they all have the same mask at this point
fn='/scratch/users/apines/p50_mice/proc/20200228/thy1gc6s_0p3mgkg_m2000_preLSD0p3mgkg_1/masked_dff_DS_BP_Smoothed.h5'
# Open the HDF5 file
with h5py.File(fn, 'r') as f:
    # Read the dataset named '/mask'
    group_mask = f['/mask'][:]

# that's all fine and dandy, now lets load the .mat structure for the NMF solution
nmfFn='/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/init.mat'
nmfSol=sp.io.loadmat(nmfFn)
initV=nmfSol['initV']

# explicitly denote NMF components here
NMFcomponents = [1,10]

# initialize components matrix in output space: 70,67
components_matrix = np.zeros((ds_dims[0], ds_dims[1]))

# set original max: should be 1 for any given component
original_max=1

# for each DMN component, cut out <.3 loadings, downsample, pass gaussian over, and add to intiialized matrix
for component_index in NMFcomponents:
    # make a grid to populate components into xy space
    xygrid=np.zeros((208,200))
    # -1 to account for difference in python vs. matlab indexing
    component = initV[:, component_index - 1]
    # put it into a grid according to mask used	
    xygrid[group_mask]=component
    # equivalent thresholding as in human analyses
    xygrid[xygrid < 0.3] = 0
    # downsample the component: 3 is downsampling factor (6 from OG images)
    data = downscale_local_mean(xygrid, (3,3))
    # gaussian smooth
    data = gaussian_filter(data, sigma=(8, 8))
    # add to DMN component matrix
    components_matrix += data

# Calculate the maximum value of the downscaled data
downscaled_max = np.max(components_matrix)
# Calculate the scale factor
scale_factor = original_max / downscaled_max
# Apply the scale factor to the downscaled data
components_matrix = components_matrix * scale_factor

# equivalent downsampling to mask
original_max = np.max(group_mask)
downsampled_mask = downscale_local_mean(group_mask, downscale)
# Calculate the maximum value of the downscaled data
downscaled_max = np.max(downsampled_mask)
# Calculate the scale factor
scale_factor = original_max / downscaled_max
# Apply the scale factor to the downscaled data
downsampled_mask = downsampled_mask * scale_factor
# save out DMN component - image format
outputfilename='/home/users/apines/Mouse_DMN_DS.png';

# plotting code
plt.imshow(downsampled_mask, cmap='gray', interpolation='none', alpha=0.7)
plt.imshow(components_matrix, cmap='inferno', interpolation='none', alpha=0.7)
plt.axis('off')
plt.savefig('/home/users/apines/Mouse_DMN_DS.png')
plt.close()

# save out components matrix for calculating angular distance
# Define the path for the .mat file
mat_filename = '/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat'

# Create a dictionary with the data to save
mat_data = {'components_matrix': components_matrix, 'mask': downsampled_mask}

# Save the dictionary to a .mat file
savemat(mat_filename, mat_data)
print(f"Components matrix saved to {mat_filename}")
