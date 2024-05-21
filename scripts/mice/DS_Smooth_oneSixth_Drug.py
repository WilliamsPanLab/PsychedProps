import os
import numpy as np
import scipy as sp
import h5py
import pandas as pd
import matplotlib.pyplot as plt
from skimage.transform import downscale_local_mean, resize
from skimage.measure import block_reduce as downscale_local_mean
from scipy.ndimage import gaussian_filter

data_dir = '/scratch/users/apines/p50_mice/proc'
subjlist = pd.read_csv('/scratch/users/apines/p50_mice/demo/ipynb/sess_preLSD.txt', delimiter=' ', names=['date', 'session'])
fname = 'masked_dff.h5'
data_key = 'vid'

# set downscale factors
downscale = 6
ds_dims = (70, 67)

# saveout - change FN to match post-bp gaus (IS)
output_suffix = '_Gro_Masked_Sml_BP_Smoothed_Sml'

# load in hard-coded mask filepath from predrug session. One subject (mouse) because they all have the same mask at this point
fn='/scratch/users/apines/p50_mice/proc/20200228/thy1gc6s_0p3mgkg_m2000_preLSD0p3mgkg_1/masked_dff_DS_BP_Smoothed.h5'
# Open the HDF5 file
with h5py.File(fn, 'r') as f:
    # Read the dataset named '/mask'
    group_mask = f['/mask'][:]

# Iterate over your data array, assuming 'all_data' is a list of numpy arrays
for index, (date, sess) in subjlist.iterrows():
    full_fname = f"{data_dir}/{date}/{sess}/{fname}"
    with h5py.File(full_fname, 'r') as h5:
        data = h5[data_key][()]
        # apply equiv downsample
        data = downscale_local_mean(data, (1,) + 2*(downscale,))
        data = data[:, :ds_dims[0], :ds_dims[1]]
        # bandpass
        sos = sp.signal.butter(2, [0.5, 4], 'bandpass', fs=15, output='sos')
        data = sp.signal.sosfiltfilt(sos, data, axis=0)
        # gaussian
        data = gaussian_filter(data, sigma=(0, 2, 2))
        # 893 frames to use lowest common denom
        data=data[:893,:,:]
	# Save the processed data
        new_filename = f"{os.path.splitext(full_fname)[0]}{output_suffix}.h5"
        with h5py.File(new_filename, 'w') as h5file:
            h5file.create_dataset('processed_data', data=data)
            h5file.create_dataset('mask', data=group_mask)
            print(f"Saved Masked data to {new_filename}")



