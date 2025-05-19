import os
import numpy as np
import scipy as sp
import h5py
import pandas as pd
from scipy.ndimage import gaussian_filter

# Data directory and file paths
data_dir = '/oak/stanford/groups/leanew1/users/apines/p50_mice/proc2/proc/'
subjlist = pd.read_csv('/scratch/users/apines/p50_mice/demo/ipynb/sess_preDrug.txt', delimiter=' ', names=['date', 'session'])
fname = 'masked_dff.h5'
data_key = 'vid'

# New saveout suffix for full filtered & smoothed data
output_suffix_filtered = '_Full_BP_Smoothed_NoDS'

# Load the mask
fn = '/oak/stanford/groups/leanew1/users/apines/p50_mice/proc2/proc/20200228/thy1gc6s_0p3mgkg_m2000_preLSD0p3mgkg_1/masked_dff_Gro_Masked_Sml.h5'
with h5py.File(fn, 'r') as f:
    group_mask = f['/mask'][:]

# Iterate over each subject/session
for index, (date, sess) in subjlist.iterrows():
    full_fname = f"{data_dir}/{date}/{sess}/{fname}"
    with h5py.File(full_fname, 'r') as h5:
        data = h5[data_key][()]  # Load full data (no downsampling)

        # Bandpass filter (1-4 Hz)
        sos = sp.signal.butter(2, [1, 4], 'bandpass', fs=15, output='sos')
        data = sp.signal.sosfiltfilt(sos, data, axis=0)

        # Gaussian smoothing
        data = gaussian_filter(data, sigma=(0, 2, 2))

        # Save the full filtered & smoothed data (no downsampling)
        new_filename_filtered = f"{os.path.splitext(full_fname)[0]}{output_suffix_filtered}.h5"
        with h5py.File(new_filename_filtered, 'w') as h5file:
            h5file.create_dataset('processed_data', data=data)
            h5file.create_dataset('mask', data=group_mask)
            print(f"Saved full BP & smoothed data to {new_filename_filtered}")

