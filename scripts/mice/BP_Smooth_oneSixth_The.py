# note this needed to be ran in a notebook, initialized from a conda environment (testpines) through sherlock's web interface and some initial configuration and tweaking. Likely easier to try and get libraries loaded through ml (module load) if on sherlock
# update 5/3/24: it is not easier to try and get libraries to load w/ ml
# script DS_BP_Smooth.py
import os
import numpy as np
import scipy as sp
import h5py
import pandas as pd 
import matplotlib.pyplot as plt
from skimage.transform import downscale_local_mean, resize
from skimage.measure import block_reduce as downscale_local_mean

data_dir = '/scratch/users/apines/p50_mice/proc'
subjlist = pd.read_csv('/scratch/users/apines/p50_mice/demo/ipynb/sess_preDrug.txt', delimiter=' ', names=['date', 'session'])
fname = 'masked_dff_Gro_Masked_Sml.h5'
data_key = 'data'

# for gaussian smooth
from scipy.ndimage import gaussian_filter

# saveout  fp
output_suffix = '_BP_Smoothed_Sml_The'

# loop over every preDrugsession
for index, (date, sess) in subjlist.iterrows():
    full_fname = f"{data_dir}/{date}/{sess}/{fname}"
    # set new filename
    new_filename = f"{os.path.splitext(full_fname)[0]}{output_suffix}.h5"
    # load it, print out info
    with h5py.File(full_fname, 'r') as h5:
        print(sess, h5.keys())
        print(index)
        mask = h5['mask'][()]
        data = h5[data_key][()]
	# Theta frequency
        sos = sp.signal.butter(2, [4, 7.4], 'bandpass', fs=15, output='sos')
        data = sp.signal.sosfiltfilt(sos, data, axis=0)
	# gaus smooth
        data = gaussian_filter(data, sigma=(0, 2, 2))
        print(data.shape)
        # add gaussian kernel pass
        processed_data=data[:893,:,:]
        # save processed data
        with h5py.File(new_filename, 'w') as h5file:
            h5file.create_dataset('processed_data', data=processed_data)
            h5file.create_dataset('mask', data=mask)
            print(f"Saved processed data to {new_filename}")

