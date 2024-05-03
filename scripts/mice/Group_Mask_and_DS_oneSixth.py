import os
import numpy as np
import scipy as sp
import h5py
import pandas as pd 
import matplotlib.pyplot as plt
from skimage.transform import downscale_local_mean, resize
from sklearn.decomposition import NMF
from skimage.measure import block_reduce as downscale_local_mean

data_dir = '/scratch/users/apines/p50_mice/proc'
subjlist = pd.read_csv('/scratch/users/apines/p50_mice/demo/ipynb/sess_preDrug.txt', delimiter=' ', names=['date', 'session'])
fname = 'masked_dff.h5'
data_key = 'vid'

downscale = 6
ds_dims = (70, 67)
group_mask = np.ones(ds_dims, dtype=bool)
fname_demo = '/scratch/users/apines/p50_mice/demo/ipynbdemo_data_v35.npz'

# split up into a loop to just create group_mask

# initialize array
if not os.path.isfile(fname_demo):
    #all_data = []
    for index, (date, sess) in subjlist.iterrows():
        full_fname = f"{data_dir}/{date}/{sess}/{fname}"
        with h5py.File(full_fname, 'r') as h5:
            print(sess, h5.keys())
            print(index)
            mask = h5['atlas'][()] > 0
            data = h5[data_key][()]
            
            if downscale is not None:
                mask = downscale_local_mean(mask.astype(float), 2*(downscale,)) > 0.5
                data = downscale_local_mean(data, (1,) + 2*(downscale,))
                mask = mask[:ds_dims[0], :ds_dims[1]]
                data = data[:, :ds_dims[0], :ds_dims[1]]

            group_mask &= mask

# saveout 
output_suffix = '_Gro_Masked_Sml'

# Iterate over your data array, assuming 'all_data' is a list of numpy arrays
for index, (date, sess) in subjlist.iterrows():
    full_fname = f"{data_dir}/{date}/{sess}/{fname}"
    with h5py.File(full_fname, 'r') as h5:
        data = h5[data_key][()]
        # apply equiv downsample
        data = downscale_local_mean(data, (1,) + 2*(downscale,))
        data = data[:, :ds_dims[0], :ds_dims[1]]
        # Save the processed data
        new_filename = f"{os.path.splitext(full_fname)[0]}{output_suffix}.h5"
        with h5py.File(new_filename, 'w') as h5file:
            h5file.create_dataset('data', data=data)
            h5file.create_dataset('mask', data=group_mask)
            print(f"Saved Masked data to {new_filename}")

