import os
from skimage.measure import block_reduce
import numpy as np
import scipy as sp
import h5py
import pandas as pd 
import matplotlib
import matplotlib.pyplot as plt
from skimage.transform import downscale_local_mean, resize
from skimage.measure import block_reduce as downscale_local_mean
from scipy.stats import mode
matplotlib.use('Agg')
import matplotlib.colors as mcolors


data_dir = '/oak/stanford/groups/leanew1/users/apines/p50_mice/proc2/proc/'
#subjlist = pd.read_csv('/scratch/users/apines/p50_mice/demo/ipynb/sess_preDrug.txt', delimiter=' ', names=['date', 'session'])
subjlist = pd.read_csv('/scratch/users/apines/p50_mice/demo/ipynb/sess_postDrug.txt', delimiter=' ', names=['date', 'session'])
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



# now we just need to take care of downsampling the atlas with a non-interpolative approach
def downsample_atlas_nearest(atlas, downscale):
    block_size = (downscale, downscale)
    return block_reduce(atlas, block_size, func=np.median).astype(np.int32)

# and make a group composite for reference
def create_group_composite_atlas(subjlist, data_dir, fname, group_mask, downscale, ds_dims):
    atlases = []
    for index, (date, sess) in subjlist.iterrows():
        full_fname = f"{data_dir}/{date}/{sess}/{fname}"
        with h5py.File(full_fname, 'r') as h5:
            atlas = h5['atlas'][()]
            downsampled_atlas = downsample_atlas_nearest(atlas, downscale)[:ds_dims[0], :ds_dims[1]]
            atlases.append(downsampled_atlas)
    stacked_atlases = np.stack(atlases, axis=-1)
    group_atlas, _ = mode(stacked_atlases, axis=-1, keepdims=False)
    return np.where(group_mask, group_atlas, 0)

# save atlas with mask
def save_composite_atlas(group_atlas, group_mask, output_path):
    with h5py.File(output_path, 'w') as h5file:
        h5file.create_dataset('atlas', data=group_atlas)
        h5file.create_dataset('mask', data=group_mask)
        print(f"Group composite atlas saved to {output_path}")

# and plot it out
def plot_group_composite(group_atlas, output_image_path):
    unique_labels = np.unique(group_atlas)
    cmap = plt.get_cmap('nipy_spectral', len(unique_labels))
    label_to_index = {label: i for i, label in enumerate(unique_labels)}
    indexed_atlas = np.vectorize(label_to_index.get, otypes=[int])(group_atlas)
    plt.figure(figsize=(8, 6))
    plt.imshow(indexed_atlas, cmap=cmap)
    cbar = plt.colorbar(ticks=range(len(unique_labels)))
    cbar.set_ticklabels([str(lbl) for lbl in unique_labels])
    plt.title('Group Composite Atlas (Distinct Labels)')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.grid(False)
    plt.savefig(output_image_path, dpi=300)
    plt.close()
    print(f"Group composite atlas image saved to {output_image_path}")

output_path = '/oak/stanford/groups/leanew1/users/apines/p50_mice/group_composite_atlas.h5'
output_image_path = '/oak/stanford/groups/leanew1/users/apines/p50_mice/group_composite_atlas.png'
group_atlas = create_group_composite_atlas(subjlist, data_dir, fname, group_mask, downscale, ds_dims)
save_composite_atlas(group_atlas, group_mask, output_path)
plot_group_composite(group_atlas, output_image_path)

