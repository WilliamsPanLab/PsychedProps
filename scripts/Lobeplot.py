import nibabel as nb
import numpy as np
import nilearn.plotting as plotting
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import copy
import sys
import os
import nibabel.freesurfer.io as fsio

# define workbench colormap
from matplotlib.colors import LinearSegmentedColormap

# Define base colors in 0–255
base_colors = np.array([
    [255, 255, 0],
    [255, 200, 0],
    [255, 120, 0],
    [255, 0, 0],
    [200, 0, 0],
    [150, 0, 0],
    [100, 0, 0],
    [60, 0, 0],
    [0, 0, 80],
    [0, 0, 170],
    [75, 0, 125],
    [125, 0, 160],
    [75, 125, 0],
    [0, 200, 0],
    [0, 255, 0],
    [0, 255, 255]
]) / 255.0  # scale to 0–1
# Define interpolation positions
interpsteps = np.linspace(0, 1, 16)
# Interpolate to 255 colors
interpolated = np.array([
    np.interp(np.linspace(0, 1, 255), interpsteps, base_colors[:, i])
    for i in range(3)
]).T  # shape: (255, 3)
# Flip to put yellow as high
interpolated = interpolated[::-1]
# Reduce overly bright ends: slice like MATLAB roybigbl_cm(15:240,:)
reduced = interpolated[14:240]
# Create custom colormap
custom_cmap = LinearSegmentedColormap.from_list("roybigbl_cm", reduced)
# load in lobes
lh_labels, _, lh_names = fsio.read_annot('/home/users/apines/lh.lobes.annot')
rh_labels, _, rh_names = fsio.read_annot('/home/users/apines/rh.lobes.annot')
lh_names = [n.decode('utf-8') for n in lh_names]
rh_names = [n.decode('utf-8') for n in rh_names]
# get big lobes
def get_lobe_mask(labels, names, keyword):
    return np.isin(labels, [i for i, name in enumerate(names) if keyword in name])
lh_frontal = get_lobe_mask(lh_labels, lh_names, 'frontal')
rh_frontal = get_lobe_mask(rh_labels, rh_names, 'frontal')
lh_parietal = get_lobe_mask(lh_labels, lh_names, 'parietal')
rh_parietal = get_lobe_mask(rh_labels, rh_names, 'parietal')
lh_temporal = get_lobe_mask(lh_labels, lh_names, 'temporal')
rh_temporal = get_lobe_mask(rh_labels, rh_names, 'temporal')
# import subj and sesh
subj = sys.argv[1]
sesh = sys.argv[2]
# output fp
pofp = '/scratch/users/apines/' + str(subj) + '_' + str(sesh) + '.png'
# prevent plot from trying to render midscript
matplotlib.use('agg')
# set resolution
plt.rcParams['figure.dpi'] = 500
# load in data
TSfp='/scratch/users/apines/data/mdma/' + str(subj) + '/' + str(sesh) + '/' + str(subj) + '_' + str(sesh) + '_L_rs1_TS_3k.func.gii'
gif_CL=nb.load(TSfp)
CL=np.array(gif_CL.agg_data()[:])
TSfp='/scratch/users/apines/data/mdma/' + str(subj) + '/' + str(sesh) + '/' + str(subj) + '_' + str(sesh) + '_R_rs1_TS_3k.func.gii'
gif_CR=nb.load(TSfp)
CR=np.array(gif_CR.agg_data()[:])

# load in snr mask
mask_L = np.array(nb.load('/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii').agg_data())
mask_R = np.array(nb.load('/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii').agg_data())
# first dim is mask
mask_L = mask_L[0,:]==0
mask_R = mask_R[0,:]==0

# load in head motion TS
xcpd_outdir='/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' +str(subj) + '/' + str(sesh) + '/func/'
rs1ConfTsv= xcpd_outdir + subj + '_' + sesh + '_task-rs_acq-mb_dir-pe0_run-0_motion.tsv'
rs1Conf=np.genfromtxt(rs1ConfTsv,delimiter='\t')
FD=rs1Conf[1:,6]
# load in DMN indices
DMNfpL='/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/Group_lh_Network_1_3k.func.gii'
DMNfpR='/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/Group_rh_Network_1_3k.func.gii'
dmn_L_vals = np.array(nb.load(DMNfpL).agg_data()).squeeze()
dmn_R_vals = np.array(nb.load(DMNfpR).agg_data()).squeeze()

# SNR mask it
CL = CL[:,mask_L]
CR = CR[:,mask_R]
dmn_L_vals = dmn_L_vals[mask_L]
dmn_R_vals = dmn_R_vals[mask_R]

# sort time series by DMN indices
sort_idx_L = np.argsort(dmn_L_vals)
sort_idx_R = np.argsort(dmn_R_vals)
CL_sorted = CL[:,sort_idx_L]
CR_sorted = CR[:,sort_idx_R]
dmn_L_sorted = dmn_L_vals[sort_idx_L]
dmn_R_sorted = dmn_R_vals[sort_idx_R]

# find DMN cutoff
cut_L = np.argmax(dmn_L_sorted > 0.3)
cut_R = np.argmax(dmn_R_sorted > 0.3)

# and transpose so CL and CR are over time
CL_sorted=CL_sorted.T
CR_sorted=CR_sorted.T

# for colormap
all_data = np.concatenate([CL_sorted.flatten(), CR_sorted.flatten()])
vmin = np.percentile(all_data, 1)
vmax = np.percentile(all_data, 99)

# plot time series
lobes_of_interest = ['frontal', 'parietal', 'temporal']
lobe_data = []
lobe_labels_out = []
# get the info to plot each lobe
for lobe in lobes_of_interest:
    for hemi, labels, mask, dmn_vals, ts_data in zip(
        ['Left', 'Right'],
        [lh_labels, rh_labels],
        [mask_L, mask_R],
        [dmn_L_vals, dmn_R_vals],
        [CL, CR]
    ):
        # get vertex indices for this lobe
        lobe_idx = [i for i, name in enumerate(lh_names if hemi == 'Left' else rh_names)
                    if lobe in name.lower()]
        lobe_vert = np.isin(labels, lobe_idx)
        lobe_masked = lobe_vert[mask]  # after SNR masking
        if np.sum(lobe_masked) == 0:
            continue
        ts_lobe = ts_data[:, lobe_masked]
        dmn_lobe = dmn_vals[lobe_masked]
        # Sort by DMN (dmn at "top")
        sort_idx = np.argsort(-dmn_lobe)
        ts_sorted = ts_lobe[:, sort_idx].T  # transpose for plotting
        lobe_data.append(ts_sorted)
        lobe_labels_out.append(f"{hemi} {lobe.capitalize()}")

# Now plot FD + each lobe
lobe_heights = [ts.shape[0] for ts in lobe_data]
# lobe matrices on basis of lobe heights
num_plots = len(lobe_data)
max_ratio = 900  # tweaking needed
min_ratio = 300   # min height
norm_heights = np.clip(lobe_heights, min_ratio, max_ratio)
# add 150 for FD for visibility
height_ratios = [150] + list(norm_heights[:num_plots])  # truncate just in case
# set plot sizes accordingly
fig, axes = plt.subplots(
    nrows=num_plots + 1,
    ncols=1,
    figsize=(10, 2 + sum(height_ratios)/330),
    sharex=True,
    gridspec_kw={'height_ratios': height_ratios}
)

axes[0].plot(FD, color='red', linewidth=0.4)
axes[0].axhline(0.2, color='black', linestyle='--', linewidth=0.5)
axes[0].set_ylabel('FD')
axes[0].set_xlim([0, len(FD)])

for i, (ax, lobe_ts) in enumerate(zip(axes[1:], lobe_data)):
    ax.imshow(lobe_ts, aspect='auto', cmap=custom_cmap, vmin=vmin, vmax=vmax,
              extent=[0, lobe_ts.shape[1], 0, lobe_ts.shape[0]])
    ax.set_ylabel(lobe_labels_out[i])
    ax.set_yticks([])

axes[-1].set_xlabel('Timepoints')
plt.tight_layout()

# add motion mask
motion_mask_fp = f'/scratch/users/apines/data/mdma/{subj}/{sesh}/{subj}_{sesh}_task-rs1_AllSegments.txt'
segments = np.genfromtxt(motion_mask_fp, delimiter=',', dtype=int)

# Determine valid intervals (start, stop where column 2 == 1)
valid_intervals = [(row[0]-1, row[1]-1) for row in segments if row[2] == 1]
all_timepoints = np.arange(FD.shape[0])

# Build a mask: 1 where valid, 0 where invalid
valid_mask = np.zeros_like(FD, dtype=bool)
for start, end in valid_intervals:
    valid_mask[start:end+1] = True

# Shade invalid areas across all plots
for ax in axes:
    # Get height of axis in data units (use max for imshow axes)
    ybot, ytop = ax.get_ylim()
    for t in range(len(valid_mask)):
        if not valid_mask[t]:
            ax.axvspan(t, t+1, color='gray', alpha=0.65, linewidth=0)

# Optional: add text or highlight on FD axis
axes[0].set_title('Temporal Mask of Representative Scan')

norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=custom_cmap)
sm.set_array([])

# add horizontal colorbar at bottom of figure
cbar = fig.colorbar(
    sm,
    ax=axes,
    orientation="horizontal",
    fraction=0.05,   # thickness
    pad=0.08,        # spacing from plots
    aspect=40        # length relative to thickness
)

cbar.set_label("BOLD Signal")

# save it
plt.savefig(pofp, dpi=500)
plt.close()
