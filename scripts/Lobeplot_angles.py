import nibabel as nb
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scipy.io as sio
import nibabel.freesurfer.io as fsio
import sys
import os

# Arguments
subj = sys.argv[1]
sesh = sys.argv[2]
pofp = '/scratch/users/apines/' + str(subj) + '_' + str(sesh) + '_ATS.png'
# Prevent plt from trying to render mid-script
matplotlib.use('agg')

# define angular colormap: blue (0°) -> yellow (90°) -> red (180°). Same as schematic (fig 1)
angular_colors = [
    (0.290, 0.353, 0.659),  # 0° = #4a5aa8 (blue)
    (0.502, 0.765, 0.251),  # 45° = #80c340 (green)
    (0.820, 0.867, 0.153),  # 90° = #d1dd27 (yellow)
    (0.988, 0.694, 0.090),  # 135° = #fcb117 (orange)
    (0.929, 0.110, 0.141)   # 180° = #ed1c24 (red)
]

custom_cmap = LinearSegmentedColormap.from_list("angular_cmap", angular_colors, N=256)

# File paths
childfp = '/scratch/users/apines/data/mdma/' + str(subj) + '/' + str(sesh)

# Load angular timeseries
ATSl = childfp + '/' + subj + '_' + sesh + '_rs1_k1_Prop_TS_dmn_L.csv'
ATSr = childfp + '/' + subj + '_' + sesh + '_rs1_k1_Prop_TS_dmn_R.csv'
ts = np.genfromtxt(ATSl, delimiter=',')
CL = np.array(ts[:])
ts = np.genfromtxt(ATSr, delimiter=',')
CR = np.array(ts[:])

# load in DMN
DMNfpL = '/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/Group_lh_Network_1_3k.func.gii'
DMNfpR = '/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/Group_rh_Network_1_3k.func.gii'
dmn_L_vals = np.array(nb.load(DMNfpL).agg_data()).squeeze()
dmn_R_vals = np.array(nb.load(DMNfpR).agg_data()).squeeze()

# Load surface geometry for fs4
surfL = '/oak/stanford/groups/leanew1/users/apines/surf/lh.sphere'
surfR = '/oak/stanford/groups/leanew1/users/apines/surf/rh.sphere'
verts_l, faces_l = nb.freesurfer.read_geometry(surfL)
verts_r, faces_r = nb.freesurfer.read_geometry(surfR)

# load in master mask to get original 5120-propTS correspondence
master_mask_L = np.genfromtxt(os.path.expanduser('~/MasterMask_L_1.csv'), delimiter=',').astype(bool)
master_mask_R = np.genfromtxt(os.path.expanduser('~/MasterMask_R_1.csv'), delimiter=',').astype(bool)
# set num of timepoints
num_timepoints = CL.shape[1]
# Initialize 5120 x time arrays with NaN
CL_5120 = np.full((5120, num_timepoints), np.nan)
CR_5120 = np.full((5120, num_timepoints), np.nan)
CL_5120[master_mask_L,:]=CL;
CR_5120[master_mask_R,:]=CR;

# great. NOW, convert it to vertices
num_vertices_L = faces_l.max() + 1
CL_2562 = np.full((num_vertices_L, num_timepoints), np.nan)
vertex_counts_L = np.zeros(num_vertices_L)
for face_idx, face_verts in enumerate(faces_l):
    # if it has non NA values in angular time series
    if not np.isnan(CL_5120[face_idx, 0]):
	# for each vertex touching this face
        for vert_idx in face_verts:
            if np.isnan(CL_2562[vert_idx, 0]):
                CL_2562[vert_idx, :] = CL_5120[face_idx, :]
                vertex_counts_L[vert_idx] = 1
            else:
                CL_2562[vert_idx, :] += CL_5120[face_idx, :]
                vertex_counts_L[vert_idx] += 1
# now normalize by count of vertices assoc. with each face
for vert_idx in range(num_vertices_L):
    if vertex_counts_L[vert_idx] > 0:
        CL_2562[vert_idx, :] /= vertex_counts_L[vert_idx]
# Convert RIGHT hemisphere faces to vertices
num_vertices_R = faces_r.max() + 1
CR_2562 = np.full((num_vertices_R, num_timepoints), np.nan)
vertex_counts_R = np.zeros(num_vertices_R)
for face_idx, face_verts in enumerate(faces_r):
    if not np.isnan(CR_5120[face_idx, 0]):
        for vert_idx in face_verts:
            if np.isnan(CR_2562[vert_idx, 0]):
                CR_2562[vert_idx, :] = CR_5120[face_idx, :]
                vertex_counts_R[vert_idx] = 1
            else:
                CR_2562[vert_idx, :] += CR_5120[face_idx, :]
                vertex_counts_R[vert_idx] += 1

for vert_idx in range(num_vertices_R):
    if vertex_counts_R[vert_idx] > 0:
        CR_2562[vert_idx, :] /= vertex_counts_R[vert_idx]

# set to same name as vertex time series in equiv. BOLD script
CL=CL_2562.T
CR=CR_2562.T
# load in snr mask
mask_L = np.array(nb.load('/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii').agg_data())
mask_R = np.array(nb.load('/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii').agg_data())
# first dim is mask
mask_L = mask_L[0,:]==0
mask_R = mask_R[0,:]==0

# load in lobe annot
lh_labels, _, lh_names = fsio.read_annot('/home/users/apines/lh.lobes.annot')
rh_labels, _, rh_names = fsio.read_annot('/home/users/apines/rh.lobes.annot')
lh_names = [n.decode('utf-8') for n in lh_names]
rh_names = [n.decode('utf-8') for n in rh_names]
# get big lobes
# load in DMN indices
DMNfpL='/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/Group_lh_Network_1_3k.func.gii'
DMNfpR='/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/Group_rh_Network_1_3k.func.gii'
dmn_L_vals = np.array(nb.load(DMNfpL).agg_data()).squeeze()
dmn_R_vals = np.array(nb.load(DMNfpR).agg_data()).squeeze()

# subject CL and CR to snr mask
CLm = CL[:, mask_L]      # (T x V_keepL)
CRm = CR[:, mask_R]      # (T x V_keepR)

# subject lobe indices to SNR mask
lh_labels_m = lh_labels[mask_L]   # (V_keepL,)
rh_labels_m = rh_labels[mask_R]   # (V_keepR,)

# subject DMN vals to SNR mask
dmn_L_m = dmn_L_vals[mask_L]      # (V_keepL,)
dmn_R_m = dmn_R_vals[mask_R]      # (V_keepR,)

# boolean masks of lobes by hemi after SNR mask
frontal_L = np.isin(lh_labels_m, [i for i, n in enumerate(lh_names) if "frontal"  in n.lower()])
frontal_R = np.isin(rh_labels_m, [i for i, n in enumerate(rh_names) if "frontal"  in n.lower()])

parietal_L = np.isin(lh_labels_m, [i for i, n in enumerate(lh_names) if "parietal" in n.lower()])
parietal_R = np.isin(rh_labels_m, [i for i, n in enumerate(rh_names) if "parietal" in n.lower()])

temporal_L = np.isin(lh_labels_m, [i for i, n in enumerate(lh_names) if "temporal" in n.lower()])
temporal_R = np.isin(rh_labels_m, [i for i, n in enumerate(rh_names) if "temporal" in n.lower()])

# extract frontal lobe of CL/CR, DMN
CL_fr = CLm[:, frontal_L]
CR_fr = CRm[:, frontal_R]
dmnL_fr = dmn_L_m[frontal_L]
dmnR_fr = dmn_R_m[frontal_R]

# extract parietal lobe of CL/CR, DMN
CL_pa = CLm[:, parietal_L]
CR_pa = CRm[:, parietal_R]
dmnL_pa = dmn_L_m[parietal_L]
dmnR_pa = dmn_R_m[parietal_R]

# extract temporal lobe of CL/CR, DMN
CL_te = CLm[:, temporal_L]
CR_te = CRm[:, temporal_R]
dmnL_te = dmn_L_m[temporal_L]
dmnR_te = dmn_R_m[temporal_R]

# sort within lobe by DMN values and threshold at .3
thr = 0.3
# --- FRONTAL ---
# find where dmn > .3
keepL_fr = dmnL_fr >= thr
keepR_fr = dmnR_fr >= thr
# threshold frontal
CL_fr_thr = CL_fr[:, keepL_fr]
CR_fr_thr = CR_fr[:, keepR_fr]
# same to DMN for equiv sorting indices
dmnL_fr_thr = dmnL_fr[keepL_fr]
dmnR_fr_thr = dmnR_fr[keepR_fr]
# no nan vertices 
keep_notnanL_fr = ~np.all(np.isnan(CL_fr_thr), axis=0)
CL_fr_thr = CL_fr_thr[:, keep_notnanL_fr]
dmnL_fr_thr = dmnL_fr_thr[keep_notnanL_fr]
keep_notnanR_fr = ~np.all(np.isnan(CR_fr_thr), axis=0)
CR_fr_thr = CR_fr_thr[:, keep_notnanR_fr]
dmnR_fr_thr = dmnR_fr_thr[keep_notnanR_fr]

sortL_fr = np.argsort(-dmnL_fr_thr)
sortR_fr = np.argsort(-dmnR_fr_thr)

CL_fr_sorted = CL_fr_thr[:, sortL_fr].T
CR_fr_sorted = CR_fr_thr[:, sortR_fr].T


# --- PARIETAL ---
keepL_pa = dmnL_pa >= thr
keepR_pa = dmnR_pa >= thr

CL_pa_thr = CL_pa[:, keepL_pa]
CR_pa_thr = CR_pa[:, keepR_pa]
dmnL_pa_thr = dmnL_pa[keepL_pa]
dmnR_pa_thr = dmnR_pa[keepR_pa]

# no nan vertices
keep_notnanL_pa = ~np.all(np.isnan(CL_pa_thr), axis=0)
CL_pa_thr = CL_pa_thr[:, keep_notnanL_pa]
dmnL_pa_thr = dmnL_pa_thr[keep_notnanL_pa]
keep_notnanR_pa = ~np.all(np.isnan(CR_pa_thr), axis=0)
CR_pa_thr = CR_pa_thr[:, keep_notnanR_pa]
dmnR_pa_thr = dmnR_pa_thr[keep_notnanR_pa]


sortL_pa = np.argsort(-dmnL_pa_thr)
sortR_pa = np.argsort(-dmnR_pa_thr)

CL_pa_sorted = CL_pa_thr[:, sortL_pa].T
CR_pa_sorted = CR_pa_thr[:, sortR_pa].T


# --- TEMPORAL ---
keepL_te = dmnL_te >= thr
keepR_te = dmnR_te >= thr

CL_te_thr = CL_te[:, keepL_te]
CR_te_thr = CR_te[:, keepR_te]
dmnL_te_thr = dmnL_te[keepL_te]
dmnR_te_thr = dmnR_te[keepR_te]


# no nan vertices
keep_notnanL_te = ~np.all(np.isnan(CL_te_thr), axis=0)
CL_te_thr = CL_te_thr[:, keep_notnanL_te]
dmnL_te_thr = dmnL_te_thr[keep_notnanL_te]
keep_notnanR_te = ~np.all(np.isnan(CR_te_thr), axis=0)
CR_te_thr = CR_te_thr[:, keep_notnanR_te]
dmnR_te_thr = dmnR_te_thr[keep_notnanR_te]

sortL_te = np.argsort(-dmnL_te_thr)
sortR_te = np.argsort(-dmnR_te_thr)

CL_te_sorted = CL_te_thr[:, sortL_te].T
CR_te_sorted = CR_te_thr[:, sortR_te].T

# for colormap
vmin = 0
vmax = 180

lobe_data = []
lobe_labels_out = []

# Frontal
lobe_data.append(CL_fr_sorted); lobe_labels_out.append("Left Frontal")
lobe_data.append(CR_fr_sorted); lobe_labels_out.append("Right Frontal")

# Parietal
lobe_data.append(CL_pa_sorted); lobe_labels_out.append("Left Parietal")
lobe_data.append(CR_pa_sorted); lobe_labels_out.append("Right Parietal")

# Temporal
lobe_data.append(CL_te_sorted); lobe_labels_out.append("Left Temporal")
lobe_data.append(CR_te_sorted); lobe_labels_out.append("Right Temporal")

print(lobe_data[0].shape)
print(lobe_data[1].shape)
# Now plot each lobe
lobe_heights = [ts.shape[0] for ts in lobe_data]
# lobe matrices on basis of lobe heights
num_plots = len(lobe_data)
max_ratio = 900  # tweaking needed
min_ratio = 100   # min height
norm_heights = np.clip(lobe_heights, min_ratio, max_ratio)
# list out height ratios
height_ratios = list(norm_heights[:num_plots])  # truncate just in case
# set plot sizes accordingly
fig, axes = plt.subplots(
    nrows=num_plots,
    ncols=1,
    figsize=(10, 2 + sum(height_ratios)/100),
    sharex=True,
    gridspec_kw={'height_ratios': height_ratios}
)

for i, (ax, lobe_ts) in enumerate(zip(axes[0:], lobe_data)):
    ax.imshow(lobe_ts, aspect='auto', cmap=custom_cmap, vmin=vmin, vmax=vmax,
              extent=[0, lobe_ts.shape[1], 0, lobe_ts.shape[0]])
    ax.set_ylabel(lobe_labels_out[i])
    ax.set_yticks([])

axes[-1].set_xlabel('Timepoints')
plt.tight_layout()

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

cbar.set_label("Angle (degrees)")
cbar.set_ticks([0, 45, 90, 135, 180])  # optional

# save it
plt.savefig(pofp, dpi=800)
plt.close()

