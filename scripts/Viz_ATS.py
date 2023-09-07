import nibabel as nb
#import hcp_utils as hcp
import numpy as np
import nilearn.plotting as plotting
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import copy
import sys
import os
from PIL import Image
# import arguments
subj = sys.argv[1]
sesh = sys.argv[2]
# output fp
ofp='/oak/stanford/groups/leanew1/users/apines/data/p50/' + str(subj) + '/' + str(sesh) + '/figs'
# prevent plt from trying to render mid-script, only upon saveout
matplotlib.use('agg')
# child filepath
childfp='/scratch/users/apines/data/mdma/' + str(subj) + '/' + str(sesh)
# timeseries filepath
ATSl=childfp + '/' + subj + '_' + sesh + '_Prop_TS_dmn_L.csv'
ATSr=childfp + '/' + subj + '_' + sesh + '_Prop_TS_dmn_R.csv'
# import the left cortex time series
ts=np.genfromtxt(ATSl,delimiter=',')
CL=np.array(ts[:])
# right cortex
ts=np.genfromtxt(ATSr,delimiter=',')
CR=np.array(ts[:])
# load in DMN
import scipy.io as sio

# Load Networks from .mat files
networks_data = sio.loadmat('/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/group_Nets_fs5.mat')
nets_LH = networks_data['nets']['Lnets'][0, 0]
nets_RH = networks_data['nets']['Rnets'][0, 0]

# Initialize matrix for each face over each of k=4 networks to save out to scratch
faceMatrix = np.zeros((len(g_noMW_combined_L) + len(g_noMW_combined_R), 4))

# Network of interest
k = 1  # 1th (2nd) network is dmn
dm_L = nets_LH[:, k]
dm_R = nets_RH[:, k]
# convert to face form
# read in surface geometry
subjects_folder = '/oak/stanford/groups/leanew1/users/apines/fs5surf'
surf_l = subjects_folder + '/lh.sphere'
surf_r = subjects_folder + '/rh.sphere'
# Load the surface data
verts_l,faces_l = nb.freesurfer.read_geometry(surf_l)
verts_r,faces_r = nb.freesurfer.read_geometry(surf_r)

# Initialize output face matrices for PG and DMN
face_dmn_l = np.zeros((faces_l.shape[0]))
face_dmn_r = np.zeros((faces_r.shape[0]))
# Interpolate DMN data onto faces for the left hemisphere
for i, face in enumerate(faces_l):
    vertex_indices = face
    face_dmn_l[i] = np.mean(dm_L[vertex_indices])

# Interpolate DMN data onto faces for the right hemisphere
for i, face in enumerate(faces_r):
    vertex_indices = face
    face_dmn_r[i] = np.mean(dm_R[vertex_indices])

# get gmedialwall mask from extract_relativeAngles script
medialwall_data = sio.loadmat('/oak/stanford/groups/leanew1/users/apines/fs5surf/medial_wall_vectors.mat')
g_noMW_combined_L=medialwall_data['g_noMW_combined_L']
g_noMW_combined_R=medialwall_data['g_noMW_combined_R']
# python adjust: these are indices not scalars
g_noMW_combined_L=g_noMW_combined_L[0]-1
g_noMW_combined_R=g_noMW_combined_R[0]-1
# Use medial wall mask as a common starting point (from which to mask both opfl vecs and net grads further)
ng_L = face_dmn_l[g_noMW_combined_L]
ng_R = face_dmn_r[g_noMW_combined_R]

# load in nagrad vertices
n_a_gradverts = sio.loadmat('/oak/stanford/groups/leanew1/users/apines/fs5surf/medial_wall_nullGrad_vectors.mat')
InclLeft=n_a_gradverts['InclLeft']
InclRight=n_a_gradverts['InclRight']
# python adjust
InclLeft=InclLeft-1
InclRight=InclRight-1
face_dmn_l=ng_L[InclLeft]
face_dmn_r=ng_R[InclRight]
# to organize vertices in terms of their position on the PG
DMindicesL = np.argsort(np.ravel(face_dmn_l))
DMindicesR = np.argsort(np.ravel(face_dmn_r))
# back to parallelism with vis_ITS.py script

# load in motion mask (temporal mask)
motMaskFp= childfp + '/' + str(subj) + '_' + str(sesh) + '_task-rs_ValidSegments_Trunc.txt'
motMask=np.genfromtxt(motMaskFp,delimiter=',')
# sort timeseries by pg
CL_sorted = CL[DMindicesL, :]
CR_sorted = CR[DMindicesR, :]
# get dmn indices for coloration
face_dmn_l_sorted=face_dmn_l[DMindicesL]
face_dmn_r_sorted=face_dmn_r[DMindicesR]
# make carpetplots
cmap = 'jet'
vmin, vmax = 0, 180
fig, ax = plt.subplots(figsize=(20, 30))
# Add vertical lines at the start of each time segment from motMask
for i, segment_start in enumerate(motMask[:, 0]): 
    # Subtract the loop iteration number from the x-value
    # -1 because python thinks 0 is 1
    # -i because opflow segments are b/w frames (only inclusive on sequence inside of head movement frames) 
    x_position = (segment_start - 1) - i
    ax.axvline(x=x_position, color='black', linestyle='--', alpha=0.5)

# denote dmn vs nondmn faces 
first_dmn_row = np.argmax(face_dmn_l_sorted > 0.3)
# make a distinct DMN vector to simplify things
ax.imshow(CL_sorted, aspect='auto', interpolation='nearest', cmap=cmap, vmin=vmin, vmax=vmax,alpha=.7)
ax.axhline(y=first_dmn_row, color='blue', linestyle='-')
# save
outfp='/oak/stanford/groups/leanew1/users/apines/data/p50/' + subj + '/' + sesh + '/figs/ACL_left.png'
plt.savefig(outfp,bbox_inches='tight')
# print out a random version for comparison (randomly organized across y axis)
# Create a random version for comparison (randomly shuffle along y-axis)
np.random.shuffle(CL_sorted)
ax.imshow(CL_sorted, aspect='auto', interpolation='nearest', cmap=cmap, vmin=vmin, vmax=vmax,alpha=.7)
outfp='/oak/stanford/groups/leanew1/users/apines/data/p50/' + subj + '/' + sesh + '/figs/ACL_left_null.png'
plt.savefig(outfp,bbox_inches='tight')
# now create the same for the same for the angular time series! 
