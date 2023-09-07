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
ITSl= childfp + '/' + str(subj) + '_' + str(sesh) + '_task-rs_unmasked_interp_L_faces.csv'
ITSr= childfp + '/' + str(subj) + '_' + str(sesh) + '_task-rs_unmasked_interp_R_faces.csv'
# import the left cortex time series
ts=np.genfromtxt(ITSl,delimiter=',')
CL=np.array(ts[:])
# drop ghost dimension
CL=np.squeeze(CL)
# convert to SD
Lstds=np.std(CL)
CL=CL/Lstds;
# right cortex
ts=np.genfromtxt(ITSr,delimiter=',')
CR=np.array(ts[:])
# drop ghost dimension
CR=np.squeeze(CR)
# convert to SD
Rstds=np.std(CR)
CR=CR/Rstds;
# load in DMN
dmFP_L='/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/group_L_AggNets_10k.func.gii'
dmFP_R='/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/group_R_AggNets_10k.func.gii'
dm_L=nb.load(dmFP_L).agg_data()[1]
dm_R=nb.load(dmFP_R).agg_data()[1]
### interpolate them to faces instead of vertices
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

# to organize vertices in terms of their position on the PG
DMindicesL=np.argsort(face_dmn_l)
DMindicesR=np.argsort(face_dmn_r)
# load in motion mask (temporal mask)
motMaskFp= childfp + '/' + str(subj) + '_' + str(sesh) + '_task-rs_ValidSegments_Trunc.txt'
motMask=np.genfromtxt(motMaskFp,delimiter=',')
# sort timeseries by pg
CL_sorted=CL[DMindicesL,:]
CR_sorted=CR[DMindicesR,:]
# get dmn indices for coloration
face_dmn_l_sorted=face_dmn_l[DMindicesL]
face_dmn_r_sorted=face_dmn_r[DMindicesR]
# make carpetplots
cmap = 'gray'
vmin, vmax = -3, 3
fig, ax = plt.subplots(figsize=(20, 30))
# Add vertical lines at the start of each time segment from motMask
for i, segment_start in enumerate(motMask[:, 0]): 
    # Subtract the loop iteration number from the x-value
    # -1 because python thinks 0 is 1
    # -i because opflow segments are b/w frames (only inclusive on sequence inside of head movement frames) 
    x_position = (segment_start - 1) - i
    ax.axvline(x=x_position, color='red', linestyle='--', alpha=0.3)

# denote dmn vs nondmn faces 
first_dmn_row = np.argmax(face_dmn_l_sorted > 0.3)
# make a distinct DMN vector to simplify things
ax.imshow(CL_sorted, aspect='auto', interpolation='nearest', cmap=cmap, vmin=vmin, vmax=vmax)
ax.axhline(y=first_dmn_row, color='blue', linestyle='-')
# save
outfp='/oak/stanford/groups/leanew1/users/apines/data/p50/' + subj + '/' + sesh + '/figs/CL_left.png'
plt.savefig(outfp,bbox_inches='tight')
# print out a random version for comparison (randomly organized across y axis)
# Create a random version for comparison (randomly shuffle along y-axis)
np.random.shuffle(CL_sorted)
ax.imshow(CL_sorted, aspect='auto', interpolation='nearest', cmap=cmap, vmin=vmin, vmax=vmax)
outfp='/oak/stanford/groups/leanew1/users/apines/data/p50/' + subj + '/' + sesh + '/figs/CL_left_null.png'
plt.savefig(outfp,bbox_inches='tight')

# repeat for right cortex ###

fig, ax = plt.subplots(figsize=(20, 30))
# Add vertical lines at the start of each time segment from motMask
for i, segment_start in enumerate(motMask[:, 0]):
    # Subtract the loop iteration number from the x-value
    # -1 because python thinks 0 is 1
    # -i because opflow segments are b/w frames (only inclusive on sequence inside of head movement frames) 
    x_position = (segment_start - 1) - i
    ax.axvline(x=x_position, color='red', linestyle='--', alpha=0.3)

# denote dmn vs nondmn faces 
first_dmn_row = np.argmax(face_dmn_r_sorted > 0.3)
# make a distinct DMN vector to simplify things
ax.imshow(CR_sorted, aspect='auto', interpolation='nearest', cmap=cmap, vmin=vmin, vmax=vmax)
ax.axhline(y=first_dmn_row, color='blue', linestyle='-')
# save
outfp='/oak/stanford/groups/leanew1/users/apines/data/p50/' + subj + '/' + sesh + '/figs/CR_right.png'
plt.savefig(outfp,bbox_inches='tight')
# print out a random version for comparison (randomly organized across y axis)
# Create a random version for comparison (randomly shuffle along y-axis)
np.random.shuffle(CR_sorted)
ax.imshow(CL_sorted, aspect='auto', interpolation='nearest', cmap=cmap, vmin=vmin, vmax=vmax)
outfp='/oak/stanford/groups/leanew1/users/apines/data/p50/' + subj + '/' + sesh + '/figs/CR_right_null.png'
plt.savefig(outfp,bbox_inches='tight')

