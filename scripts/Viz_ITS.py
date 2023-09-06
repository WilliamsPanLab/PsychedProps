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
ITSl= childfp + '/' + str(subj) + '_' + str(sesh) + '_task-rs_p2mm_masked_interp_L_faces.csv'
ITSr= childfp + '/' + str(subj) + '_' + str(sesh) + '_task-rs_p2mm_masked_interp_R_faces.csv'
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
# load PG
pgFP_L='/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients_10k_L.dscalar.func.gii'
pgFP_R='/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients_10k_R.dscalar.func.gii' 
pg_L=nb.load(pgFP_L).agg_data()[0]
pg_R=nb.load(pgFP_R).agg_data()[0]
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
face_pg_l = np.zeros((faces_l.shape[0]))
face_pg_r = np.zeros((faces_r.shape[0]))
face_dmn_l = np.zeros((faces_l.shape[0]))
face_dmn_r = np.zeros((faces_r.shape[0]))
# Interpolate PG data onto faces for the left hemisphere
for i, face in enumerate(faces_l):
    vertex_indices = faces_l[face,]
    face_pg_l[i] = np.mean(pg_L[vertex_indices])

# Interpolate PG data onto faces for the right hemisphere
for i, face in enumerate(faces_r):
    vertex_indices = face
    face_pg_r[i] = np.mean(pg_R[vertex_indices])

# Interpolate DMN data onto faces for the left hemisphere
for i, face in enumerate(faces_l):
    vertex_indices = face
    face_dmn_l[i] = np.mean(dm_L[vertex_indices])

# Interpolate DMN data onto faces for the right hemisphere
for i, face in enumerate(faces_r):
    vertex_indices = face
    face_dmn_r[i] = np.mean(dm_R[vertex_indices])

# to organize vertices in terms of their position on the PG
PGindicesL=np.argsort(face_pg_l)
PGindicesR=np.argsort(face_pg_r)
# load in motion mask
motMaskFp= childfp + '/' + str(subj) + '_' + str(sesh) + '_task-rs_ValidSegments_Trunc.txt'
motMask=np.genfromtxt(motMaskFp,delimiter=',')
# sort timeseries by pg
# make carpetplots
# add vertical lines at motion mask segments
# set opacity < 1 for non-dmn regions
