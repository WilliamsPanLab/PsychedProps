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
pg_L=nb.load(pgFP_l)
pg_R=nb.load(pgFP_r)
# load in DMN
dmFP_L='/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/group_L_AggNets_10k.func.gii'
dmFP_R='/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/group_R_AggNets_10k.func.gii'
dm_L=nb.load(dmFP_L)
dm_R=nb.load(dmFP_R)
# interpolate them to faces instead of vertices

# to organize vertices in terms of their position on the PG
#PGindicesL=np.argsort(pgL.agg_data()[0])
#PGindicesR=np.argsort(pgR.agg_data()[0])
# load in motion mask
motMaskFp=

