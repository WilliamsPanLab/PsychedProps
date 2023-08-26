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
# import arguments
subj = sys.argv[1]
sesh = sys.argv[2]
# output fp
pofp = sys.argv[3]
#pofp = '/cbica/projects/pinesParcels/PWs/scripts/'
# prevent plt from trying to render mid script, only saveout
matplotlib.use('agg')
# set resolution
plt.rcParams['figure.dpi'] = 1000
# set parent fp
parentfp='/scratch/users/apines/data/mdma/' + str(subj) + '/' +str(sesh)
### load in lobe indices
#leftLubs=np.genfromtxt('/cbica/projects/pinesParcels/data/lobe_label_lh.csv')
#rightLubs=np.genfromtxt('/cbica/projects/pinesParcels/data/lobe_label_rh.csv')

###### import FC map
# to organize vertices in terms of their position on the PG
#PGindicesL=np.argsort(pgL.agg_data()[0])
#PGindicesR=np.argsort(pgR.agg_data()[0])
# to organize vertices in terms of their position in FC
#FCindicesL=np.argsort(L_FC)
#FCindicesR=np.argsort(R_FC)

###### import interp TS
TSfp= parentfp + '/' + str(subj) + '_' + str(sesh) + '_task-rs_p2mm_masked_interp_L.mgh'
ts=nb.load(TSfp)
CL=np.array(ts.dataobj[:])
# drop ghost dimension
CL=np.squeeze(CL)
# convert to SD
Lstds=np.std(CL)
CL=CL/Lstds;
# and transpose so time is on x-axis
CL=np.transpose(CL)
# import right hemisphere
TSfp= parentfp + '/' + str(subj) + '_' + str(sesh) + '_task-rs_p2mm_masked_interp_R.mgh'
ts=nb.load(TSfp)
CR=np.array(ts.dataobj[:])
# drop ghost dimension
CR=np.squeeze(CR)
# convert to SD
Rstds=np.std(CR)
CR=CR/Rstds;
# and transpose so time is on x-axis
CR=np.transpose(CR)

### import confounds file for GS
fmriprep_outdir='/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' +str(subj) + '/' + str(sesh) + '/func/'
# and figure dir for later
xcpd_figdir='/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' +str(subj) + '/figures/'
# Define the file paths for the confound TSV files
confFilepath1 = fmriprep_outdir + str(subj) + '_' + str(sesh) + '_task-rs_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv' 
confFilepath2 = fmriprep_outdir + str(subj) + '_' + str(sesh) + '_task-rs_acq-mb_dir-pe1_run-0_desc-confounds_timeseries.tsv' 
# Load the TSV files using numpy
conf1 = np.genfromtxt(confFilepath1, delimiter='\t', names=True, dtype=None)
conf2 = np.genfromtxt(confFilepath2, delimiter='\t', names=True, dtype=None)
# Extract the 'global_signal' columns
GS1 = conf1['global_signal']
GS2 = conf2['global_signal']
# Combine the arrays to match concatenated neuroimages
gs=np.concatenate((rs1Conf,rs2Conf),axis=0)
### import motion-masked gs
GS_file=parentfp + '/' + str(subj) + '_' + str(sesh) + '_GS_p2mm.csv';


# get the dorsomedial thalamus time series
SubcortTS1fp=parentfp + str(subj) + '_' + str(sesh) + '_rs_SubCortROIS.txt'
SubcortTS2fp=parentfp + str(subj) + '_' + str(sesh) + '_rs2_SubCortROIS.txt'
# concatenate
SubcortTS=np.concatenate((np.genfromtxt(SubcortTS1fp,delimiter='\t'),np.genfromtxt(SubcortTS2fp,delimiter='\t')),axis=1)
# take 34 and 9th row for left and right DMThal
DMThalL=SubcortTS[33,]
DMThalR=SubcortTS[8,]

# for time series (plot) and PSD (psd)
import matplotlib.pyplot as plt
# plt.psd(s, 512, 1 / diff)

# each vertex normalized w/r/t global SD
# R
Rstds=np.std(CR)
CR=CR/Rstds;
# sort cortical data /w/r/t FC
CL=CL[:,FCindicesL]
CR=CR[:,FCindicesR]
# transpose for grayplots
CLt=np.transpose(CL)
CRt=np.transpose(CR)

# plot power spectral density of percentile bins of FC-ranked cortical time series

###### plot

# Define the layout of subplots using plt.subplots()
fig, axs = plt.subplots(6, 1, figsize=(30, 10),gridspec_kw={'height_ratios': [1, 1, 1, 5, 1, 5]})

# Plot gsdt on the first subplot
axs[0].plot(gsdt)
axs[i].autoscale(enable=True, axis='both', tight=True)

# Plot gs on the second subplot
axs[1].plot(gs)
axs[i].autoscale(enable=True, axis='both', tight=True)

# Plot DMThalL on the third subplot
axs[2].plot(DMThalL)
axs[i].autoscale(enable=True, axis='both', tight=True)

# Plot the matrix of imshow on the fourth subplot
cmap = 'gray'
vmin, vmax = -3, 3
im = axs[3].imshow(CLt, aspect='auto', interpolation='nearest',cmap=cmap,vmin=vmin,vmax=vmax)
a=fig.colorbar(im)
a.remove()

# Plot DMThalR on the fifth subplot
axs[4].plot(DMThalR)
axs[i].autoscale(enable=True, axis='both', tight=True)

# Plot CRt on the sixth subplot
im = axs[5].imshow(CRt, aspect='auto', interpolation='nearest',cmap=cmap,vmin=vmin,vmax=vmax)
fig.colorbar(im, ax=axs[5])
a=fig.colorbar(im)
a.remove()

# Set the titles for each subplot
axs[0].set_title('gsdt')
axs[1].set_title('gs')
axs[2].set_title('DMThalL')
axs[3].set_title('CLt')
axs[4].set_title('DMThalR')
axs[5].set_title('CRt')

# Set the x-axis label for the matrix of imshow subplot
axs[3].set_xlabel('Time')

# Adjust the spacing between subplots to avoid overlap
plt.subplots_adjust(hspace=0.5)

# Save the figure
plt.savefig('my_figure.png', dpi=300)


# calculate data volume:100 aspect ratio correction number
ARCN=(CL.shape[0])/100
# Left Hemi
plt.subplot(3,2,1)
plt.matshow(CLt,cmap='gray',vmin=-2, vmax=2,extent=[0,100,0,1], aspect=100,fignum=False)
# Right Hemi
plt.subplot(3,2,4)
plt.matshow(CRt,cmap='gray',vmin=-2, vmax=2,extent=[0,100,0,1], aspect=100,fignum=False)
###### save plot to subj folder
plt.savefig(pofp + str(subj) + '_DMThal_tall.png',bbox_inches='tight')

