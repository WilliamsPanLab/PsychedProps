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
parentfp='/oak/stanford/groups/leanew1/users/apines/data/mdma/' + str(subj) + '/' +str(sesh) + '/'
### load in lobe indices
#leftLubs=np.genfromtxt('/cbica/projects/pinesParcels/data/lobe_label_lh.csv')
#rightLubs=np.genfromtxt('/cbica/projects/pinesParcels/data/lobe_label_rh.csv')

###### import FC map
#PGfp= parentfp2 + '_gradients_L_3k.func.gii'
#pgL=nb.load(PGfp)
#PGfp= parentfp2 + '_gradients_R_3k.func.gii'
#pgR=nb.load(PGfp)
L_FCfp=parentfp + str(subj) + '_' + str(sesh) + '_L_ROIfc.csv'
R_FCfp=parentfp + str(subj) + '_' + str(sesh) + '_R_ROIfc.csv'
L_FC=np.genfromtxt(L_FCfp,delimiter=',')
R_FC=np.genfromtxt(R_FCfp,delimiter=',')

# to organize vertices in terms of their position on the PG
#PGindicesL=np.argsort(pgL.agg_data()[0])
#PGindicesR=np.argsort(pgR.agg_data()[0])
# to organize vertices in terms of their position in FC
FCindicesL=np.argsort(L_FC)
FCindicesR=np.argsort(R_FC)

###### import TS
TSfp='/scratch/users/apines/' + str(subj) + '_' + str(sesh) + '_L_AggTS_3k.func.gii'
gif_CL=nb.load(TSfp)
CL=np.array(gif_CL.agg_data()[:])
TSfp='/scratch/users/apines/' + str(subj) + '_' + str(sesh) + '_R_AggTS_3k.func.gii'
gif_CR=nb.load(TSfp)
CR=np.array(gif_CR.agg_data()[:])

### import OpFl amplitude

### import confounds file for GS
xcpd_outdir='/scratch/groups/leanew1/xcpd_outMDMA_36p_despike_bp/xcp_d/' +str(subj) + '/' + str(sesh) + '/func/'

rs1ConfTsv= xcpd_outdir + subj + '_' + sesh + '_task-rs_acq-mb_dir-pe0_run-0_design.tsv'
rs2ConfTsv= xcpd_outdir + subj + '_' + sesh + '_task-rs_acq-mb_dir-pe1_run-0_design.tsv'
rs1Conf=np.genfromtxt(rs1ConfTsv,delimiter='\t')
rs2Conf=np.genfromtxt(rs2ConfTsv,delimiter='\t')
# concatenate the two time series
rsConf=np.concatenate((rs1Conf,rs2Conf),axis=0)
# get the global signal
# column 27 is gs, 30 is dt_gs
gs=rsConf[:,26]
gsdt=rsConf[:,29]
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
# L
Lstds=np.std(CL)
CL=CL/Lstds;
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

