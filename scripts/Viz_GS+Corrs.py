import nibabel as nb
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
# sesh = sys.argv[2]
# output fp
#pofp = '/cbica/projects/pinesParcels/PWs/scripts/'
# prevent plt from trying to render mid script, only saveout
matplotlib.use('agg')
# set resolution
plt.rcParams['figure.dpi'] = 1000
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
# initialize GS array for each session
gsArr=np.zeros((4,846))
# initialize PG x GSxBOLD correlation array
PGxGS=np.zeros(4)
# for each session: sesh = ses-00, ses-01, ses-02, ses-03
for ses in range(0,4):
    # convert numerical value to ses-00 format
    if ses == 0:
        sesh = 'ses-00'
    elif ses == 1:
        sesh = 'ses-01'
    elif ses == 2:
        sesh = 'ses-02'
    elif ses == 3:
        sesh = 'ses-03'
    # set parent fp
    parentfp='/scratch/users/apines/data/mdma/' + str(subj) + '/' +str(sesh)
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
    # Concatenate the global signal arrays
    gs = np.concatenate((GS1,GS2),axis=0)
    # load in 
    gsArr[ses,:]=gs;
    # load in niftis to calculate regional correlation with GS
    niftiPath1 = fmriprep_outdir + str(subj) + '_' + str(sesh) + '_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_bold.dtseries.nii'
    niftiPath2 = fmriprep_outdir + str(subj) + '_' + str(sesh) + '_task-rs_acq-mb_dir-pe1_run-0_space-fsLR_den-91k_bold.dtseries.nii'
    # Load the nifti files using nibabel
    nifti1 = nb.load(niftiPath1)
    nifti2 = nb.load(niftiPath2)
    # Extract the data from the nifti files
    nifti1Data = nifti1.get_fdata()
    nifti2Data = nifti2.get_fdata()
    # Combine the arrays to match concatenated neuroimages
    niftiData = np.concatenate((nifti1Data,nifti2Data),axis=0)
    # Calculate the correlation between the global signal and the neuroimage
    # Create an empty array to store the correlation coefficients
    correlations = np.zeros(niftiData.shape[1])
    for i in range(niftiData.shape[1]):
        correlation = np.corrcoef(gs, niftiData[:, i])[0, 1]
        correlations[i] = correlation
    ### make a cifti of the correlation values at each grayordinate
    # create a cifti image object from example dscalar
    exampleCiftifp='/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients.dscalar.nii'
    exampleCifti = nb.load(exampleCiftifp)
    # annoying cifti code
    axes = [exampleCifti.header.get_axis(i) for i in range(exampleCifti.ndim)]
    time_axis, brain_model_axis = axes
    # extract cifti data
    ciftiData=exampleCifti.dataobj[:,:]
    # get pg
    PG=ciftiData[0,:]
    # get correlation of PG with GS correlations for this session
    corrPG=np.corrcoef(PG,correlations)[0,1]
    PGxGS[ses]=corrPG
    # replace PG with correlation values
    ciftiData[0,:]=correlations;
    # replace first map with correlation values
    new_cifti = nb.Cifti2Image(ciftiData,header=(time_axis, brain_model_axis),nifti_header=exampleCifti.nifti_header)
    # save out the cifti
    ciftiOutfp=parentfp + '/gsCorr_' + str(sesh) + '.dscalar.nii'
    nb.save(new_cifti, ciftiOutfp)

###### plot the four gsArr time series on the same plot (GS), insert value of PGxGS for each session in-plot

# Plot the GS time series
plt.figure(figsize=(12, 6))
for ses in range(4):
    plt.plot(gsArr[ses], label=f'GS {ses}')
    plt.text(0, gsArr[ses][0], 'PGxGS: {:.2f}'.format(PGxGS[ses]), fontsize=12, color='black')

plt.xlabel('Time')
plt.ylabel('Global Signal (GS)')
plt.legend()
plt.title('Global Signal (GS) Time Series with PGxGS Values')
plt.grid()
pltname='/home/users/apines/' + str(subj) + '_gs.png'
plt.savefig(pltname,bbox_inches='tight')

