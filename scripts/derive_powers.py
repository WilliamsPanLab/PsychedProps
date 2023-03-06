# for each task, load in GS, filter, load in despiked time series, find peaks in filt GS, delineate delay profile and magnitude of amygdalae relative to GS
import scipy
import nibabel as nb
import numpy as np
from scipy import stats
from scipy import signal
from scipy.signal import find_peaks
from scipy.signal import butter, filtfilt
from scipy import signal
from numpy import genfromtxt
import sys
import os.path
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

# Initialize time series as 2d matrix with 8 rows: with 1 row for each roi (+1 for global +4 amyg), 1 row for subj, 1 row for task, and 1 row for session, and undeclared temporal frames (TRs)
TimeSeriesMat=np.zeros((9,1))
# Initialize by 3 spatial frames (slices), undeclared temporal frames (TRs), along nifti dimension 1 (length dimension) and nifti dimension 2 (width dimension)
CiftiMat=np.zeros((3,1,91,109))
# initialize TR counter
counter=0
prevCountIndex=0
# initialize second counter
secondsDeep=0
### define filter parameters
# sampling freq
fs=1/.71
# lowpass
lowpass=0.10
# highpass
highpass=0.01
# Subject is set to the first passed argument
subj = sys.argv[1]
sname = subj
# sesh is set to be the second passed argument
sesh = sys.argv[2]
# row number of ROIs in time series .txt
sROIs=[18,19,43,44,0,1,2,3,25,26,27,28]
# set to match https://github.com/yetianmed/subcortex/blob/master/Group-Parcellation/3T/Subcortex-Only/Tian_Subcortex_S3_3T_label.txt (-1 bc python)
cROIs=[26,76,46,22,21,93,71]
# set to match https://github.com/ThomasYeoLab/CBIG/blob/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti/Schaefer2018_100Parcels_17Networks_order_info.txt (-1 bc python)
# set out filepath
childfp='/oak/stanford/groups/leanew1/users/apines/data/p50/' + str(sname) + '/' + str(sesh) + '/' 
# for each task (except other RS with opposite phase encoding: write in sep. loop
tasks=['rs','rs2','emotion','wm','gambling'];
for T in range(len(tasks)):
	#initialize pulse count matrix
	powers=np.zeros((51,19))
	task=tasks[T]
	# load in GS
	confFilepath='/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' + str(sname) + '/' + str(sesh) + '/func/' + str(sname) + '_' + str(sesh) + '_task-' + str(task) + '_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv'
	# subcortical time series
	sTSfp=childfp + str(sname) + '_' + str(sesh) + '_' + str(task) + '_SubCortROIS.txt'
	# cortical time series
	cTSfp=childfp + str(sname) + '_' + str(sesh) + '_' + str(task) + '_CortROIS.txt'
	# allow flexibility to load in rs2 filepath with slightly different filenames
	if task == 'rs2':
		confFilepath='/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' + str(sname) + '/' + str(sesh) + '/func/' + str(sname) + '_' + str(sesh) + '_task-rs_acq-mb_dir-pe1_run-0_desc-confounds_timeseries.tsv'
	confFile=np.genfromtxt(confFilepath,delimiter="\t")
	# drop header row from conf file
	confFile=confFile[1:,]
	# transpose so time is x-axis (columns)
	confTS=np.transpose(confFile)
	# extract FD
	FDs=confTS[16,]
	# drop na cell
	FDs=FDs[1:len(FDs)]
	# extract mean FD
	FD=np.mean(FDs)	
	# extract GS
	GS=confTS[0,]
	# load in despiked subcortical TS
	sTS=np.genfromtxt(sTSfp)
	# and cortical
	cTS=np.genfromtxt(cTSfp)
	# loop over each subcortical ROI
	for sROI in range(len(sROIs)):
		# extract ROI bold
		ROIbold=sTS[sROIs[sROI],]
		# extract PSD
		(f, S)= scipy.signal.welch(ROIbold, 1/.71,nperseg=100)
		# plop into output df
		powers[:,sROI]=S
	# now for cortical
	for cROI in range(len(cROIs)):
		# extract ROI bold
		ROIbold=cTS[cROIs[cROI],]
		# psd
		(f, S)= scipy.signal.welch(ROIbold, 1/.71,nperseg=100)
		# Into df
		powers[:,(cROI+sROI)]=S
	saveFN_pulseC=childfp + str(subj) + '_' + str(tasks[T]) + '_PSDs.csv'
	np.savetxt(saveFN_pulseC,powers,delimiter=",")
	print(task)
