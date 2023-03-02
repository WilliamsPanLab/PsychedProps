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
import sklearn
from sklearn import linear_model
import os.path
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
### define butterwoth filter
def butter_bandpass(data, fs, lowpass, highpass, order=2):
	b, a = butter(
		order / 2,
		[highpass, lowpass],
		btype="bandpass",
		output="ba",
		fs=fs
	)
	filtered_data = np.zeros_like(data)
	filtered_data = filtfilt(b,a,data,padtype="constant",padlen=data.shape[0]-1)
	return filtered_data


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
# set out filepath
childfp='/oak/stanford/groups/leanew1/users/apines/data/p50/' + str(sname) + '/' + str(sesh) + '/' 
# for each task (except other RS with opposite phase encoding: write in sep. loop
tasks=['rs','rs2','wm','gambling','emotion'];
for T in range(len(tasks)):
	task=tasks[T]
	# load in GS
	confFilepath='/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' + str(sname) + '/' + str(sesh) + '/func/' + str(sname) + '_' + str(sesh) + '_task-' + str(task) + '_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv'
	sTSfp=childfp + str(sname) + '_' + str(sesh) + '_' + str(task) + '_SubCortROIS.txt'
	# allow flexibility to load in rs2 filepath with slightly different filenames
	if task == 'rs2':
		confFilepath='/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' + str(sname) + '/' + str(sesh) + '/func/' + str(sname) + '_' + str(sesh) + '_task-rs_acq-mb_dir-pe1_run-0_desc-confounds_timeseries.tsv'
	confFile=np.genfromtxt(confFilepath,delimiter="\t")
	# drop header row from conf file
	confFile=confFile[1:,]
	# transpose so time is x-axis (columns)
	confTS=np.transpose(confFile)
	# extract GS
	GS=confTS[0,]
	# filter it
	GSf=butter_bandpass(GS,fs,lowpass,highpass,order=2)
	# load in despiked subcortical TS
	sTS=np.genfromtxt(sTSfp)
	# extract right amyg 
	lAMY_r=sTS[18,]
	mAMY_r=sTS[19,]
	# left
	lAMY_l=sTS[43,]
	mAMY_l=sTS[44,]
    # convert to percent signal change
    lAMY_r=(((GSf+lAMY_r)*100)/GSf)-100
    mAMY_r=(((GSf+mAMY_r)*100)/GSf)-100
    lAMY_l=(((GSf+lAMY_l)*100)/GSf)-100
    mAMY_l=(((GSf+mAMY_l)*100)/GSf)-100
    # .3% threshold
    thresh=.3
    # fings peaks in all ROIs
    peaks_lAMY_r=find_peaks(lAMY_r,height=thresh)
    peaks_mAMY_r=find_peaks(mAMY_r,height=thresh)
    peaks_lAMY_l=find_peaks(lAMY_l,height=thresh)
    peaks_mAMY_l=find_peaks(mAMY_l,height=thresh)
    print(task)
    print(len(peaks_lAMY_r[0]))
    print(len(peaks_mAMY_r[0]))
    print(len(peaks_lAMY_l[0]))
    print(len(peaks_mAMY_l[0]))
