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
tasks=['rs','rs2','mid','wm'];
for T in range(len(tasks)):
	task=tasks[T]
	# load in GS
	confFilepath='/oak/stanford/groups/leanew1/SHARED_DATASETS/private/psilocybin/data/derivatives/fmriprep-20.2.3/fmriprep/' + str(sname) + '/' + str(sesh) + '/func/' + str(sname) + '_' + str(sesh) + '_task-' + str(task) + '_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv'
	sTSfp=childfp + str(sname) + '_' + str(sesh) + '_' + str(task) + '_SubCortROIS.txt'
	# allow flexibility to load in rs2 filepath with slightly different filenames
	if task == 'rs2':
		confFilepath='/oak/stanford/groups/leanew1/SHARED_DATASETS/private/psilocybin/data/derivatives/fmriprep-20.2.3/fmriprep/' + str(sname) + '/' + str(sesh) + '/func/' + str(sname) + '_' + str(sesh) + '_task-rs_acq-mb_dir-pe1_run-0_desc-confounds_timeseries.tsv'
	confFile=np.genfromtxt(confFilepath,delimiter="\t")
	# drop header row from conf file
	confFile=confFile[1:,]
	# transpose so time is x-axis (columns)
	confTS=np.transpose(confFile)
	# extract GS
	GS=confTS[0,]
	# filter it
	GSf=butter_bandpass(GS,fs,lowpass,highpass,order=2)
	# normalize it for negativity checking later
	GAvg=np.mean(GSf)
	GSD=np.std(GSf)
	GSf=((GSf-GAvg)/GSD)
	# load in despiked subcortical TS
	sTS=np.genfromtxt(sTSfp)
	# extract right amyg 
	lAMY_r=sTS[18,]
	mAMY_r=sTS[19,]
	# left
	lAMY_l=sTS[43,]
	mAMY_l=sTS[44,]
	# find troughs between peaks in GS
	GS_troughs, _ = find_peaks(-GSf, distance=8)
	# calculate GS troughs with negative find_peaks
	GS_troughs = GS_troughs[GSf[GS_troughs]<0]
	# set troughs number to -1 because last is not relevant to inter trough interval
	troughsNum=len(GS_troughs)-1
	print(task)
	print(troughsNum)
	# initialize delay matrix
	delayMatrix=np.zeros((4,troughsNum))
	# and magnitude
	magMatrix=np.zeros((4,troughsNum))
	# if there are at least two troughs, we can look at delay in the peak b/w them
	if len(GS_troughs) > 1:
		# for each inter trough interval
		for t in range(troughsNum):
			tstart=GS_troughs[t]
			tend=GS_troughs[t+1]
			# use those markers to get amygdalar time series in this interval
			t_lAMY_r=lAMY_r[(tstart-1):(tend+1)]
			t_mAMY_r=mAMY_r[(tstart-1):(tend+1)]
			t_lAMY_l=lAMY_l[(tstart-1):(tend+1)]
			t_mAMY_l=mAMY_l[(tstart-1):(tend+1)]
			t_GS=GSf[(tstart-1):(tend+1)]
			# find peaks within this bin (subcort)
			peak_t_lAMY_r, _ =find_peaks(t_lAMY_r,distance=((tend+1)-(tstart-1)))
			peak_t_mAMY_r, _ =find_peaks(t_mAMY_r,distance=((tend+1)-(tstart-1)))
			peak_t_lAMY_l, _ =find_peaks(t_lAMY_l,distance=((tend+1)-(tstart-1)))
			peak_t_mAMY_l, _ =find_peaks(t_mAMY_l,distance=((tend+1)-(tstart-1)))
			# find peaks in filtered GS (GS)
			gs_peak, _ = find_peaks(t_GS,distance=((tend+1)-(tstart-1)))
			# get delay of amygdalae
			# and write catches if no peak detected
			if (len(peak_t_lAMY_r) !=0):
				delayMatrix[0,t]=gs_peak-peak_t_lAMY_r
				magMatrix[0,t]=t_lAMY_r[peak_t_lAMY_r]
			else:
				delayMatrix[0,t]=999
				magMatrix[0,t]=999
			if (len(peak_t_mAMY_r) !=0):
				delayMatrix[1,t]=gs_peak-peak_t_mAMY_r
				magMatrix[1,t]=t_mAMY_r[peak_t_mAMY_r]
			else:
				delayMatrix[1,t]=999
				magMatrix[1,t]=999
			if (len(peak_t_lAMY_l) !=0):
				delayMatrix[2,t]=gs_peak-peak_t_lAMY_l
				magMatrix[2,t]=t_lAMY_l[peak_t_lAMY_l]
			else:
				delayMatrix[2,t]=999
				magMatrix[2,t]=999
			if (len(peak_t_mAMY_l) !=0):
				delayMatrix[3,t]=gs_peak-peak_t_mAMY_l
				magMatrix[3,t]=t_mAMY_l[peak_t_mAMY_l]
			else:
				delayMatrix[3,t]=999
				magMatrix[3,t]=999
	# save out for this task and session
	# delay mat
	saveFNDM=childfp + str(subj) + '_' + str(tasks[T]) + '_delayMat.csv'
	np.savetxt(saveFNDM,delayMatrix,delimiter=",")
	# mag mat
	saveFNMM=childfp + str(subj) + '_' + str(tasks[T]) + '_MagMat.csv'
	np.savetxt(saveFNMM,magMatrix,delimiter=",")
	# gs with troughs marked plot
	plt.plot(GSf,color='black')
	plt.plot(GS_troughs,GSf[GS_troughs],"x")
	figName=childfp+str(subj)+'_'+str(task) + 'GS_Troughs.png'
	plt.savefig(figName,bbox_inches='tight')
	plt.close()
