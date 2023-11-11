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
    # filter out non-pulse
    # set non-interest TR index
    # index out non-interest TRs
    delayMatrix=delayMatrix[:,interestTRs]
    magMatrix=magMatrix[:,interestTRs]
    # populate time series mat along TRs used
    TRsUsed=range(prevCountIndex,len(GSf)) # replace this with TRs used in indexing below for clarity, index into TRs used with gs troughs to get binary trough/notrough marker in same time series matrix
    # add in subject
    TimeSerisMat(1,TRsUsed))=subj
    # add in session
    TimeSeriesMat(2,TRsUsed)=sesh
    # add in task
    TimeSeriesMat(3,TRsUsed)=tasks[T]
    # add in global signal 
    timeSeriesMat(4,TRsUsed)=GSf
    # add in amygdala subfield signals
    timeSeriesMat(5,TRsUsed)=lAMY_r
    timeSeriesMat(6,TRsUsed)=mAMY_r
    timeSeriesMat(7,TRsUsed))=lAMY_l
    timeSeriesMat(8,TRsUsed))=mAMY_l
    # add in binary variable to mark GS troughs used
    timeSereiesMat(9,TRsUsed[GS_troughs])=1
    # populate cifti mat # later, after parsing spikes of interest via signal change
    # load functional cifti
    # extract slices of interest
    # populate ciftimat in each slice
    #### ADD DUMMY VOLUME TO TIME SERIES AND CIFTIS FOR SYNCHRONIZED VISUALIZATION
    # update counter: + 2 to get by end of TS for this segment and dummy volume
    counter=counter+2
    # update start and end index
    prevCountIndex=prevCountIndex+counter
    # update seconds deep
    secondsDeep=secondsDeep+.71
# end of loop over tasks
# save out time series mat for this subject
plt.plot(GS_troughs,GSf[GS_troughs],"x")
	figName=childfp+str(subj)+'_'+str(task) + 'GS_Troughs.png'
	plt.savefig(figName,bbox_inches='tight')
	plt.close()

fig, ax = plt.subplots(figsize=(15, 6))
ax.plot(timeSeriesMat[4,],color='black')
# conserve horizontal axis
ax2=ax.twinx()
ax2.plot(timeSeriesMat[5,],color='red')
ax2.plot(timeSeriesMat[6,],c,color='orange')
ax2.plot(timeSeriesMat[7,],c,color='red')
ax2.plot(timeSeriesMat[8,],c,color='orange')
# note they should share TR dimension and dummy volumes
# save out cifti mat for this subject
# save out time series mat for this subject
