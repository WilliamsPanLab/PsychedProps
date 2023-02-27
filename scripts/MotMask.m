function apply_motion_mask(subj,sesh)
% initialize empty vector for average length
TRvecNum=[];
sname=subj;
% filepath for output
ProjectFolder = '/scratch/groups/leanew1/xcpd_outMDMA_36p_nodespike/xcp_d/sub-MDMA001/ses-00/func';
ResultantFolder = ['/cbica/projects/hcpd/data/motMasked_contSegs/'];

% note two resting states with .71 TR
tasks=['rs','rs2','wm','mid','gambling','emotion'];

% for each task
for t=1:length(tasks)
	% for this task
	task=tasks[t];
	% WRITE CATCH FOR RS2: DIFFERENT FILEPATHS
	if
	else
	% retrieve confounds file
	confFilepath=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' sname '/' sesh '/func/' sname '_' sesh '_task-' task '_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv'];
	% retrieve data
	fpParent=['/scratch/groups/leanew1/xcpd_outMDMA_36p_nodespike/xcp_d/' sname '/' sesh '/files/MNINo/'];
	fp=[fpParent 'task-' task '_DCANBOLDProc_v4.0.0_Atlas.dtseries.nii'];
