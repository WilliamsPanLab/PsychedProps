function Extract_DMNSeg_Spun(subj,sesh,task)
% set parent directory
parentfp=['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' subj '/' sesh '/func'];

% define some paths 
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(Paths{1}))

% re-adjust for rsfmri naming conventions
if string(task)=="rs1"
	fp=[parentfp '/' subj '_' sesh '_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii'];
	C=read_cifti(fp);
elseif string(task)=="rs2"
	fp=[parentfp '/' subj '_' sesh '_task-rs_acq-mb_dir-pe1_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii'];
	C=read_cifti(fp);
else
% read in time series
C=read_cifti([parentfp '/' subj '_' sesh '_task-' task '_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii']);
end

% extract time series
C_timeseries=C.cdata;

% load in temporal mask
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh];
tmaskfp=[childfp '/' subj '_' sesh '_task-' task '_AllSegments.txt'];
tmask=load(tmaskfp);
% just get ts length from fmriprep outputs to be safe/consistent
confFilepath1=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subj '/' sesh '/func/' subj '_' sesh '_task-' task '_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv'];
% adapt if it is resting state
if string(task)=="rs1"
	confFilepath1=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subj '/' sesh '/func/' subj '_' sesh '_task-rs_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv'];
elseif string(task)=="rs2"
	confFilepath1=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subj '/' sesh '/func/' subj '_' sesh '_task-rs_acq-mb_dir-pe1_run-0_desc-confounds_timeseries.tsv'];
else
end
conf1=readtable(confFilepath1,"FileType","text",'Delimiter', '\t');
% extract FD columns
FD=table2array(conf1(:,'framewise_displacement'));

% make binary mask for continuous segments
TRwise_mask_cont=zeros(1,length(FD));
% Loop through each row in Absolut
for row = 1:size(tmask, 1)
	if tmask(row, 3) == 1
        	% Extract the start and end values from the current row
                startValue = tmask(row, 1);
                endValue = tmask(row, 2);
                % change TRwise_mask_cont to 1 where this sequence of continuous good TRs occurs
                TRwise_mask_cont(startValue:endValue)=1;
        else
        end
end

% get correlation matrix of full time series... just cortex
C_timeseries=C_timeseries(1:59412,logical(TRwise_mask_cont));
% Get valid surface vertex indices used in CIFTI
models = C.diminfo{1}.models;
vl_L = models{1}.vertlist+1;        % 0-based indexing
vl_R = models{2}.vertlist+1;        % 0-based indexing
start_R = models{2}.start;        % CIFTI index (starts at 29697)

% loop over 10k spun VIS
SpunDMNs=load('/oak/stanford/groups/leanew1/users/apines/VISspins_32k.mat');
% initialize outputs
DMNSpunSegs=zeros(1,2000);
for k=1:2000
	% get spun DMN
	SpunDMN_L=SpunDMNs.bigrotl_32k(k,:);
	SpunDMN_R=SpunDMNs.bigrotr_32k(k,:);
	% omit spun mw
	SpunDMN_L(SpunDMN_L>1)=0;
	SpunDMN_R(SpunDMN_R>1)=0;
	% omit nans
	SpunDMN_L(isnan(SpunDMN_L))=0;
	SpunDMN_R(isnan(SpunDMN_R))=0;
	% match cifti indexing
	SpunDMN_L=SpunDMN_L(vl_L);
	SpunDMN_R=SpunDMN_R(vl_R);
	%combined DMN
	combinedMask= false(59412, 1);
	% based of cifti indices
	combinedMask(1:29696) = SpunDMN_L;
	combinedMask(29697:59412) = SpunDMN_R;
	DMNInds=find(combinedMask>.3);
	% get non-dmn indices
	nonDMNinds=setdiff(1:59412,DMNInds);
	DMN_mat=C_timeseries(DMNInds,:);
	nonDMN_mat=C_timeseries(nonDMNinds,:);
	% avoid making full correlation matrix due to memory demands
	correlation_matrix = 1 - pdist2(DMN_mat, nonDMN_mat, 'correlation');  
	BWFC=mean(mean(correlation_matrix,'omitnan'),'omitnan');
	DMNSpunSegs(k)=BWFC;
end
% create rownames
stringVec = compose("Spin%d", 1:2000);
% save out
T=table(DMNSpunSegs','RowNames',stringVec);
outFP=['/scratch/users/apines/data/mdma/' subj '/' sesh];
writetable(T,[outFP '/' subj '_' sesh '_' task '_VISSeg_Spun.csv'],'WriteRowNames',true)


