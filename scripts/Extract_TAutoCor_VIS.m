function Extract_TAutoCor(subj,sesh,task)
% set parent directory
parentfp=['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' subj '/' sesh '/func'];

% define some paths 
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(Paths{1}))

% re-adjust for rsfmri naming conventions
if string(task)=="rs1"
	fp=[parentfp '/' subj '_' sesh '_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii'];
	C=ft_read_cifti_mod(fp);
elseif string(task)=="rs2"
	fp=[parentfp '/' subj '_' sesh '_task-rs_acq-mb_dir-pe1_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii'];
	C=ft_read_cifti_mod(fp);
else
% read in time series
C=ft_read_cifti_mod([parentfp '/' subj '_' sesh '_task-' task '_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii']);
end

% extract time series
C_timeseries=C.data;

% load in DMN
DMN=ft_read_cifti_mod('/oak/stanford/groups/leanew1/users/apines/maps/Network2_fslr.dscalar.nii');
DMNInds=find(DMN.data>.3);

% mask timeseries
C_timeseries=C_timeseries(DMNInds,:);

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

% read in good segments indicator
CSIfp=[childfp '/' subj '_' sesh '_task-' task '_AllSegments.txt'];
CSI = importdata(CSIfp);
numSegments=sum(CSI(:,3)==1);
ValidSegs=CSI(CSI(:,3)==1,:);

% initialize an autocor array
ACarray=zeros(numSegments,length(DMNInds));

% loop over each segement
for n=1:numSegments
	% pull out TR numbers for this segment
	segStart=ValidSegs(n,1);
	segEnd=ValidSegs(n,2);
	% loop over each vertex within DMN
	for v=1:length(DMNInds);
		% extract data
		dataInSegAndVert=C_timeseries(v,segStart:segEnd);
		% OG time series (-1)
		OGTS=dataInSegAndVert(1:end-1);
		ShiftTS=dataInSegAndVert(2:end);
		% calculate autocor (t-1 vs. t)
		correlationOfInterest=corrcoef(OGTS,ShiftTS);
		% note this returns a correlation matrix: only interested in the between-timeseries measurement (off-diag)
		ACarray(n,v)=correlationOfInterest(1,2);
	% end vertex loop
	end	
	% end segment loop
end
% omit nans induced from very slight extra medial wall exclusion
cols_with_nan = any(isnan(ACarray),1);
col_indices_with_nan = find(cols_with_nan);
ACarray=ACarray(:,~cols_with_nan);
% get mean autocor per segment
MAC_PS=mean(ACarray,2);

% get proportion of total TRs
totalTRs=sum(ValidSegs(:,4));
Proportions=ValidSegs(:,4)./totalTRs;
% weight temporal autocors: autocor for segment i * proportion of total time series comprised of segment i
weighted_PS=MAC_PS.*Proportions;
% save out weighted autocor for DMN
avTA=sum(weighted_PS);

T=table(avTA,'RowNames',"Row1");
outFP=['/scratch/users/apines/data/mdma/' subj '/' sesh];
writetable(T,[outFP '/' subj '_' sesh '_' task '_TemporalAutoCor_VIS.csv'],'WriteRowNames',true)


