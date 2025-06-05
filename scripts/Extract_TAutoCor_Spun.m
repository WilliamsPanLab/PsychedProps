function Extract_TAutoCor_Spun(subj,sesh,task)
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
% read in good segments indicator
CSIfp=[childfp '/' subj '_' sesh '_task-' task '_AllSegments.txt'];
CSI = importdata(CSIfp);
numSegments=sum(CSI(:,3)==1);
ValidSegs=CSI(CSI(:,3)==1,:);

% need some extra cifti data to match indices from .mat
models = C.diminfo{1}.models;
% Get valid surface vertex indices used in CIFTI
vl_L = models{1}.vertlist+1;        % 0-based indexing
vl_R = models{2}.vertlist+1;        % 0-based indexing
start_R = models{2}.start;        % CIFTI index (starts at 29697)

% loop over 10k spun dmns
SpunDMNs=load('/oak/stanford/groups/leanew1/users/apines/DMNspins_32k.mat');
% initialize outputs
DMNSpunTAs=zeros(1,2000);

% for each spin quantify autocor in spun mask
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
	% mask timeseries
	DMN_timeseries=C_timeseries(DMNInds,:);
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
			dataInSegAndVert=DMN_timeseries(v,segStart:segEnd);
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
	% get mean autocor per segment
	MAC_PS=mean(ACarray,2);
	% get proportion of total TRs
	totalTRs=sum(ValidSegs(:,4));
	Proportions=ValidSegs(:,4)./totalTRs;
	% weight temporal autocors: autocor for segment i * proportion of total time series comprised of segment i
	weighted_PS=MAC_PS.*Proportions;
	% save out weighted autocor for DMN
	avTA=sum(weighted_PS);
	% add to spun values
	DMNSpunTAs(k)=avTA;
end	
% set rownames
stringVec = compose("Spin%d", 1:2000);
% saveout
T=table(DMNSpunTAs','RowNames',stringVec);
outFP=['/scratch/users/apines/data/mdma/' subj '/' sesh];
writetable(T,[outFP '/' subj '_' sesh '_' task '_TemporalAutoCor_Spun.csv'],'WriteRowNames',true)


