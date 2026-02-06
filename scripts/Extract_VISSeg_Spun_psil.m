function Extract_DMNSeg_Spun_psil(subj,sesh,task)
% set parent directory
parentfp=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/WashU_psilocybin/' subj '/' subj '_' sesh '/func/' subj '_' sesh '_' task '.dtseries.nii'];
% define some paths 
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(Paths{1}))

% read in cifti
C=read_cifti(parentfp);
% extract time series
C_timeseries=C.cdata;

% load in temporal mask
childfp=['/scratch/users/apines/data/psil/' subj '/' sesh];
tmaskfp=[childfp '/' subj '_' sesh '_task-' task '_AllSegments.txt'];
tmask=load(tmaskfp);
% make binary mask for continuous segments
TRwise_mask_cont=zeros(1,size(C_timeseries,2));
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

models = C.diminfo{1}.models;
% Get valid surface vertex indices used in CIFTI
vl_L = models{1}.vertlist+1;        % 0-based indexing
vl_R = models{2}.vertlist+1;        % 0-based indexing
start_R = models{2}.start;        % CIFTI index (starts at 29697)

% loop over 10k spun dmns
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
outFP=['/scratch/users/apines/data/psil/' subj '/' sesh];
writetable(T,[outFP '/' subj '_' sesh '_' task '_VISSeg_Spun.csv'],'WriteRowNames',true)


