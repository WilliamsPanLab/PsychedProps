function Extract_TAutoCor_Spun_lsd(subj,sesh,task)
% set parent directory
parentfp=['/scratch/users/apines/LSD_ICL/rest_proc/' subj '/' subj '_' sesh '_' task '_filt.dtseries.nii'];

% define some paths 
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(Paths{1}))

% read in citi
C=read_cifti(parentfp);

% extract time series
C_timeseries=C.cdata;

% load in temporal mask
childfp=['/scratch/users/apines/LSD_ICL/rest_proc/' subj];
tmaskfp=[childfp '/' subj '_' sesh '_task-' task '_AllSegments.txt'];
tmask=load(tmaskfp);
% read in good segments indicator
CSI = importdata(tmaskfp);
numSegments=sum(CSI(:,3)==1);
ValidSegs=CSI(CSI(:,3)==1,:);
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

% need some extra cifti data to match indices from .mat
models = C.diminfo{1}.models;
% Get valid surface vertex indices used in CIFTI
vl_L = models{1}.vertlist+1;        % 0-based indexing
vl_R = models{2}.vertlist+1;        % 0-based indexing
start_R = models{2}.start;        % CIFTI index (starts at 29697)

% loop over 10k spun MOT
SpunDMNs=load('/oak/stanford/groups/leanew1/users/apines/MOTspins_32k.mat');
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
        % match cifti indexing`
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
	% add to spun values
	DMNSpunTAs(k)=avTA;
end	
% set rownames
stringVec = compose("Spin%d", 1:2000);
% saveout
T=table(DMNSpunTAs','RowNames',stringVec);
outFP=['/scratch/users/apines/LSD_ICL/rest_proc/' subj];
writetable(T,[outFP '/' subj '_' sesh '_' task '_TemporalAutoCor_Spun_MOT.csv'],'WriteRowNames',true)


