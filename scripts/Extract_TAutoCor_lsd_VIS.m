function Extract_TAutoCor_lsd(subj,sesh,task)
% set parent directory
parentfp=['/scratch/users/apines/LSD_ICL/rest_proc/' subj '/' subj '_' sesh '_' task '_filt.dtseries.nii'];

% define some paths 
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(Paths{1}))

% load in VIS
DMN=ft_read_cifti_mod('/oak/stanford/groups/leanew1/users/apines/maps/Network2_fslr.dscalar.nii');
DMNInds=find(DMN.data>.3);

% read in cifti
C=ft_read_cifti_mod(parentfp);

% extract time series
C_timeseries=C.data;

% read in good segments indicator
childfp=['/scratch/users/apines/LSD_ICL/rest_proc/' subj];
CSIfp=[childfp '/' subj '_' sesh '_task-' task '_AllSegments.txt'];
CSI = importdata(CSIfp);
numSegments=sum(CSI(:,3)==1);
ValidSegs=CSI(CSI(:,3)==1,:);

% mask timeseries
C_timeseries=C_timeseries(DMNInds,:);

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
% one scan from subj 14 and all from subj three is missing a tiny amount of vertices which is throwing this for a loop: omit NA rows just for them
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
outFP=['/scratch/users/apines/LSD_ICL/rest_proc/' subj];
writetable(T,[outFP '/' subj '_' sesh '_' task '_TemporalAutoCor_VIS.csv'],'WriteRowNames',true)


