function MotMask(subj,sesh,task)
% this function motion masks the concatenated runs 
% note that both validsegcell_full and validsegcell_trunc are saved out: only _Trunc used in this workflow
% validsegcell_full refers to the TR and span of segment in the outlier FD masked images
% validsegcell truncated refers to the TR and span of segment in the outlierFD masked AND short segment masked images
% minproc and FD leave this script matching truncated, proc TS and GS leave this script matching full
% (full for BandPass, Trunc for OpFl)

% add paths for cifti stuff
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% filepath for output
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh];

% mask concatenated data
fpL=[childfp '/' subj '_' sesh '_L_' task '_TS_3k.mgh'];
fpR=[childfp '/' subj '_' sesh '_R_' task '_TS_3k.mgh'];

% setting continuous frame threshold to 8 TRs in a row
Threshold=8			

% if file exists, run it
if isfile(fpL)
	% read in mgh
	dataL=MRIread(fpL);
	dataR=MRIread(fpR);
	% squeeze to get rid of extra dimensions
	tsl=squeeze(dataL.vol);
	tsr=squeeze(dataR.vol);
	% get size
	tsSize=size(tsl);
	numTRs=tsSize(2);
	% load in mask
	confFilepath1=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/connectome/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subj '/' sesh '/func/' subj '_' sesh '_task-' task '_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv'];
	% adapt if it is resting state
	if string(task)=="rs1"
		confFilepath1=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/connectome/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subj '/' sesh '/func/' subj '_' sesh '_task-rs_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv'];
	elseif string(task)=="rs2"
		confFilepath1=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/connectome/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subj '/' sesh '/func/' subj '_' sesh '_task-rs_acq-mb_dir-pe1_run-0_desc-confounds_timeseries.tsv'];
	else
	end
	conf1=readtable(confFilepath1,"FileType","text",'Delimiter', '\t');
	% extract FD columns
	FD=table2array(conf1(:,'framewise_displacement'));
	% extract GS columns
	GS=table2array(conf1(:,'global_signal'));
	% combine them to match concatenated neuroimages
	% get to FD_thresh of .2 mm
	TRwise_mask=FD>.2;
	% building in sanity check. Motion removed frames + interrupted sequence removed frames + ending frames should = start frame #
	motRemovFrames=sum(TRwise_mask)
	% 1 indicates flagged for FD over selected threshold, reverse 'em so 0 indicates removal
	TRwise_mask=~TRwise_mask;
	% find changepoints in binary bask
	d = [true, diff(TRwise_mask') ~= 0];
	% index of changepoints
	dInd=find(d);
	% make segment lists in terms of absolute TR (startingTR:endingTR:include?) and relative TR (starting TR out of continuous segments: length of segment)
	% -1 because we are looking at segments between the bookends, but +1 again becasue d doesn't return the final bookend
	Absolut=cell(length(dInd),4);
	% populate 1st column in absolute with starting point of each segment
	Absolut(:,1)=num2cell(dInd(1:(length(dInd))));
	% populate 2nd column with ending point of each segment
	Absolut(:,2)=[num2cell(dInd(2:(length(dInd)))) numTRs];
	% and subtract 1 from ending point because the ending point is the impacted volume
	Absolut(:,2)=num2cell(cell2mat(Absolut(:,2))-1);
	% (except for final endpoint, which does end at very last volume)
	Absolut(end,2)=num2cell(Absolut{end,2}+1);
	% find difference in indices of changepoints (span of mask/non-mask epochs)
	n = diff([dInd, numTRs]); 
	% find which segments correspond to non-mask
	maskValAtChange=TRwise_mask(dInd);
	% make last column in Absolut boolean include based on motion
	Absolut(:,3)=num2cell(maskValAtChange);
	% to do so, first add a duration column to absolut
	Absolut(:,4)=cellfun(@minus, Absolut(:,2), Absolut(:,1), 'UniformOutput', false);
       	% and plus 1 because the TR range is inclusive: i.e., frames 1-3 is actually 3 frames, not 3 - 1 frames
	Absolut(:,4)=num2cell(cell2mat(Absolut(:,4))+1) 
	% extract continuous segments as those where boolean include is = 1 in absolut
	continuousSegments = Absolut([Absolut{:, 3}] == 1, 4);
	% create a relative list of starting TR and duration of segments uninterupt. by combined mask
	UTSegSize=size(continuousSegments);
	UTSegNum=UTSegSize(1);
	UTSegCell=cell(UTSegNum,2);
	% plant in duration of clean segments
	for i=1:UTSegNum
		UTSegCell(i,2)=continuousSegments(i);
	end
	% make 1st column start position in .2mm outlier masked sequence
	% (just the start where prev. segment left off, no masked TRs in gaps b/w)
	UTSegCell(1,1)=num2cell(1);
	for i=2:UTSegNum
		UTSegCell(i,1)=num2cell(UTSegCell{i-1,1}+UTSegCell{i-1,2});
	end
	% find segments with more continuous TRs than threshold
        OverThreshSegments=find(cell2mat(continuousSegments)>=Threshold);
        % sanity check for TRs excluded for being in interrupted segments
	UnderThreshSegments=find(cell2mat(continuousSegments)<(Threshold));
	ExcludedTRs=sum(cell2mat(continuousSegments(UnderThreshSegments)));
	% sum remaining segments to get included TRs if this thresh chosen
        RemainingTRs=sum(cell2mat(continuousSegments(OverThreshSegments)))
	% index of which TR valid segments start at
	TRStarts = [Absolut{[Absolut{:, 3}] == 1,1}];	
	ValidTRStarts=TRStarts(OverThreshSegments);
	% index out segments greater than TR thresh from UnThreshSegmentCellstruct
	ValidSegCell=UTSegCell(OverThreshSegments,:);
	% adjust Absolut boolean to reflect minimum continous threshold inclusion criterion
	Absolut(cell2mat(Absolut(:, 4)) < Threshold, 3) = num2cell(0);
	% adjust ValidSegCell_Trunc to describe in terms of retained TRs only
	% in other words, epoch start points should be one "TR" after 
	ValidSegCell_Trunc=ValidSegCell;
	ValidSegCell_Trunc(1,1)=num2cell(1);
	for s=2:length(OverThreshSegments)
		prevStart=ValidSegCell_Trunc(s-1,1);
		prevLength=ValidSegCell_Trunc(s-1,2);
		NewStart=prevStart{:}+prevLength{:};
		ValidSegCell_Trunc(s,1)=num2cell(NewStart);
	end	
	% save 2-column df indicating start of valid segments and length: matches masked, relative dtseries
	segmentfnTr=[childfp '/' subj '_' sesh '_task-' task '_ValidSegments_Trunc'];
	writetable(cell2table(ValidSegCell_Trunc),segmentfnTr,'WriteVariableNames',0)
	% write 4 column, absolute segment output
	asegmentfnTr=[childfp '/' subj '_' sesh '_task-' task '_AllSegments.txt'];
	writetable(cell2table(Absolut),asegmentfnTr,'WriteVariableNames',0)
	% make binary mask for continuous segments
	TRwise_mask_cont=zeros(1,length(FD));
	% Loop through each row in Absolut
	for row = 1:size(Absolut, 1)
	    if Absolut{row, 3} == 1
	        % Extract the start and end values from the current row
	        startValue = Absolut{row, 1};
	        endValue = Absolut{row, 2};        
    		% change TRwise_mask_cont to 1 where this sequence of continuous good TRs occurs
		TRwise_mask_cont(startValue:endValue)=1;
	    else
	    end
	end
	% apply to GS
	GSmc=GS(logical(TRwise_mask_cont));
	% apply to mgh
	masked_trs_cont_l=tsl(:,logical(TRwise_mask_cont));
	masked_trs_cont_r=tsr(:,logical(TRwise_mask_cont));
	% insert back into OG structure to saveout in equivalent structure
	dataL.vol=masked_trs_cont_l;
	dataR.vol=masked_trs_cont_r;
	% set output filepath
	ofpl=[childfp '/' subj '_' sesh '_task-' task '_p2mm_masked_L.mgh'];
	ofpr=[childfp '/' subj '_' sesh '_task-' task '_p2mm_masked_R.mgh'];
        % implement sanity check. FD > thresh removed + noncontinuous segment removed + remaining = OG
	assert(motRemovFrames + ExcludedTRs + sum(TRwise_mask_cont) == numTRs, 'Assertion failed: Over FD trs + noncontinuous TRs + remaining TRs is not equal to original TRs');
	% writeout global signal
	GSP=[childfp '/' subj '_' sesh '_' task '_GS_p2mm.csv'];
	csvwrite(GSP,GSmc)
	% saveout time series
	MRIwrite(dataL,ofpl)
	MRIwrite(dataR,ofpr)
end

