function MotMask(subj,sesh)
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
fpL=[childfp '/' subj '_' sesh '_L_AggTS_10k.mgh'];
fpR=[childfp '/' subj '_' sesh '_R_AggTS_10k.mgh'];

% if file exists, run it
if isfile(fpL)
	% read in mgh
	fpL
	dataL=MRIread(fpL);
	dataR=MRIread(fpR);
	% squeeze to get rid of extra dimensions
	tsl=squeeze(dataL.vol);
	tsr=squeeze(dataR.vol);
	% get size
	tsSize=size(tsl);
	numTRs=tsSize(2);
	% load in mask
	confFilepath1=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subj '/' sesh '/func/' subj '_' sesh '_task-rs_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv'];
	conf1=readtable(confFilepath1,"FileType","text",'Delimiter', '\t');
	confFilepath2=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subj '/' sesh '/func/' subj '_' sesh '_task-rs_acq-mb_dir-pe1_run-0_desc-confounds_timeseries.tsv'];
	conf2=readtable(confFilepath2,"FileType","text",'Delimiter', '\t');
	% extract FD columns
	FD1=table2array(conf1(:,'framewise_displacement'));
	FD2=table2array(conf2(:,'framewise_displacement'));
	% extract GS columns
	GS1=table2array(conf1(:,'global_signal'));
	GS2=table2array(conf2(:,'global_signal'));
	% combine them to match concatenated neuroimages
	FD=vertcat(FD1,FD2);
	GS=vertcat(GS1,GS2);	
	% get to FD_thresh of .2 mm
	TRwise_mask=FD>.2;
	% building in sanity check. Motion removed frames + interrupted sequence removed frames + ending frames should = start frame #
	motRemFrames=sum(TRwise_mask)
	% length of mask corresponds to number of TRs
	% 1 indicates flagged for FD over selected threshold, reverse 'em so 0 indicates removal
	TRwise_mask=~TRwise_mask;
	% use it to mask GS in parallel with time series
	GSm=GS(TRwise_mask);
	% remove TRs with corresp. flag
	masked_trs_l=tsl(:,TRwise_mask);
	masked_trs_r=tsr(:,TRwise_mask);
	% reconfig cifti metadata to reflect new number of TRs
	newSize=size(masked_trs_l);
	newTRnum=newSize(2);
	% setting continuous frame threshold to 8 TRs in a row
	Threshold=8			
	% find changepoints in binary bask
	d = [true, diff(TRwise_mask') ~= 0];
	% index of changepoints
	dInd=find(d);
	% find difference in indices of changepoints (span of mask/non-mask epochs)
	n = diff([dInd, numTRs]); 
	% find which segments correspond to non-mask
	maskValAtChange=TRwise_mask(dInd);
	ContSegments=n(:,maskValAtChange);
	% create list of starting TR and duration of segments uninterupt. by combined mask
	UTSegSize=size(ContSegments);
	UTSegNum=UTSegSize(2);
	UTSegCell=cell(UTSegNum,2);
	% plant in duration of clean segments
	for i=1:UTSegNum
		UTSegCell(i,2)=num2cell(ContSegments(i));
	end
	% make 1st column start position in .2mm outlier masked sequence
	% (just the start where prev. segment left off, no masked TRs in gaps b/w)
	UTSegCell(1,1)=num2cell(1);
	for i=2:UTSegNum
		UTSegCell(i,1)=num2cell(UTSegCell{i-1,1}+UTSegCell{i-1,2});
	end
	% now we need to artificially insert the discontinuity between scans. That is, if one segment spans the end of the first scan and into the second, it needs to be broken up to reflect the true discontinuity
      	% find any segment that runs over end of scan 1 into scan 2
       	CrossSegment = find(cell2mat(UTSegCell(:,1)) < length(FD1) & cell2mat(UTSegCell(:,1)) + cell2mat(UTSegCell(:,2)) > length(FD1));
	% Iterate over the extended segments and adjust the duration
	for i = 1:length(CrossSegment)
    		% Calculate the new duration for the segment
    		newSegmentDuration = (length(FD1) - UTSegCell{CrossSegment, 1})+1;
    		% record where the new segment will starts
		newSegmentStart = length(FD1)+1;
		% record where the new segment will end
		remainingSegmentDuration = (UTSegCell{CrossSegment, 2} - newSegmentDuration)+1;
		% Update the duration in UTSegCell
    		UTSegCell{CrossSegment, 2} = newSegmentDuration;
    		% Insert a new row for the remaining portion of the extended segment
    		UTSegCell = [UTSegCell(1:CrossSegment, :); {newSegmentStart, remainingSegmentDuration}; UTSegCell(CrossSegment+1:end, :)];
	end	
	% update ContSegments
	ContSegments=cell2mat(UTSegCell(:,2));
	% find segments with more continuous TRs than threshold
        OverThreshSegments=find(ContSegments>Threshold);
        % sanity check for TRs excluded for being in interrupted segments
	UnderThreshSegments=find(ContSegments<(Threshold+1));
	ExcludedTRs=sum(ContSegments(UnderThreshSegments));
	% sum remaining segments to get included TRs if this thresh chosen
        RemainingTRs=sum(ContSegments(OverThreshSegments));
	% index of which TR valid segments start at
        ValidTRStarts=dInd(maskValAtChange);
	% index out segments greater than TR thresh from UnThreshSegmentCellstruct
	ValidSegCell=UTSegCell(OverThreshSegments,:);
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
	% save 2-column df indicating start of valid segments and length: matches dtseries at derive_wave
	segmentfnTr=[childfp '/' subj '_' sesh '_task-rs_ValidSegments_Trunc'];
	writetable(cell2table(ValidSegCell_Trunc),segmentfnTr,'WriteVariableNames',0)
	% make binary mask for continuous segments
	TRwise_mask_cont=zeros(1,newTRnum);
	% need to start double checking here through assertion
	% comprised of 1 at starting TR through length of valid segment
	for seg=1:length(OverThreshSegments);
		TRwise_mask_cont(ValidSegCell{seg,1}:(ValidSegCell{seg,1}+ValidSegCell{seg,2}-1))=1;
	end
	% apply to GS
	GSmc=GSm(logical(TRwise_mask_cont));
	% apply to mgh
	masked_trs_cont_l=masked_trs_l(:,logical(TRwise_mask_cont));
	masked_trs_cont_r=masked_trs_r(:,logical(TRwise_mask_cont));
	% insert back into OG structure to saveout in equivalent structure
	dataL.vol=masked_trs_cont_l;
	dataR.vol=masked_trs_cont_r;
	% set output filepath
	ofpl=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_L.mgh'];
	ofpr=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_R.mgh'];
        % implement sanity check. FD > thresh removed + noncontinuous segment removed + remaining = OG
	assert(motRemFrames + ExcludedTRs + sum(TRwise_mask_cont) == numTRs, 'Assertion failed: Over FD trs + noncontinuous TRs + remaining TRs is not equal to original TRs');
	% saveout
	MRIwrite(dataL,ofpl)
	MRIwrite(dataR,ofpr)
end
