function AngDist_to_FreqDur(subj,sesh,task)

% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))

% load data
inFP=['/scratch/users/apines/data/mdma/' subj '/' sesh];
AngDistFP=[inFP '/' subj '_' sesh '_' task '_k1_AngDistMat.mat'];
AngDist=load(AngDistFP).AngDist;

% load in remaining TRs/continuous segments
parentfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
CSIfp=[parentfp '/' subj '_' sesh '_task-' task '_ValidSegments_Trunc.txt'];
CSI = importdata(CSIfp);

% assure that TR count is the same between time series and valid segments txt
SegNum=size(CSI);
SegNum=SegNum(1);
% get # of between-TR frames in valid segments: -1 because frames 1:4 have 3 OpFl estimates
BWTRFrameCounts=CSI(:,2)-1;
% number of between frames in 
mr_ts_bwFs=size(AngDist.Left,2);
mr_ts_Faces_L=size(AngDist.Left,1);
mr_ts_Faces_R=size(AngDist.Right,1);

% 1:1 mapping check of number of opfl obs and number of valid TRs
if sum(BWTRFrameCounts) ~= mr_ts_bwFs
        disp('TRs from Valid Segments txt and AngDist do not match. Fix it.')
        return
end

% initialize BUP counts and durations
BUPcount=0;
BUPdurations=[];
TDcount=0;
TDdurations=[];

% for each face LEFT
for F=1:mr_ts_Faces_L
	% for each continuous segment
	for S=1:SegNum
		% get adapted start-frame: starts at 1 for segment 1, -1*(F-1) for subsequent segments (same reason as BWTRFrameCounts)
		SF=CSI(S,1)-(S-1);
		% get indices of this segment: need (
		SegStart=SF;
		SegEnd=SegStart+BWTRFrameCounts(S)-1;
		% extract OpFl in this segment in this face
		OFMeasures=AngDist.Left(F,SegStart:SegEnd);
		%%%% derive number of BUP episodes
		% find indices where the value crosses below and above thresholds for BUP and TD
		% note it is not 90 in this iteration
		above_90 = OFMeasures > 90;
		below_90 = OFMeasures <= 90;
		% get start/end indices for each episode
		BUPstarts=[];
		BUPends=[];
		% for each timepoint (other than first, because episode definiton requires a shift from TD to BUP)
		for tp=2:(BWTRFrameCounts(S)-1)
			BUPstart=below_90(tp)==0 && below_90(tp+1)==1;
			BUPend=below_90(tp)==0 && below_90(tp-1)==1;
			if BUPstart==1
				BUPstarts=[BUPstarts (tp+1)];
			end
			if BUPend==1
				BUPends=[BUPends (tp-1)];
			end
		end
		% if any episodes found (choppy scans segements might not meet criteria
		if length(BUPends) ~=0 && length(BUPstarts) ~=0
		% if the first episode is already ongoing by the start of the segment (precedes BUP start criteria), omit it
		if BUPends(1)<BUPstarts(1)
			BUPends=BUPends(2:end);
		end
		% if there's no endpoint for the last bup episode, omit it
		if length(BUPstarts)==(length(BUPends)+1)
			BUPstarts=BUPstarts(1:end-1);
		end
		% for each BUP episode
		for B=1:length(BUPstarts)
			% +1 to count
			BUPcount=BUPcount+1;
			% calculate duration
			BUPst=BUPstarts(B);
			BUPen=BUPends(B);
			% +1 for inclusivity: episode starting at TR x and ending at the same TR is a length of 1, not 0 (and ep. from 2:3 is 2 TRs long)
			BUPdur=BUPen-BUPst+1;
			% record duration
			BUPdurations=[BUPdurations BUPdur];
		end
		% end conditional of if this segment has identifiable episodes
		end
		% for each TD episode
		% get start/end indices for each episode
                TDstarts=[];
                TDends=[];
		% for each timepoint (other than first, because episode definiton requires a shift from TD to BUP)
                for tp=2:(BWTRFrameCounts(S)-1)
                        TDstart=above_90(tp)==0 && above_90(tp+1)==1;
                        TDend=above_90(tp)==0 && above_90(tp-1)==1;
                        if TDstart==1
                                TDstarts=[TDstarts (tp+1)];
                        end
                        if TDend==1
                                TDends=[TDends (tp-1)];
                        end
                end
		% if any episodes found
		if length(TDends) ~=0 && length(TDstarts) ~=0
		 % if the first episode is already ongoing by the start of the segment (precedes TD start criteria), omit it
                if TDends(1)<TDstarts(1)
                        TDends=TDends(2:end);
                end
                % if there's no endpoint for the last bup episode, omit it
                if length(TDstarts)==(length(TDends)+1)
                        TDstarts=TDstarts(1:end-1);
                end
		% for each TD episode
                for T=1:length(TDstarts)
                        % +1 to count
                        TDcount=TDcount+1;
                        % calculate duration
                        TDst=TDstarts(T);
                        TDen=TDends(T);
                        % +1 for inclusivity: episode starting at TR x and ending at the same TR is a length of 1, not 0 (and ep. from 2:3 is 2 TRs long)
                        TDdur=TDen-TDst+1;
                        % record duration
                        TDdurations=[TDdurations TDdur];
                end
		% end conditional of if any episodes found
		end
	% end each segment
	end
% end each face
end

% for each face RIGHT
for F=1:mr_ts_Faces_R
        % for each continuous segment
        for S=1:SegNum
                % get adapted start-frame: starts at 1 for segment 1, -1*(F-1) for subsequent segments (same reason as BWTRFrameCounts)
                SF=CSI(S,1)-(S-1);
                % get indices of this segment: need (
                SegStart=SF;
                SegEnd=SegStart+BWTRFrameCounts(S)-1;
                % extract OpFl in this segment in this face
                OFMeasures=AngDist.Right(F,SegStart:SegEnd);
                %%%% derive number of BUP episodes
                % find indices where the value crosses below and above 90
                above_90 = OFMeasures > 90;
                below_90 = OFMeasures <= 90;
                % get start/end indices for each episode
                BUPstarts=[];
                BUPends=[];
                % for each timepoint (other than first, because episode definiton requires a shift from TD to BUP)
                for tp=2:(BWTRFrameCounts(S)-1)
                        BUPstart=below_90(tp)==0 && below_90(tp+1)==1;
                        BUPend=below_90(tp)==0 && below_90(tp-1)==1;
                        if BUPstart==1
                                BUPstarts=[BUPstarts (tp+1)];
                        end
                        if BUPend==1
                                BUPends=[BUPends (tp-1)];
                        end
                end
		% if any episodes found (choppy scans segements might not meet criteria
                if length(BUPends) ~=0 && length(BUPstarts) ~=0
                % if the first episode is already ongoing by the start of the segment (precedes BUP start criteria), omit it
                if BUPends(1)<BUPstarts(1)
                        BUPends=BUPends(2:end);
                end
                % if there's no endpoint for the last bup episode, omit it
                if length(BUPstarts)==(length(BUPends)+1)
                        BUPstarts=BUPstarts(1:end-1);
                end
                % for each BUP episode
                for B=1:length(BUPstarts)
                        % +1 to count
                        BUPcount=BUPcount+1;
                        % calculate duration
                        BUPst=BUPstarts(B);
                        BUPen=BUPends(B);
                        % +1 for inclusivity: episode starting at TR x and ending at the same TR is a length of 1, not 0 (and ep. from 2:3 is 2 TRs long)
                        BUPdur=BUPen-BUPst+1;
                        % record duration
                        BUPdurations=[BUPdurations BUPdur];
                end
		% end criteria of needing instances of starts and ends, as above
		end
                % for each TD episode
                % get start/end indices for each episode
                TDstarts=[];
                TDends=[];
                % for each timepoint (other than first, because episode definiton requires a shift from TD to BUP)
                for tp=2:(BWTRFrameCounts(S)-1)
                        TDstart=above_90(tp)==0 && above_90(tp+1)==1;
                        TDend=above_90(tp)==0 && above_90(tp-1)==1;
                        if TDstart==1
                                TDstarts=[TDstarts (tp+1)];
                        end
                        if TDend==1
                                TDends=[TDends (tp-1)];
                        end
                end
		% if any episodes found
                if length(TDends) ~=0 && length(TDstarts) ~=0
                 % if the first episode is already ongoing by the start of the segment (precedes TD start criteria), omit it
                if TDends(1)<TDstarts(1)
                        TDends=TDends(2:end);
                end
                % if there's no endpoint for the last bup episode, omit it
                if length(TDstarts)==(length(TDends)+1)
                        TDstarts=TDstarts(1:end-1);
                end
                % for each TD episode
                for T=1:length(TDstarts)
                        % +1 to count
                        TDcount=TDcount+1;
                        % calculate duration
                        TDst=TDstarts(T);
                        TDen=TDends(T);
                        % +1 for inclusivity: episode starting at TR x and ending at the same TR is a length of 1, not 0 (and ep. from 2:3 is 2 TRs long)
                        TDdur=TDen-TDst+1;
                        % record duration
                        TDdurations=[TDdurations TDdur];
                end
		% end conditional of if episodes exist
		end
        % end each segment
        end
% end each face
end

% Define the bins (1 to 40)
bins = 0.5:1:50.5; % Edges of bins
% histc durations BUP
BU_dur_counts = histcounts(BUPdurations, bins);
TD_dur_counts = histcounts(TDdurations, bins);
% save durations BUP
outFP=[inFP '/' subj '_' sesh '_' task '_BU_dur_counts.mat'];
save(outFP,'BU_dur_counts');
outFP=[inFP '/' subj '_' sesh '_' task '_TD_dur_counts.mat'];
save(outFP,'TD_dur_counts');
% save count BUP
outFP=[inFP '/' subj '_' sesh '_' task '_BU_counts.mat'];
save(outFP,'BUPcount');
% save count TD
outFP=[inFP '/' subj '_' sesh '_' task '_TD_counts.mat'];
save(outFP,'TDcount');
% save out mean duration
meanBUPdur=mean(BUPdurations);
outFP=[inFP '/' subj '_' sesh '_' task '_BU_meanDur.mat'];
save(outFP,'meanBUPdur');
meanTDdur=mean(TDdurations);
outFP=[inFP '/' subj '_' sesh '_' task '_TD_meanDur.mat'];
save(outFP,'meanTDdur');

