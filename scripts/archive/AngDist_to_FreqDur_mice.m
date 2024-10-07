function AngDist_to_FreqDur_mice(subj,sesh)
% load data
inFP=['/scratch/users/apines/data/mouse/'];
AngDistFP=[inFP subj '_' num2str(sesh) '_Prop_TS_dmn.csv']
AngDist=readmatrix(AngDistFP);

% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))

% initialize BUP counts and durations
BUPcount=0;
BUPdurations=[];
TDcount=0;
TDdurations=[];

% for each pixel
for F=1:size(AngDist,1)
	% extract OpFl in this segment in this face
	OFMeasures=AngDist(F,:);
	%%%% derive number of BUP episodes
	% find indices where the value crosses below and above thresholds for BUP and TD
	% note it is not 90 in this iteration
	above_90 = OFMeasures > 90;
	below_90 = OFMeasures <= 90;
	% get start/end indices for each episode
	BUPstarts=[];
	BUPends=[];
	% for each timepoint (other than first, because episode definiton requires a shift from TD to BUP)
	for tp=2:size(AngDist,1)
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
        for tp=2:size(AngDist,1)
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
% end each face
end

% conversion for saveout
sesh=num2str(sesh);

% Define the bins (1 to 40)
bins = 0.5:1:50.5; % Edges of bins
% histc durations BUP
BU_dur_counts = histcounts(BUPdurations, bins);
TD_dur_counts = histcounts(TDdurations, bins);
% save durations BUP
outFP=[inFP '/' subj '_' sesh '_BU_dur_counts.mat'];
save(outFP,'BU_dur_counts');
outFP=[inFP '/' subj '_' sesh '_TD_dur_counts.mat'];
save(outFP,'TD_dur_counts');
% save count BUP
outFP=[inFP '/' subj '_' sesh '_BU_counts.mat'];
save(outFP,'BUPcount');
% save count TD
outFP=[inFP '/' subj '_' sesh '_TD_counts.mat'];
save(outFP,'TDcount');
% save out mean duration
meanBUPdur=mean(BUPdurations);
outFP=[inFP '/' subj '_' sesh '_BU_meanDur.mat'];
save(outFP,'meanBUPdur');
meanTDdur=mean(TDdurations);
outFP=[inFP '/' subj '_' sesh '_TD_meanDur.mat'];
save(outFP,'meanTDdur');

