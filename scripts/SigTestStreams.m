function SigTestStreams(subj,sesh)
% see if there are vertices more connected by streamlines than expected

% load in subj streamlines
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
fn=[childfp '/' subj '_' sesh '_streamConnectivity_L.mat'];
subjStreams=load(fn);

% load in number of TRs
CSIfp=[childfp '/' subj '_' sesh '_task-rs_ValidSegments_Trunc.txt'];
CSI = importdata(CSIfp);
% trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
numTRsVS=CSI(SegNum,1)+CSI(SegNum,2)-1;

% correct by # of TRs: relative to 846 for simulated
correctFactor=
subjStreams.x=subjStreams.x./numTRsVS

% get indices of where streams exist
streamIndices = find(subjStreams.x > 0);

% initialize count matrix of equal size, which will record # of times simulated connectivity exceeds real

% set simulation parent fp
simfp=['/scratch/users/apines/'];
% for each of 100 simulated streamlines, see if more fake streamlines exist between vertices than real data
for i=1:100
	% simulated matrix name
	sfn=[simfp 'SimStreams/' seed '_streamConnectivity_L.mat'];
	% load in simulated matrix
	simmat=load(sfn);
	% check simulated vs. observed in each stream index
	for v=streamIndices
		% record if observed streams exceeds simulated
		
	% end each sim loop
end
% consider visvertvec of row with most remaining streams

% consider printing out adjacency matrix with same row highlighted

% consider combining vertex-wise values across rows connected by common sig. threads?
