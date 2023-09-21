function SigTestStreams(subj,sesh)
% see if there are vertices more connected by streamlines than expected
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% load in subj streamlines
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
fn=[childfp '/' subj '_' sesh '_streamConnectivity_L.mat'];
subjStreams=load(fn);

% load in medial wall to mask stream connectivity
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/fs5surf';
% use native freesurfer command for mw mask indices
surfML = [SubjectsFolder '/lh.Medial_wall.label'];
mwIndVec_l = read_medial_wall_label(surfML);
% make "isn't medial wall" indices
nonMW_L=setdiff(1:10242,mwIndVec_l);
subjStreams.AdjMatrix_L=subjStreams.AdjMatrix_L(nonMW_L,nonMW_L);

% load in number of TRs
CSIfp=[childfp '/' subj '_' sesh '_task-rs_ValidSegments_Trunc.txt'];
CSI = importdata(CSIfp);
% trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
numTRsVS=sum(CSI(:,2));

% correct by # of TRs: relative to 846 for simulated
correctFactor=1+(1-(numTRsVS/846))
subjStreams.AdjMatrix_L=subjStreams.AdjMatrix_L .* correctFactor;

% get indices of where streams exist
streamIndices = find(subjStreams.AdjMatrix_L > 0);
[streamInd_row, streamInd_col]=ind2sub(size(subjStreams.AdjMatrix_L), streamIndices);

% initialize count matrix of equal size, which will record # of times simulated connectivity exceeds real
countMat=zeros?:??
% set simulation parent fp
simfp=['/scratch/users/apines/'];
% for each of 100 simulated streamlines, see if more fake streamlines exist between vertices than real data
for i=1:100
	% simulated matrix name
	sfn=[simfp 'SimStreams/' num2str(seed) '_streamConnectivity_L_sparse.mat'];
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
