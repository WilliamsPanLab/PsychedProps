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
correctFactor=846/numTRsVS;
subjStreams.AdjMatrix_L=subjStreams.AdjMatrix_L .* correctFactor;

% get indices of where streams exist
streamIndices = find(subjStreams.AdjMatrix_L > 0);
[streamInd_row, streamInd_col]=ind2sub(size(subjStreams.AdjMatrix_L), streamIndices);
% real location
real_locs = [streamInd_row, streamInd_col];
% will have to initialize countmat at full size: can limit within loop to just test where either simulated or real has streams
% use nonMW x nonMW to get "full" number
% initialize count matrix of equal size, which will record # of times simulated connectivity exceeds real
countMat=zeros(length(nonMW_L),length(nonMW_L));

% set simulation parent fp
simfp=['/scratch/users/apines/'];

% for each of 100 simulated streamlines, see if more fake streamlines exist between vertices than real data
for i=1:100
	i
	tic
	% simulated matrix name
	sfn=[simfp 'SimStreams/' num2str(i) '_streamConnectivity_L_sparse.mat'];
	if exist(sfn,'file');
		% load in simulated matrix
		simmat=load(sfn);
		[simstream_row, simstream_col]=find(simmat.sparseSimMat);
		sim_locs = [simstream_row, simstream_col];
		% get union of where real and simulated streams exist to test 
		Union_locs=union(real_locs,sim_locs,'rows');
		% check simulated vs. observed in each stream index
		for v=1:length(Union_locs)
			% extract locations
			coords=Union_locs(v,:);
			% observed
			Obs=subjStreams.AdjMatrix_L(coords(1),coords(2));
			% simulated
			[~,~,Sim]=find(simmat.sparseSimMat(coords(1),coords(2)));
			% record if observed streams exceeds simulated
			if Obs>Sim
				countMat(coords(1),coords(2))=countMat(coords(1),coords(2))+1;
			elseif Sim>Obs
				countMat(coords(1),coords(2))=countMat(coords(1),coords(2))-1;
			end
		end
		% end each sim loop
		toc
	% end conditional check for file
	else
	end
end

% find where real streams exceeded all simulated
exceedingStreams=countMat > 99;
% find where real streams fell below all simulated
belowStreams=countMat < -99;
% sparse versions
exceedingStreams_s=sparse(exceedingStreams);
belowStreams_s=sparse(belowStreams);
% save out processed info
save([childfp '/exceedingStreams.mat'],'exceedingStreams_s')
save([childfp '/belowStreams.mat'],'belowStreams_s')
% print out some strong examples
% Sort rows in descending order of sum for exceedingStreams
[b, sorted_exceeding_indices] = sort(sum(exceedingStreams, 2), 'descend');
% Select the top 5 rows
top_exceeding_rows = sorted_exceeding_indices(1:5);

% Sort rows in descending order of sum for belowStreams
[b, sorted_below_indices] = sort(sum(belowStreams, 2), 'ascend');
% Select the top 5 rows
top_below_rows = sorted_below_indices(1:5);

% Iterate over the top rows for exceedingStreams
for i = 1:numel(top_exceeding_rows)
    current_row = top_exceeding_rows(i);
    % Process the current row
    Vis_VertvecFs5(log(subjStreams.AdjMatrix_L(current_row, :)), zeros(1, 10242), [childfp '/AboveLocus_' num2str(i) '.png']);
end

% Iterate over the bottom rows for belowStreams (lowest # of streams arrived)
for i = 1:numel(top_below_rows)
    current_row = top_below_rows(i);
    % Process the current row
    Vis_VertvecFs5(log(subjStreams.AdjMatrix_L(current_row, :)), zeros(1, 10242), [childfp '/BelowLocus_' num2str(i) '.png']);
end
