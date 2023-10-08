function SigTestStreams(subj,sesh)
% see if there are vertices more connected by streamlines than expected
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))

% for each task
for task = 'rs1'
% load in subj streamlines
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
fn=[childfp '/' subj '_' sesh '_' task '_streamConnectivity_L.mat'];
subjStreams=load(fn);
% load in null streamlines
nullchildfp=['/scratch/users/apines/SimStreams'];
nullSteams=load([nullchildfp 'Simulated_OpFl_' subj '_' sesh '_' task '.mat']);
% load in medial wall to mask stream connectivity
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
% use native freesurfer command for mw mask indices
surfML = [SubjectsFolder '/lh.Medial_wall.label'];
mwIndVec_l = read_medial_wall_label(surfML);
% make "isn't medial wall" indices
nonMW_L=setdiff(1:2562,mwIndVec_l);
% extract non-medial wall
subjStreams.AdjMatrix_L=subjStreams.AdjMatrix_L(nonMW_L,nonMW_L);
nullStreams.AdjMatrix_L=nullStreams.AdjMatrix_L(nonMW_L,nonMW_L);
% load in euclidean distances
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
% for surface data
surfL = [SubjectsFolder '/lh.pial'];
% surface topography
surfL = read_surf(surfL);
% initialize distance matrix
bdsml=zeros(length(surfL),length(surfL));
% for all vertices
for V=1:length(bdsml);
		% vertex props (x,y,z coords)
		initVertL=surfL(V,:);
		xVL=initVertL(1);
		yVL=initVertL(2);
		zVL=initVertL(3);
        
		% search through all of them for eucl. dist. calc.
		for i=1:length(bdsml);
			xL=surfL(i,1);
			yL=surfL(i,2);
			zL=surfL(i,3);
			eucld_L=sqrt((xL-xVL)^2+(yL-yVL)^2+(zL-zVL)^2);
			bdsml(V,i)=eucld_L;
		end	
end
% same non-medial wall extraction
bdsml=bdsml(nonMW_L,nonMW_L);



scatter(bdsml(:), nullStreams.AdjMatrix_L(:), 'b.');


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
