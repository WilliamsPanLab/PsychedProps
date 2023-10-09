function SigTestStreams(subj,sesh,task)
% see if there are vertices more connected by streamlines than expected
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% for each task
%for task = 'rs1'
% load in subj streamlines
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
fn=[childfp '/' subj '_' sesh '_' task '_streamConnectivity_L.mat'];
subjStreams=load(fn);
% load in null streamlines
nullchildfp=['/scratch/users/apines/SimStreams'];
nullStreams=load([nullchildfp '/' subj '_' sesh '_' task '_streamConnectivity_L.mat']);
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


scatter(bdsml(:), subjStreams.AdjMatrix_L(:), 'r.');
alpha(.01)
% set equiv x and y lim
xlim([min(bdsml(:)), max(bdsml(:))]);
ylim([min(subjStreams.AdjMatrix_L(:)), max(subjStreams.AdjMatrix_L(:))]);
print(['~/' subj 'RealStreams.png'],'-dpng')
figure
scatter(bdsml(:), nullStreams.AdjMatrix_L(:), 'b.');
alpha(.01)
% set equiv x and y lim
xlim([min(bdsml(:)), max(bdsml(:))]);
ylim([min(subjStreams.AdjMatrix_L(:)), max(subjStreams.AdjMatrix_L(:))]);
print(['~/' subj 'NullStreams.png'],'-dpng')

% get indices of where streams exist
streamIndices = find(subjStreams.AdjMatrix_L > 0);
[streamInd_row, streamInd_col]=ind2sub(size(subjStreams.AdjMatrix_L), streamIndices);
% real location
real_locs = [streamInd_row, streamInd_col];
% will have to initialize countmat at full size: can limit within loop to just test where either simulated or real has streams
% use nonMW x nonMW to get "full" number

% initialize sigStreams
sig_Streams = subjStreams.AdjMatrix_L;

% for each nonero location
for L = 1:length(real_locs)
	L
	% get euclidean distance of this pair
	EucD=bdsml(real_locs(L,1),real_locs(L,2));
	% see what null range of connections is for this distance (within 1 mm)
	ThisRangeEucD=find(bdsml > (EucD-1) & bdsml < (EucD+1));
	% get 99.9th pecentile value for null connectiosn in this range
	percentile_value = prctile(nullStreams.AdjMatrix_L(ThisRangeEucD), 99.9);
	trueValue=subjStreams.AdjMatrix_L(real_locs(L,1),real_locs(L,2));
	% INSERT IF!
	if trueValue>percentile_value
		% convert to 0 if real doesnt meet this threshold
		sig_Streams(real_locs(L,1),real_locs(L,2))=0;
	else
	end
end

% save something out
save(['/scratch/users/apines/SimStreams/' subj '_' sesh '_' task  '_sigStreams.mat'],'sig_Streams','-v7.3');

% find instances 


% set diagonal to 0
%subjStreams.AdjMatrix_L_z=subjStreams.AdjMatrix_L;
%subjStreams.AdjMatrix_L_z = subjStreams.AdjMatrix_L_z - diag(diag(subjStreams.AdjMatrix_L_z));

% Calculate the similarity matrix (you may need to choose an appropriate distance metric)
similarity_matrix = pdist(subjStreams.AdjMatrix_L_z, 'euclidean');
% Perform hierarchical clustering
cluster_tree = linkage(similarity_matrix, 'average'); % 'average' linkage method
% set threshold
threshold=9000
% Create a dendrogram to visualize the clustering results
dendrogram(cluster_tree, 'ColorThreshold', threshold);
% If you want to get the cluster assignments based on the threshold
cluster_assignments = cluster(cluster_tree, 'cutoff', threshold, 'criterion', 'distance');
print('~/testClust.png','-dpng')

% add code to get vertex/count list of streams to visualize: ideally those that are connected

