% to combine simulated streamline matrices and compress them with sparse + single
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% load in medial wall
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/fs5surf';
% use native freesurfer command for mw mask indices
surfML = [SubjectsFolder '/lh.Medial_wall.label'];
mwIndVec_l = read_medial_wall_label(surfML);
% make "isn't medial wall" indices
nonMW_L=setdiff(1:10242,mwIndVec_l);
 
% initialize matrix
M = length(nonMW_L); 
N = length(nonMW_L);
P = 100;   % Number of matrices to stack
stacked_matrix = zeros(M, N, P);

% for all n matrices
for i=1:100
	i
	filepath=['/scratch/users/apines/SimStreams/' num2str(i) '_streamConnectivity_L.mat'];
	% test if file exists
	if exist(filepath,'file')
		% load in one
		simMat=load(filepath);
		% mask it
		simMat=simMat.AdjMatrix_L(nonMW_L,nonMW_L);
		% make it float (single?)
		simMat=single(simMat);
		stacked_matrix(:,:,i)=simMat;
	% print something out if it doesnt exist
	else
		fprintf('File %s does not exist.\n', filepath);
	end
end
save('/scratch/users/apines/Left_AllSims.mat','stacked_matrix','-v7.3')
% make sparse
%stacked_sparse_matrix=sparse(stacked_matrix);
% save out
%save('/scratch/users/apines/Left_AllSims_Sparse.mat','stacked_sparse_matrix','-v7.3')
