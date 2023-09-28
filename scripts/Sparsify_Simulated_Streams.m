% to combine simulated streamline matrices and compress them with sparse + single
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% load in medial wall
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/fs5surf';
% use native freesurfer command for mw mask indices
surfML = [SubjectsFolder '/lh.Medial_wall.label'];
mwIndVec_l = read_medial_wall_label(surfML);
% make "isn't medial wall" indices
nonMW_L=setdiff(1:10242,mwIndVec_l);

% for all n matrices
for i=1:100
	i
	filepath=['/scratch/users/apines/SimStreams/' num2str(i) '_streamConnectivity_L.mat'];
	outpath=['/scratch/users/apines/SimStreams/' num2str(i) '_streamConnectivity_L_sparse.mat'];
	% test if file exists
	if exist(filepath,'file')
		% and if sparse matrix doesnt already exist
		if ~exist(outpath, 'file')
			% load in one
			simMat=load(filepath);
			% mask it
			simMat=simMat.AdjMatrix_L(nonMW_L,nonMW_L);
			sparseSimMat=sparse(simMat);
			% save it
			save(['/scratch/users/apines/SimStreams/' num2str(i) '_streamConnectivity_L_sparse.mat'],'sparseSimMat','-v7.3');	
		else
			fprintf('File %s already exists.\n', outpath);
		end
	% print something out if it doesnt exist
	else
		fprintf('File %s does not exist.\n', filepath);
	end
end
