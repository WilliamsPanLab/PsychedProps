function Vert_ASL(subj,sesh)
% compute the vertex-wise anterior/superior/left vector from SigStreams
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% load in sigStreams
sigStreams=load(['/scratch/users/apines/SimStreams/' subj '_' sesh '_sigStreams_b.mat']).sig_Streams_Below;
% load in surface for euclidean locations
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfL = [SubjectsFolder '/lh.inflated'];
[vx_l, faces_l] = read_surf(surfL);
% medial wall
surfML = [SubjectsFolder '/lh.Medial_wall.label'];
mwIndVec_l = read_medial_wall_label(surfML);
% make "isn't medial wall" indices
nonMW_L=setdiff(1:2562,mwIndVec_l);
% initialize vector score matrix (2562 x 3 or subset of 2562, x 3)
vectorScores=zeros(length(nonMW_L),3);
% for each vertex
for v=1:length(vx_l(nonMW_L))
	% get euclidean location of this vertex (x,y,z)
	eucLocV=vx_l(nonMW_L(v),:);
	% get row of interest from sigstreams
	sigStreamRow=sigStreams(v,:);
	% initialize vector score
	vectorScore=zeros(1,3);
	% find each vertex with a nonzerovalue in this row
	nonzvs=find(sigStreamRow);
	% initialize a coordinate set to get the mean location of
	coordSet=zeros(length(nonzvs),3);
	for v2=nonzvs
		% get current cell
		currcell=sigStreamRow(v2);
		% get euclidean location of cell
		eucLocV2=vx_l(nonMW_L(v2),:);
		% get vector from v to v2
		vectorV2V=eucLocV2-eucLocV;
		% multiply by counts, can later ammend to sqrt of counts
		vectorV2V=vectorV2V*currcell;
		% add derived vector to vector score
		vectorScore=vectorScore+vectorV2V;
	end
	% populate vector score matrix
	vectorScores(v,:)=vectorScore;
end
% make a medial wall version
MWvectorscores=zeros(2562,3);
MWvectorscores(nonMW_L,:)=vectorScores;

% saveout vector score matrix
save(['/scratch/users/apines/SimStreams/' subj '_' sesh '_vectorScores.mat'],'MWvectorscores')
end
