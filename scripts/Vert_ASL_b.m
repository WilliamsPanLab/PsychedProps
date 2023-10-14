function Vert_ASL_b(subj,sesh)
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
% initialize vector location matrix (2562 x 3 or subset of 2562, x 3)
vectorLoc=zeros(length(nonMW_L),3);
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
		% plug into coordinate set
		coordSet(find(nonzvs==v2),:)=eucLocV2;
	end
	% get mean location
	meanLoc=mean(coordSet);
	% populate vector score matrix
	vectorLoc(v,:)=meanLoc;
end
% make a medial wall version
MWvectorlocations=zeros(2562,3);
MWvectorlocations(nonMW_L,:)=vectorLoc;

% saveout vector score matrix
save(['/scratch/users/apines/SimStreams/' subj '_' sesh '_vectorLocs.mat'],'MWvectorlocations')
