function OpFlStreamlines_Sim(seed)
% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
childfp=['/scratch/users/apines/'];
% load in optical flow output
data=load([childfp 'Simulated_OpFl_rs_fs5_' seed '.mat']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/fs5surf';
% for surface data
surfL = [SubjectsFolder '/lh.sphere'];
surfR = [SubjectsFolder '/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% normalize verts to unit sphere
% left
numV=length(vx_l);
vx_l(numV+1:end, :) = VecNormalize(vx_l(numV+1:end, :));
% right
numV=length(vx_r);
vx_r(numV+1:end, :) = VecNormalize(vx_r(numV+1:end, :));

% calculate normative directionality streamlines for this subject and session
n = size(faces_l, 1);
T = TriRep(faces_l, vx_l);
P_L = T.incenters;
n = size(faces_r, 1);
T = TriRep(faces_r, vx_r);
P_R = T.incenters;
%%% L as left Hemi
L=data.us.vf_left;
%%% R as right
R=data.us.vf_right;
% use native freesurfer command for mw mask indices
surfML = [SubjectsFolder '/lh.Medial_wall.label'];
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = [SubjectsFolder '/rh.Medial_wall.label'];
mwIndVec_r = read_medial_wall_label(surfMR);
% make "isn't medial wall" indices
nonMW_L=setdiff(1:10242,mwIndVec_l);
nonMW_R=setdiff(1:10242,mwIndVec_r);
%%% get length of time serires
tsLength=length(L);

% initialize vector field to be plotted
VF_L=zeros(n,3,tsLength);
VF_R=zeros(n,3,tsLength);
% loop over to extract from unfortunate cell structure
for i=1:tsLength;
	VF_L(:,:,i)=L{i};
	VF_R(:,:,i)=R{i};
end
% correct spatial coordinates to be on scale of vector field
P_L=P_L./100;
P_R=P_R./100;
vx_l=vx_l./100;
vx_r=vx_r./100;
% initialize SC matrix equivalent
AdjMatrix_L=zeros(length(vx_l),length(vx_l));
AdjMatrix_R=zeros(length(vx_r),length(vx_r));
% initialize threeVFs
%threeVFs=zeros(3,3);
pool=parpool('local',9);
% for each vertex
parfor v=1:length(vx_l)
	% print v
	v
	% initialize row for this vertex to index into
	adjacency_row = zeros(1, length(vx_l));
	% for each TR
	for t=1:tsLength
		% plant a new seed
		CurrSeed=vx_l(v,:);
		% set endpoint of tracking: if t + 30 > length of time series, cant track beyond ts
		endpoint=30;
		if (endpoint+t)>tsLength
			endpoint=tsLength-t;
		else
		end
		% for 30 timepoints (unless condition above met)
		for t2=1:endpoint
			% find the three nearest faces
			[nearFaces,nearDistances]=find_nearest_faces(CurrSeed,P_L);
			% get minimum distance for weighting
			minDist=min(nearDistances);
			% get the 3 vector fields, original TR + 30 timepoints (-1 because 1 + 1 is iter 1)
			threeVFs=VF_L(nearFaces,:,(t+t2-1));
			% get a weighting vector
			Wvec1=minDist/nearDistances(1);
			Wvec2=minDist/nearDistances(2);
			Wvec3=minDist/nearDistances(3);
			CorrectedVec1=threeVFs(1,:)*Wvec1;
			CorrectedVec2=threeVFs(2,:)*Wvec2;
			CorrectedVec3=threeVFs(3,:)*Wvec3;	
			% get a weighted average directionality across the 3 faces
			CorrectedAvg=(CorrectedVec1+CorrectedVec2+CorrectedVec3)/(Wvec1+Wvec2+Wvec3);
			% propagate ongoing seeds one step
			CurrSeed=CurrSeed+CorrectedAvg;
			% find which vertex the resultant propagating seed is at
			nearest_vertex=find_nearest_vertex(CurrSeed,vx_l);
			% add them to adjacency row at (v,site of ongoing seed)
			adjacency_row(nearest_vertex)=adjacency_row(nearest_vertex)+1;
		end
	end
	% add adjacency row to adjacency matrix
	AdjMatrix_L(v,:)=adjacency_row;
end
delete(pool);
% save out matrices for this participant
fn=[childfp 'SimStreams/' seed '_streamConnectivity_L.mat'];
save(fn,'AdjMatrix_L','-v7.3');
