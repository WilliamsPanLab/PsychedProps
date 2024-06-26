function OpFl_toVerts(subj,sesh,task)
% take angles and interpolate them to vertices. Then combine with BOLD time series.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% load in angles
parentFP=['/scratch/users/apines/data/psil/' subj '/' sesh];
datafp=[parentFP '/' subj '_' sesh '_' task '_OpFl.mat'];
data=load(datafp);
Angles_L=data.us.vf_left;
Angles_R=data.us.vf_right;

%%%%%% load in surfaces
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfL = [SubjectsFolder '/lh.sphere'];
surfR = [SubjectsFolder '/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;

% get length of time series
lenOpFl=size(Angles_L);
lenOpFl=lenOpFl(2);

% initialize vertexwise vectors
vertWise_Vecs_l=zeros(2562,lenOpFl,3);
vertWise_Vecs_r=zeros(2562,lenOpFl,3);

% left hemi: for each vertex, resampling opflow directions (xyz)
for vertInd=1:2562
        % get faces intertwined with this vertex
        [InvolvedFaces_l,~]=find(faces_l==vertInd);
	% for each timepoint
	for tp=1:lenOpFl
		% get corresponding vector fields
		OGVecs_L=Angles_L{tp};
        	% get mean vector directionality at these faces (x y z)
        	VertVec_L=mean(OGVecs_L(InvolvedFaces_l,:),1);
		% insert into new df
		vertWise_Vecs_l(vertInd,tp,:)=VertVec_L;
	end
end
% right hemi: same thing
for vertInd=1:2562
        % get faces intertwined with this vertex
        [InvolvedFaces_r,~]=find(faces_r==vertInd);
        % for each timepoint
        for tp=1:lenOpFl
                % get corresponding vector fields
                OGVecs_R=Angles_R{tp};
                % get mean vector directionality at these faces (x y z)
                VertVec_R=mean(OGVecs_R(InvolvedFaces_r,:),1);
                % insert into new df
                vertWise_Vecs_r(vertInd,tp,:)=VertVec_R;
        end
end

% save out
fpl=[parentFP '/' subj '_' sesh '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
fpr=[parentFP '/' subj '_' sesh '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
save(fpl, 'vertWise_Vecs_l');
save(fpr, 'vertWise_Vecs_r');
