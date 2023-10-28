% load in best community assignment to make groups of normative group-derived vectors, save out reference streams
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% best community assignments
bestcomms=load('~/best_solution.mat').best_solution;
% load in normative group-derived vectors of above-expected streams
above = load('/scratch/users/apines/SimStreams/Group_sigStreams_a.mat').sig_Streams_Above;

% get stream directionality
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
        sigStreamRow=above(v,:);
        % initialize vector score
        vectorScore=zeros(1,3);
        % find each vertex with a nonzerovalue in this row
        nonzvs=find(sigStreamRow);
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

% mask vector scores by instances where bestcomms (~=0) exist
commMask=find(bestcomms==0);
% mask it
vectorScores(commMask,:)=0;

% make a medial wall version
MWvectorscores=zeros(2562,3);
MWvectorscores(nonMW_L,1)=vectorScores(:,1);
MWvectorscores(nonMW_L,2)=vectorScores(:,2);
MWvectorscores(nonMW_L,3)=vectorScores(:,3);

% make a full vector with MW zeroed
vertvec=zeros(2562,1);
vertvec(nonMW_L)=bestcomms;
% set comms to 1
vertvec(find(vertvec))=1;

% print out a version
Vis_Surf_n_Vecfield(vertvec,zeros(2562,1),MWvectorscores,zeros(2562,3),'~/bestcomms_vecs.png')
