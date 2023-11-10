% Load data
above = load('/scratch/users/apines/SimStreams/Group_sigStreams_a.mat');
below = load('/scratch/users/apines/SimStreams/Group_sigStreams_b.mat');

adjmat_above = above.sig_Streams_Above;
adjmat_below = below.sig_Streams_Below;

% Combine the z-score matrices
A = adjmat_above - adjmat_below;

% 100 iterations of louvain
num_iterations = 1000;
% the number of nodes
num_nodes = size(A, 1); 

% Initialize the co-occurrence matrix
co_occurrence_matrix = zeros(num_nodes);
% keep track of q across iterations
qs=zeros(1,1000);
commAffs=zeros(num_nodes,1000);
% for each louvain iteration
for iter = 1:num_iterations
    % Run Louvain community detection
    [CommAff,q] = community_louvain(A, 1, [], 'negative_asym');
    qs(iter)=q;
    commAffs(:,iter)=CommAff;
    % Update the co-occurrence matrix based on community assignments
    for i = 1:num_nodes
        for j = 1:num_nodes
            if CommAff(i) == CommAff(j)
                co_occurrence_matrix(i, j) = co_occurrence_matrix(i, j) + 1;
            end
        end
    end
end

% Find the solution with the maximum sum of MI (most consistent with others)
[best_q, best_solution_idx] = max(qs);
best_solution = commAffs(:, best_solution_idx);





% consider a transparent cortical surface with edges between vertices drawn (color-coded by 


% Define a threshold for co-occurrence
threshold = 0.80 * num_iterations;

% Create a vector to store labels for groups of nodes
zerolabs=zeros(num_nodes,1);
% for each node
for node = 1:num_nodes
	% get this nodes community indices in best_solution
	nodeComm=best_solution(node);
	commInds=find(best_solution==nodeComm);
	% get nodal value of co-occurence in communities in this matrix
	nodalCoOc=sum(co_occurrence_matrix(node,commInds));
	% retain thoose only greater than threshold
	if nodalCoOc<(length(commInds)*threshold)
		zerolabs(node)=1;
	end
end
% update winning community structure
best_solution(logical(zerolabs))=0;

% load in a medial wall set of vertices 
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
% for surface data
surfL = [SubjectsFolder '/lh.sphere'];
surfR = [SubjectsFolder '/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% use native freesurfer command for mw mask indices
surfML = [SubjectsFolder '/lh.Medial_wall.label'];
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = [SubjectsFolder '/rh.Medial_wall.label'];
mwIndVec_r = read_medial_wall_label(surfMR);
% make "isn't medial wall" indices
nonMW_L=setdiff(1:2562,mwIndVec_l);
nonMW_R=setdiff(1:2562,mwIndVec_r);

% initialize a visualization vector
VertVec=zeros(1,2562);
% plop in communities
VertVec(nonMW_L)=best_solution;
Vis_Vertvec(VertVec,zeros(1,2562),'~/comms_80_coOc.png')
% save out community structure
save('~/best_solution.mat','best_solution')

