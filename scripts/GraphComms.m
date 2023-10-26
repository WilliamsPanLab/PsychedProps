% Load your data
above = load('/scratch/users/apines/SimStreams/Group_sigStreams_a.mat');
below = load('/scratch/users/apines/SimStreams/Group_sigStreams_b.mat');

adjmat_above = above.sig_Streams_Above;
adjmat_below = below.sig_Streams_Below;

% Combine the z-score matrices
A = adjmat_above - adjmat_below;

% 100 iterations of louvain
num_iterations = 100;
% the number of nodes
num_nodes = size(A, 1); 

% Initialize the co-occurrence matrix
co_occurrence_matrix = zeros(num_nodes);

% for each louvain iteration
for iter = 1:num_iterations
    % Run Louvain community detection
    CommAff = community_louvain(A, 1, [], 'negative_asym');
    
    % Update the co-occurrence matrix based on community assignments
    for i = 1:num_nodes
        for j = 1:num_nodes
            if CommAff(i) == CommAff(j)
                co_occurrence_matrix(i, j) = co_occurrence_matrix(i, j) + 1;
            end
        end
    end
end




% consider a transparent cortical surface with edges between vertices drawn (color-coded by 





% Define a threshold for co-occurrence
threshold = 0.90 * num_iterations;

% Create a vector to store labels for groups of nodes
group_labels = zeros(1, num_nodes);
% group multi-membership
group_mm = zeros(1,num_nodes);
% Initialize group counter
group_counter = 1;

% for each node
for node = 1:num_nodes
	% Count the number of links for the current node
        num_links = sum(co_occurrence_matrix(node, :) >= threshold);
        % if it has enough links to resemble a community
        if num_links >= 20
		% note if the vertex has a membership already
    		if group_labels(node) == 0
        		% Assign the current node to a new group
       			group_labels(node) = group_counter;
	
	            	% Check other nodes for consistency and assign them to the same group
	            	for other_node = (node + 1):num_nodes
	            	    if co_occurrence_matrix(node, other_node) >= threshold
	            	        group_labels(other_node) = group_counter;
	            	    end
	            	end

	            	% Increment the group counter
	            	group_counter = group_counter + 1;
        	end
    	end
end
% note where group labels are still zero
zerolabs=find(group_labels==0);

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
VertVec(nonMW_L)=group_labels;


Vis_Vertvec(VertVec,zeros(1,2562),'~/comms.png')


% figure out a way to multi-label 
