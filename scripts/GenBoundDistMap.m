% Define the subject ID
subject = 'sub-MDMA001';  % Update with the specific subject ID

% Define the file path to the network maps (change ses-00 if needed)
network_maps_path = ['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subject '/ses-00/pfm/Bipartite_PhysicalCommunities+SpatialFiltering.dtseries.nii'];

% cifti-sep the files into giftis for each hemi (func.gii)
# separate hemispheres - left
wb_command -cifti-separate $AgNet COLUMN -metric CORTEX_LEFT ${childfp}/${subj}_L_AggNets.func.gii 
# right hemi
wb_command -cifti-separate $AgNet COLUMN -metric CORTEX_RIGHT ${childfp}/${subj}_R_AggNets.func.gii

% Load the network maps
network_maps = read_cifti(network_maps_path);
% get community membership from coarse solution
community_assignment = network_maps.cdata(:, 1); 
% get community assignment L
% community assignment R


% Define the file paths to the surface information
surfL_fp = '~/surfs/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii'
surfR_fp = '~/surfs/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii'
surfL=gifti(surfL_fp);
surfR=gifti(surfR_fp);

% extract vertices 
verticesL = surfL.vertices;

% Calculate the minimal distances for vertices within the community of interest (DMN is community 1)
community_of_interest = 1;
community_vertices = verticesL(community_assignment == community_of_interest, :);
distances_to_boundary = calculate_distances_to_boundary(community_vertices, verticesL);

% Create a map representing the minimal distances
minimal_distances_map = zeros(size(verticesL, 1), 1);
minimal_distances_map(community_assignment == community_of_interest) = distances_to_boundary;

