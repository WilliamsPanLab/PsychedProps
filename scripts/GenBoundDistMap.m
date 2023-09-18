% add paths for cifti stuff
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% Define the subject ID
subject = 'sub-MDMA001';  % Update with the specific subject ID

% define childfp
childfp=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subject '/ses-00/pfm'];

% Define the file path to the network maps
network_maps_path = [childfp '/Bipartite_PhysicalCommunities+SpatialFiltering.dtseries.nii'];

% cifti-sep the files into giftis for each hemi (func.gii)
% separate hemispheres - left
cmd_l = ['wb_command -cifti-separate ' network_maps_path ' COLUMN -metric CORTEX_LEFT ' childfp '/' subject '_L_AggNets.func.gii']; 
cmd_r = ['wb_command -cifti-separate ' network_maps_path ' COLUMN -metric CORTEX_RIGHT ' childfp '/' subject '_R_AggNets.func.gii'];
system(cmd_l)
system(cmd_r)

% Load the network maps
network_maps_L = gifti([childfp '/' subject '_L_AggNets.func.gii']);
network_maps_R = gifti([childfp '/' subject '_R_AggNets.func.gii']);

% Define the file paths to the surface information
surfL_fp = '~/surfs/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii'
surfR_fp = '~/surfs/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii'
surfL=gifti(surfL_fp);
surfR=gifti(surfR_fp);

% extract vertices 
verticesL = surfL.vertices;
verticesR = surfR.vertices;

% read in network_maps as reference cifti
RefCifti=ft_read_cifti_mod(network_maps_path);

% find cortical vertices on surface cortex (not medial wall)
LH_idx = RefCifti.brainstructure(1:length(surfL.vertices))~=-1;
RH_idx = RefCifti.brainstructure((length(surfL.vertices)+1):(length(surfL.vertices)+length(surfR.vertices)))~=-1;

% set community of interest: consider making this an argument for subjs w/o DMN as 1
community_of_interest = 1;

% get community membership from coarse solution
community_assignment_L = network_maps_L.cdata(:, 1);
% community assignment R
community_assignment_R = network_maps_R.cdata(:, 1);
% parse non-MW items from community affiliation vectors
community_assignment_L=community_assignment_L(LH_idx);
community_assignment_R=community_assignment_R(RH_idx);

% get DMN indices
DMN_L_ind = find(community_assignment_L==1);
DMN_R_ind = find(community_assignment_R==1);
nonDMN_L_ind = find(community_assignment_L~=1);
nonDMN_R_ind = find(community_assignment_R~=1);

% initialize vertex-wise vector depicting minimal distance from edge of DMN (within DMN)
nonDMN_dist_L=zeros(sum(LH_idx),1);
nonDMN_dist_R=zeros(sum(RH_idx),1);

% load in distance matrix to index into
DistanceMat_L=smartload('/scratch/users/apines/DistanceMatrix_Geod_L.mat');
DistanceMat_R=smartload('/scratch/users/apines/DistanceMatrix_Geod_R.mat');
% get out-of-network distances for each DMN vertex
% Left
for i=1:length(DMN_L_ind)
	vertInd=DMN_L_ind(i);
	nonDMN_dist_L(vertInd)=min(DistanceMat_L(vertInd,nonDMN_L_ind));
end
% make it 32k length for easy saving again
nonDMN_dist_L_32k=zeros(32492,1);
nonDMN_dist_L_32k(LH_idx)=nonDMN_dist_L;
% Right
for i=1:length(DMN_R_ind)
        vertInd=DMN_R_ind(i);
        nonDMN_dist_R(vertInd)=min(DistanceMat_R(vertInd,nonDMN_R_ind));
end
% make it 32k length for easy saving again
nonDMN_dist_R_32k=zeros(32492,1);
nonDMN_dist_R_32k(LH_idx)=nonDMN_dist_R;

% replace data in input cifti
network_maps_L.cdata(:,1)=nonDMN_dist_L_32k;
network_maps_R.cdata(:,1)=nonDMN_dist_R_32k;
% save them out
save(network_maps_L,[childfp 'DMN_Dist_L.func.gii']);
save(network_maps_R,[childfp 'DMN_Dist_R.func.gii']);

