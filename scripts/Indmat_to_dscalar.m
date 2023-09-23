function mat_to_dscalar(subj)

% add matlab paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));

% child filepath of genbounddist
childfp=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subject '/ses-00/pfm'];
NetsLeft=[childfp '/DMN_Dist_L-smoothed.func.gii'];
NetsRight=[childfp '/DMN_Dist_R-smoothed.func.gii'];

% resampled V of interest (group or individ - CHANGE AS FIT TO MATCH NAME)
LH=gifti(NetsLeft);
RH=gifti(NetsRight);

%%%%% convert to 1:59412 vector with medial wall information
% Define the file paths to the surface information
surfL_fp = '~/surfs/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii';
surfR_fp = '~/surfs/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii';
surfL=gifti(surfL_fp);
surfR=gifti(surfR_fp);
% extract vertices 
verticesL = surfL.vertices;
verticesR = surfR.vertices;
% Define the file path to the network maps
network_maps_path = [childfp '/Bipartite_PhysicalCommunities+SpatialFiltering.dtseries.nii'];
RefCifti=ft_read_cifti_mod(network_maps_path);
% find cortical vertices on surface cortex (not medial wall)
LH_idx = RefCifti.brainstructure(1:length(surfL.vertices))~=-1;
RH_idx = RefCifti.brainstructure((length(surfL.vertices)+1):(length(surfL.vertices)+length(surfR.vertices)))~=-1;

% extract from struct
V=V{:};

% cifti to replace cdata in
HP=read_cifti('/oak/stanford/groups/leanew1/users/apines/maps/hcp.gradients.dscalar.nii');
HP.cdata(1:59412,1:4)=V(1:59412,1:4);
outputfile=[childfp '/DMNDist.dscalar.nii'];
write_cifti(HP,outputfile)
