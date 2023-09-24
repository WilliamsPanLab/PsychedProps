function Netgiis_2_mat_Ind(subj)
% convert downsampled network .giis to .mats
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% this participant's filepath
funcgiiFolder = ['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subj '/ses-00/pfm'];
netgiis_L=[funcgiiFolder '/' subj '_L_AggNets_10k.func.gii'];
netgiis_R=[funcgiiFolder '/' subj '_R_AggNets_10k.func.gii'];
% load in functional networks: Left
Lnets=gifti(netgiis_L);
% load in functional networks: Right
Rnets=gifti(netgiis_R);
% set to the distance-from-dmn-smoothed maps
Lnets=Lnets.cdata(:,1);
Rnets=Rnets.cdata(:,1);
% save out
nets=struct;
nets.Lnets=Lnets;
nets.Rnets=Rnets;
save([funcgiiFolder '/' subj '_Nets_fs5.mat'],'nets')
