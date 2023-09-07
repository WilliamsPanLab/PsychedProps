function Netgiis_2_mat(subj)
% convert downsampled network .giis to .mats
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% this participant's filepath
funcgiiFolder = ['/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/'];
netgiis_L=[funcgiiFolder '' subj '_L_AggNets_10k.func.gii'];
netgiis_R=[funcgiiFolder '' subj '_R_AggNets_10k.func.gii'];
% load in functional networks: Left
Lnets=gifti(netgiis_L);
% load in functional networks: Right
Rnets=gifti(netgiis_R);
% set to 4 networks
Lnets=Lnets.cdata(:,1:4);
Rnets=Rnets.cdata(:,1:4);
% save out
nets=struct;
nets.Lnets=Lnets;
nets.Rnets=Rnets;
save([funcgiiFolder '' subj '_Nets_fs5.mat'],'nets')