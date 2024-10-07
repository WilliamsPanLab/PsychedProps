% convert downsampled network .giis to .mats
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% set network directory
netdir='/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize'
netgii_L=[netdir '/Group_lh_Smooth_1_3k.func.gii'];
netgii_R=[netdir '/Group_rh_Smooth_1_3k.func.gii'];
% load in functional networks: Left
Lnets=gifti(netgii_L);
% load in functional networks: Right
Rnets=gifti(netgii_R);
% set to 14 networks
Lnets=Lnets.cdata;
Rnets=Rnets.cdata;
% save out
nets=struct;
nets.Lnets=Lnets;
nets.Rnets=Rnets;
save([netdir '/Smooth_Nets_fs4.mat'],'nets')
