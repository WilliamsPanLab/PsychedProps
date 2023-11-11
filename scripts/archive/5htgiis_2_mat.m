% convert downsampled network .giis to .mats
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% this participant's filepath
netgiis_L=['~/5ht2aL_3k.func.gii'];
netgiis_R=['~/5ht2aR_3k.func.gii'];
% load in functional networks: Left
Lnets=gifti(netgiis_L);
% load in functional networks: Right
Rnets=gifti(netgiis_R);
% unitary map
Lnets=Lnets.cdata(:,1);
Rnets=Rnets.cdata(:,1);
% save out
nets=struct;
nets.Lnets=Lnets;
nets.Rnets=Rnets;
save(['~/5ht2a_fs4.mat'],'nets')
