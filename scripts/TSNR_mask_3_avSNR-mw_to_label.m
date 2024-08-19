% add paths
ToolFolder = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder))
% load in medial wall
% left
left_MW=read_label([],'/oak/stanford/groups/leanew1/users/apines/fs5surf/lh.Medial_wall.label');
% right
right_MW=read_label([],'/oak/stanford/groups/leanew1/users/apines/fs5surf/rh.Medial_wall.label');
% load in fs5 snr
% left
left_SNR=gifti('/oak/stanford/groups/leanew1/users/apines/maps/mdma_L_AggTSNR_10k.func.gii').cdata(:,1);
% right
right_SNR=gifti('/oak/stanford/groups/leanew1/users/apines/maps/mdma_R_AggTSNR_10k.func.gii').cdata(:,1);
% Combine SNR vectors to deduce 15th percentile for removal
combined_SNR = [left_SNR; right_SNR];
% Calculate 15th percentile
percentile_15 = prctile(combined_SNR, 15);
% Set values below 15th percentile to 0
left_SNR(left_SNR<percentile_15)=0;
right_SNR(right_SNR<percentile_15)=0;
% use read_surf to get regular pial surf indices, xyz
RemoveIndex_L = find(~left_SNR);
RemoveIndex_R = find(~right_SNR);
% get surface data
surfL = ['/oak/stanford/groups/leanew1/users/apines/fs5surf/lh.pial'];
[vx_L, ~] = read_surf(surfL);
surfR = ['/oak/stanford/groups/leanew1/users/apines/fs5surf/rh.pial'];
[vx_R, ~] = read_surf(surfR);
xyz_Coor_L = vx_L(RemoveIndex_L, :);
xyz_Coor_R = vx_R(RemoveIndex_R, :);
RemoveIndex_L = RemoveIndex_L - 1;
RemoveIndex_R = RemoveIndex_R - 1;
% set output filepath
Mask_L = ['/oak/stanford/groups/leanew1/users/apines/fs5surf/lh.Mask_SNR.label'];
Mask_R = ['/oak/stanford/groups/leanew1/users/apines/fs5surf/rh.Mask_SNR.label'];
write_label(RemoveIndex_L, xyz_Coor_L, zeros(length(RemoveIndex_L), 1), Mask_L);
write_label(RemoveIndex_R, xyz_Coor_R, zeros(length(RemoveIndex_R), 1), Mask_R);
% same but for a gifti
% use left and right snr as templates
left_SNR=gifti('/oak/stanford/groups/leanew1/users/apines/maps/mdma_L_AggTSNR_10k.func.gii');
right_SNR=gifti('/oak/stanford/groups/leanew1/users/apines/maps/mdma_R_AggTSNR_10k.func.gii');
left_SNR.cdata(1:10242,1)=0;
left_SNR.cdata(RemoveIndex_L+1,1)=1;
right_SNR.cdata(1:10242,1)=0;
right_SNR.cdata(RemoveIndex_R+1,1)=1;
% save out
Mask_L = ['/oak/stanford/groups/leanew1/users/apines/fs5surf/lh.Mask_SNR.func.gii'];
Mask_R = ['/oak/stanford/groups/leanew1/users/apines/fs5surf/rh.Mask_SNR.func.gii'];
save(left_SNR,Mask_L);
save(right_SNR,Mask_R);
