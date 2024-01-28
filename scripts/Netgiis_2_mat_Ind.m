function Netgiis_2_mat_Ind(subj)
% convert downsampled network .giis to .mats
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% this participant's filepath
funcgiiFolder = ['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/'];
% initialize output
Lnets=zeros(2562,4);
Rnets=zeros(2562,4);
% for each of 4 networks
for k=1:4
	netgiis_L=[funcgiiFolder subj '_lh_Smooth_' num2str(k) '_3k.func.gii'];
	netgiis_R=[funcgiiFolder subj '_rh_Smooth_' num2str(k) '_3k.func.gii'];
	% load in functional networks: Left
	Lnets(:,k)=gifti(netgiis_L).cdata;
	% load in functional networks: Right
	Rnets(:,k)=gifti(netgiis_R).cdata;
	
end
% save out
nets=struct;
nets.Lnets=Lnets;
nets.Rnets=Rnets;
save([funcgiiFolder '/' subj '_Nets_fs4.mat'],'nets')
