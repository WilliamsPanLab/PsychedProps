% convert downsampled network .giis to .mats
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% group consensus filepath
funcgiiFolder = ['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/'];
Lnets=zeros(2562,14);
Rnets=zeros(2562,14);
% for each of 14 networks
for k=1:14
	netgiis_L=[funcgiiFolder 'Group_lh_14Network_' num2str(k) '_3k.func.gii'];	
	netgiis_R=[funcgiiFolder 'Group_rh_14Network_' num2str(k) '_3k.func.gii'];	
	% load in functional networks: Left
	Lnetsk=gifti(netgiis_L);
	% load in functional networks: Right
	Rnetsk=gifti(netgiis_R);
	% set to 4 networks
	Lnets(:,k)=Lnetsk.cdata;
	Rnets(:,k)=Rnetsk.cdata;
end
% save out
nets=struct;
nets.Lnets=Lnets;
nets.Rnets=Rnets;
% save out
save([funcgiiFolder 'gro_14Nets_fs4.mat'],'nets')
