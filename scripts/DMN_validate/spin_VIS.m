% add spin test functions
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/spin_test/scripts');
% add basic functions
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% set parent directory
parentfp='/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize'
% load in dmn
DMN_L=gifti([parentfp '/Group_lh_Network_2.func.gii']).cdata; 
DMN_R=gifti([parentfp '/Group_rh_Network_2.func.gii']).cdata; 
% load in medial wall/tsnr mask so spin test can label them as invalid
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs5surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs5surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
DMN_L(logical(mwAndTSNR_L))=100;
DMN_R(logical(mwAndTSNR_R))=100;
% write a table of em for spin test to work with
writetable(table(DMN_L),['~/forVISspins_L.csv'],'WriteVariableNames',0);
writetable(table(DMN_R),['~/forVISspins_R.csv'],'WriteVariableNames',0);
% spin em
SpinPermuFS(['~/forVISspins_L.csv'], ['~/forVISspins_R.csv'], 10000, '/oak/stanford/groups/leanew1/users/apines/VISspins.mat');
