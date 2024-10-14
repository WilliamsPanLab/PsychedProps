
%
% The second step of single brain parcellation, clustering of 50 group atlas to create the final group atlas
% For the toolbox of single brain parcellation, see: 
%

clear
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
ProjectFolder = '/oak/stanford/groups/leanew1/users/apines/data'
resultantFolder = [ProjectFolder '/RobustInitialization'];
%mkdir(resultantFolder);
inFile = [resultantFolder '/ParcelInit_List.txt'];
system(['rm ' inFile]);
% select output files
AllFiles = g_ls([ProjectFolder '/SingleParcellation/FlatMice_Initialization/*/Initialization_num30_comp13_S1_4122_L_2576_spaR_1_vxInfo_1_ard_0/*.mat']);
for i = 1:length(AllFiles)
  cmd = ['echo ' AllFiles{i} ' >> ' inFile];
  system(cmd);
end
% Parcellate into 13 networks
K = 13;
selRobustInit(inFile, K, resultantFolder);
