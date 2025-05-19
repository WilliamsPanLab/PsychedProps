
%
% The second step of single brain parcellation, clustering of 50 group atlas to create the final group atlas
% For the toolbox of single brain parcellation, see: 
%

clear
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
ProjectFolder = '/oak/stanford/groups/leanew1/users/apines/data/'
resultantFolder = [ProjectFolder '/RobustInitialization20'];
mkdir(resultantFolder);
inFile = [resultantFolder '/ParcelInit_List.txt'];
system(['rm ' inFile]);
AllFiles = g_ls([ProjectFolder '/SingleParcellation/fs5_Surf_Initialization/*/*comp20*/*.mat']);
for i = 1:length(AllFiles)
  cmd = ['echo ' AllFiles{i} ' >> ' inFile];
  system(cmd);
end

% Parcellate into 4 networks
K = 20;
selRobustInit(inFile, K, resultantFolder);
