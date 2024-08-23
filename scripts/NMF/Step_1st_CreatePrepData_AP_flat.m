% add matlab path for used functions
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder'));
% only specification is neighborhood distance
gNb = constructW_flat(4);
% save gNb into file for later use
prepDataName=['/oak/stanford/groups/leanew1/users/apines/data/FlatMicePrepData.mat'];
save(prepDataName, 'gNb');
