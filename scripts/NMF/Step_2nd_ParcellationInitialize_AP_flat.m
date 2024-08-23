function Step_2nd_ParcellationInitialize_AP_flat()

%clear
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

ProjectFolder = '/oak/stanford/groups/leanew1/users/apines/data/'
ParcellationFolder = [ProjectFolder '/SingleParcellation'];
InitializationFolder = [ParcellationFolder '/FlatMice_Initialization'];
mkdir(InitializationFolder);
mkdir([InitializationFolder '/Input']);
SubjectsQuantity = 30;

RawDataFolder= '/scratch/users/apines/p50_mice/proc/'
CiftiCell = g_ls([RawDataFolder '*/*pre*_1/masked_dff_DS_BP_Smoothed.h5'])

prepDataFile=['/oak/stanford/groups/leanew1/users/apines/data/FlatMicePrepData.mat'];

% set NMF parameters
spaR = 1;
vxI = 1;
ard = 0;
iterNum = 1000;
tNum = 893 ; % number of time points - avg from just resting state scans after .2 mm thresh masking
alpha = 2.2;
beta = 220;
resId = 'Initialization';
K = 13;

% Repeat 50 times
for i = 1:50
  disp(num2str(i));
  % john cena'd (this is a bad joke, all it means is I set the file to find to a nonexistant file so it always passes this check)
  ResultFile_check = dir([InitializationFolder, '/InitializationCANTSEEMERes_', num2str(i), '/*num', num2str(SubjectsQuantity),'_comp' num2str(K) '*/ini*.mat']);
  if isempty(ResultFile_check)
    SubjectsIDs = randperm(length(CiftiCell), SubjectsQuantity);
    sbjListFile = [InitializationFolder '/Input/sbjListFile_mice' num2str(i) '.txt'];
    system(['rm -rf ' sbjListFile]);
    for j = 1:length(SubjectsIDs)
      cmd = ['echo ' CiftiCell{SubjectsIDs(j)} ' >> ' sbjListFile];
      system(cmd);
    end

    outDir = [InitializationFolder '/InitializationRes_' num2str(i)];
    
    % copy template batch submission to append
    ScriptPathTemplate=[InitializationFolder '/tmp' num2str(i) '.sh'];
    Script_cp = ['cp /home/users/apines/sbatchTemplate.sh ' ScriptPathTemplate];
    system(Script_cp)
    ScriptPath = [InitializationFolder '/temp_' num2str(K) '_' num2str(i) '.sh'];

    save([InitializationFolder '/Configuration_' num2str(i) '.mat'], 'sbjListFile', 'prepDataFile', 'outDir', ...
          'spaR', 'vxI', 'ard', 'iterNum', 'K', 'tNum', 'alpha', 'beta', 'resId');
    cmd = ['matlab -nosplash -nodesktop -r ' ...
          '"addpath(genpath(''' ToolFolder ''')),load(''' ...
          InitializationFolder '/Configuration_' num2str(i) '.mat''),deployFuncInit_flat(sbjListFile, ' ...
          'prepDataFile, outDir, spaR, vxI, ard, iterNum, K, tNum, alpha, beta, resId),exit(1)">"' ...
          InitializationFolder '/ParcelInit_' num2str(K) '_'  num2str(i) '.log" 2>&1'];
    fid = fopen(ScriptPath, 'w');
    fprintf(fid, cmd);
    % copy matlab command into sbatch template
    command=['cat ' ScriptPath '>>' ScriptPathTemplate];
    system(command);
    system(['sbatch ' ScriptPathTemplate]);
 	pause(8)
 end
end

