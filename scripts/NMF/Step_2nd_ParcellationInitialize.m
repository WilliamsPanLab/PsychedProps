
%
% The first step of single brain parcellation, initialization of group atlas
% Each time resample 100 subjects, and repeat 50 times
% For the toolbox of single brain parcellation, see: 
%

clear
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

ProjectFolder = '/oak/stanford/groups/leanew1/users/apines/data/'
ParcellationFolder = [ProjectFolder '/SingleParcellation'];
InitializationFolder = [ParcellationFolder '/fs5_Surf_Initialization'];
mkdir(InitializationFolder);
mkdir([InitializationFolder '/Input']);
SubjectsQuantity = 100;
%SubjectsQuantity = 60;
%RawDataFolder = '/scratch/users/apines/data/mdma/'
RawDataFolder = '/scratch/groups/leanew1/xcpd_outP50_36p_bp/bv_rs'
LeftCell = g_ls([RawDataFolder '/*_L_concat_TS_10k.mgh']);
RightCell = g_ls([RawDataFolder '/*_R_concat_TS_10k.mgh']);
prepDataFile = ['/oak/stanford/groups/leanew1/users/apines/data/SingleParcellation/CreatePrepData.mat'];

SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/fs5surf/'

% for surface data
surfL = [SubjectsFolder '/lh.pial'];
surfR = [SubjectsFolder '/rh.pial'];
surfML = [SubjectsFolder '/lh.Mask_SNR.label'];
surfMR = [SubjectsFolder '/rh.Mask_SNR.label'];


spaR = 1;
vxI = 1;
ard = 0;
iterNum = 2000;
%tNum = 423; % number of time points
tNum=1700
alpha = 1;
beta = 10;
resId = 'Initialization';
K = 20; % numer of networks

% Repeat 50 times
for i = 1:50
  disp(num2str(i));
  % john cena'd 
  ResultFile_check = dir([InitializationFolder, '/InitCANTSEEMEializationRes_', num2str(i), '/*num', num2str(SubjectsQuantity),'_comp' num2str(K) '*/ini*.mat']);
  if isempty(ResultFile_check)
    SubjectsIDs = randperm(length(LeftCell), SubjectsQuantity);
    sbjListFile = [InitializationFolder '/Input/sbjListFile_' num2str(i) '.txt'];
    system(['rm -rf ' sbjListFile]);
    for j = 1:length(SubjectsIDs)
      cmd = ['echo ' LeftCell{SubjectsIDs(j)} ' >> ' sbjListFile];
      system(cmd);
      cmd = ['echo ' RightCell{SubjectsIDs(j)} ' >> ' sbjListFile];
      system(cmd);
    end    
    
    outDir = [InitializationFolder '/InitializationRes_' num2str(i)];

    % copy template batch submission to append
    ScriptPathTemplate=[InitializationFolder '/tmp' num2str(i) '.sh'];
    Script_cp = ['cp /home/users/apines/sbatchTemplate.sh ' ScriptPathTemplate];
    system(Script_cp)
    ScriptPath = [InitializationFolder '/temp' num2str(i) '.sh'];

    save([InitializationFolder '/Configuration_' num2str(i) '.mat'], 'sbjListFile', 'surfL', 'surfR', 'surfML', 'surfMR', 'prepDataFile', 'outDir', ...
          'spaR', 'vxI', 'ard', 'iterNum', 'K', 'tNum', 'alpha', 'beta', 'resId');
    cmd = ['matlab -nosplash -nodesktop -r ' ...
          '"addpath(genpath(''' ToolFolder ''')),load(''' ...
          InitializationFolder '/Configuration_' num2str(i) '.mat''),deployFuncInit_surf_fs(sbjListFile, surfL, surfR, ' ...
          'surfML, surfMR, prepDataFile, outDir, spaR, vxI, ard, iterNum, K, tNum, alpha, beta, resId),exit(1)">"' ...
          InitializationFolder '/ParcelInit' num2str(i) '.log" 2>&1'];
    fid = fopen(ScriptPath, 'w');
    fprintf(fid, cmd);
    % copy matlab command into sbatch template
    command=['cat ' ScriptPath '>>' ScriptPathTemplate];
    system(command);
    system(['sbatch ' ScriptPathTemplate]);
%        pause(7)
 end
end

