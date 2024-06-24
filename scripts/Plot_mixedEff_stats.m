% load in t-stats
DrugTs_L=readtable('/scratch/users/apines/taskVerts/valid_Drug_Ts_L.csv');
DrugTs_R=readtable('/scratch/users/apines/taskVerts/valid_Drug_Ts_R.csv');
TaskTs_L=readtable('/scratch/users/apines/taskVerts/valid_Task_Ts_L.csv');
TaskTs_R=readtable('/scratch/users/apines/taskVerts/valid_Task_Ts_R.csv');
DrugTaskTs_L=readtable('/scratch/users/apines/taskVerts/valid_DrugTask_Ts_L.csv');
DrugTaskTs_R=readtable('/scratch/users/apines/taskVerts/valid_DrugTask_Ts_R.csv');

% load in masking procedure (snr and MW) to match what stats were ran on
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% add TSNR mask, includes medial wall
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
mw_L=zeros(1,2562);
mw_L(mwAndTSNR_L==1)=1;
mw_R=zeros(1,2562);
mw_R(mwAndTSNR_R==1)=1;
vmwIndVec_l=find(mw_L);
vmwIndVec_r=find(mw_R);
% make medial wall vector: vertices
g_noMW_combined_L=setdiff([1:2562],vmwIndVec_l);
g_noMW_combined_R=setdiff([1:2562],vmwIndVec_r);

% plant t stats into appropriate vertices (still expecting 999's in this iteration)
% Drug plot
leftVec=zeros(1,2562);
rightVec=zeros(1,2562);
leftVec(g_noMW_combined_L)=DrugTs_L{:,1};
rightVec(g_noMW_combined_R)=DrugTs_R{:,1};
Vis_Vertvec(leftVec,rightVec,'~/Drug_Ts.png')
% task
leftVec=zeros(1,2562);
rightVec=zeros(1,2562);
leftVec(g_noMW_combined_L)=TaskTs_L{:,1};
rightVec(g_noMW_combined_R)=TaskTs_R{:,1};
Vis_Vertvec(leftVec,rightVec,'~/Task_Ts.png')
% task drug interaction
leftVec=zeros(1,2562);
rightVec=zeros(1,2562);
leftVec(g_noMW_combined_L)=DrugTaskTs_L{:,1};
rightVec(g_noMW_combined_R)=DrugTaskTs_R{:,1};
Vis_Vertvec(leftVec,rightVec,'~/DrugTaskIntrxn_Ts.png')
