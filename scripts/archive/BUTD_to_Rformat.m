function BUTD_to_Rformat(subj)
% convert to csv for R usage, gonogo and nonconfaces

% GNG
% filepath 
outFP_L=['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/' subj '/gng_BUTD_L.mat'];
outFP_R=['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/' subj '/gng_BUTD_R.mat'];
% load in L	
F_L=load(outFP_L);
% load in R
F_R=load(outFP_R);
% printout columns 1 (BUProp), 6 (BU_resvec_R), 7 (TD_resvec_R), and 8 (whole-circle resvec theta)
matL=cell2mat(F_L.OutDf_L);
matR=cell2mat(F_R.OutDf_R);
csv_L=table(matL(:,:));
csv_R=table(matR(:,:));
writetable(csv_L,['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/' subj '/gng_BUTD_L.csv'])
writetable(csv_R,['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/' subj '/gng_BUTD_R.csv'])

% NCF
% filepath 
outFP_L=['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/' subj '/ncf_BUTD_L.mat'];
outFP_R=['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/' subj '/ncf_BUTD_R.mat'];
% load in L     
F_L=load(outFP_L);
% load in R
F_R=load(outFP_R);
% printout columns 1 (BUProp), 6 (BU_resvec_R), 7 (TD_resvec_R), and 8 (whole-circle resvec theta)
matL=cell2mat(F_L.OutDf_L);
matR=cell2mat(F_R.OutDf_R);
csv_L=table(matL(:,:)); 
csv_R=table(matR(:,:));
writetable(csv_L,['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/' subj '/ncf_BUTD_L.csv'])
writetable(csv_R,['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/' subj '/ncf_BUTD_R.csv'])

% CON
% filepath 
outFP_L=['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/' subj '/con_BUTD_L.mat'];
outFP_R=['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/' subj '/con_BUTD_R.mat'];
% load in L     
F_L=load(outFP_L);
% load in R
F_R=load(outFP_R);
% printout columns 1 (BUProp), 6 (BU_resvec_R), 7 (TD_resvec_R), and 8 (whole-circle resvec theta)
matL=cell2mat(F_L.OutDf_L);
matR=cell2mat(F_R.OutDf_R);
csv_L=table(matL(:,:));
csv_R=table(matR(:,:));
writetable(csv_L,['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/' subj '/con_BUTD_L.csv'])
writetable(csv_R,['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/' subj '/con_BUTD_R.csv'])
