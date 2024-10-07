function BUTD_to_Rformat_MDMA(subj,sesh)
% convert to csv for R usage
% filepaths
childfp=['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/' subj ];
outFP_L=[childfp '/' subj '_' sesh '_rs_BUTD_L.mat'];
outFP_R=[childfp '/' subj '_' sesh '_rs_BUTD_R.mat'];
% load in L     
F_L=load(outFP_L);
% load in R
F_R=load(outFP_R);
% printout columns 1 (BUProp), 6 (BU_resvec_R), 7 (TD_resvec_R), and 8 (whole-circle resvec theta)
matL=cell2mat(F_L.OutDf_L);
matR=cell2mat(F_R.OutDf_R);
csv_L=table(matL(:,:));
csv_R=table(matR(:,:));
writetable(csv_L,['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/' subj '/' subj '_' sesh '_rs_BUTD_L.csv'])
writetable(csv_R,['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/mdma/' subj '/' subj '_' sesh '_rs_BUTD_R.csv'])

