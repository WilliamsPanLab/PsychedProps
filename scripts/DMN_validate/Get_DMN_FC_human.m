function Get_DMN_FC_human(subj,sesh,task) 
% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))

% load in DMN
networks=load(['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/gro_Nets_fs4.mat']);
Dnet_LH=networks.nets.Lnets(:,1);
Dnet_RH=networks.nets.Rnets(:,1);

% use DMN mask and fs4 time series to get average DMN time series
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
% load in data
fpL=[childfp '/' subj '_' sesh '_task-' task '_p2mm_masked_L.mgh'];
fpR=[childfp '/' subj '_' sesh '_task-' task '_p2mm_masked_R.mgh'];
dataL=MRIread(fpL).vol;
dataR=MRIread(fpR).vol;
% squeeze to get rid of extra dimensions
TRs_l_g=squeeze(dataL);
TRs_r_g=squeeze(dataR);
% make DMN mask 
DMNbool_L=Dnet_LH>.3;
DMNbool_R=Dnet_RH>.3;
% get avg DMN time series across both hemis
DMNTS_L=TRs_l_g(DMNbool_L,:);
DMNTS_R=TRs_r_g(DMNbool_R,:);
mergedDMNTS=[DMNTS_L; DMNTS_R];
meanDMNTS=mean(mergedDMNTS);
% calculate average correlation with DMN time series
DMNcor_L = zeros(size(TRs_l_g, 1), 1);
DMNcor_R = zeros(size(TRs_r_g, 1), 1);
% left hemi
for i = 1:size(TRs_l_g, 1)
    DMNcor_L(i) = corr(meanDMNTS', TRs_l_g(i, :)');
end
% right hemi
for i = 1:size(TRs_r_g, 1)
    DMNcor_R(i) = corr(meanDMNTS', TRs_r_g(i, :)');
end

% make output structure
CorStruct={};
CorStruct.DMNcor_L=DMNcor_L;
CorStruct.DMNcor_R=DMNcor_R;

% save out as a 1xn vertices vector for each subject with subj sesh task in fn for merging later
save([childfp '/' subj '_' sesh '_' task '_DMNFC_map.mat'],'CorStruct')
