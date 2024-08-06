% add basic functions
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% load the spins
spins=load('/oak/stanford/groups/leanew1/users/apines/DMNspins.mat');
% load the yeo7 dmn (all in fs5)
YeoAtlasFolder = '/oak/stanford/groups/leanew1/users/apines/maps';
% load in SNR/MW map
lh_mwSNR=gifti('/oak/stanford/groups/leanew1/users/apines/fs5surf/lh.Mask_SNR.func.gii').cdata(:,1);
rh_mwSNR=gifti('/oak/stanford/groups/leanew1/users/apines/fs5surf/rh.Mask_SNR.func.gii').cdata(:,1);
% get indices
non_mwSNR_L=find(lh_mwSNR==0);
non_mwSNR_R=find(rh_mwSNR==0);

% Yeo 7 Systems
[~, Label_lh7, ~] = read_annotation([YeoAtlasFolder '/lh.Yeo2011_7Networks_N1000.annot']);
[~, Label_rh7, names_rh7] = read_annotation([YeoAtlasFolder '/rh.Yeo2011_7Networks_N1000.annot']);
NetworkName7 = {'Visual','Motor','DA','VA','Limbic','FP','DM'};
% load in NMF DMN: 
parentfp='/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize';
L_networks=gifti([parentfp '/Group_lh_Network_1.func.gii']).cdata;
R_networks=gifti([parentfp '/Group_rh_Network_1.func.gii']).cdata;
% index mw/SNR out of both
Label_lh7=Label_lh7(non_mwSNR_L);
Label_rh7=Label_rh7(non_mwSNR_R);
L_networks=L_networks(non_mwSNR_L);
R_networks=R_networks(non_mwSNR_R);

% now get DMN indices from masked y7: note code for DMN is 5127885
DMN_indices_L=find(Label_lh7==5127885);
DMN_indices_R=find(Label_rh7==5127885);
NonDMN_indices_L=find(Label_lh7~=5127885);
NonDMN_indices_R=find(Label_rh7~=5127885);

% calculate real t-stat (DMN loadings in yeo7 dmn mask vs. outside of it
DMNloadings_in_yeo7_L=L_networks(DMN_indices_L);
DMNloadings_out_of_yeo7_L=L_networks(NonDMN_indices_L);
DMNloadings_in_yeo7_R=R_networks(DMN_indices_R);
DMNloadings_out_of_yeo7_R=R_networks(NonDMN_indices_R);
[h,p,ci,stats_true_L]=ttest2(DMNloadings_in_yeo7_L,DMNloadings_out_of_yeo7_L);
[h,p,ci,stats_true_R]=ttest2(DMNloadings_in_yeo7_R,DMNloadings_out_of_yeo7_R);

%%% now get permuted tests
% initialize output vectors
stats_perm_L=zeros(1,10000);
stats_perm_R=zeros(1,10000);
for i = 1:10000
	% pull out this rotation
	rotL=spins.bigrotl(i,:);
	rotR=spins.bigrotr(i,:);
	% apply mw mask
	rotL_masked=rotL(non_mwSNR_L);
	rotR_masked=rotR(non_mwSNR_R);
	% extract DMN values within y7 mask
	permDMNloadings_in_yeo7_L=rotL_masked(DMN_indices_L);
	permDMNloadings_out_of_yeo7_L=rotL_masked(NonDMN_indices_L);
	permDMNloadings_in_yeo7_R=rotR_masked(DMN_indices_R);
	permDMNloadings_out_of_yeo7_R=rotR_masked(NonDMN_indices_R);
	% OMIT 100's: NAs, incurred from rotating MW into non MW space
	permDMNloadings_in_yeo7_L(permDMNloadings_in_yeo7_L==100)=0;
	permDMNloadings_in_yeo7_R(permDMNloadings_in_yeo7_R==100)=0;
	permDMNloadings_out_of_yeo7_L(permDMNloadings_out_of_yeo7_L==100)=0;
	permDMNloadings_out_of_yeo7_R(permDMNloadings_out_of_yeo7_R==100)=0;
	% run t-tests
	[h,p,ci,stats_L]=ttest2(permDMNloadings_in_yeo7_L,permDMNloadings_out_of_yeo7_L);
	[h,p,ci,stats_R]=ttest2(permDMNloadings_in_yeo7_R,permDMNloadings_out_of_yeo7_R);
	stats_perm_L(i)=stats_L.tstat;
	stats_perm_R(i)=stats_R.tstat;
	disp(i)
end

% make true stat the 10,001st value in both
stats_perm_L(10001)=stats_true_L.tstat;
stats_perm_R(10001)=stats_true_R.tstat;
% save out
writetable(table(stats_perm_L'),'~/perm_vs_obs_DMN_Y7_L.csv')
writetable(table(stats_perm_R'),'~/perm_vs_obs_DMN_Y7_R.csv')
