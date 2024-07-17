% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

%%%%%% pl_80
% use g_ls to get all pl_80 left
L_pl_drug_list=g_ls(['/scratch/users/apines/data/mdma/*/pl_80_L_Mag.mat']);
R_pl_drug_list=g_ls(['/scratch/users/apines/data/mdma/*/pl_80_R_Mag.mat']);
% initialize output vectors, load in a template to do so
template_L=load(L_pl_drug_list{1});
template_R=load(R_pl_drug_list{1});
aggregate_L=zeros(length(L_pl_drug_list),length(template_L.pl_80_L));
aggregate_R=zeros(length(R_pl_drug_list),length(template_R.pl_80_R));
% aggregate
for i=1:length(L_pl_drug_list)
        aggregate_L(i,:)=load(L_pl_drug_list{i}).pl_80_L;
	aggregate_R(i,:)=load(R_pl_drug_list{i}).pl_80_R;
end
% average
avg_L=mean(aggregate_L);
avg_R=mean(aggregate_R);
% plot
Vis_Vertvec(avg_L,avg_R,'~/PL_vs_80_MagDif.png')

%%%%%% pl_120
% use g_ls to get all pl_120 left
L_pl_drug_list=g_ls(['/scratch/users/apines/data/mdma/*/pl_120_L_Mag.mat']);
R_pl_drug_list=g_ls(['/scratch/users/apines/data/mdma/*/pl_120_R_Mag.mat']);
% initialize output vectors, load in a template to do so
template_L=load(L_pl_drug_list{1});
template_R=load(R_pl_drug_list{1});
aggregate_L=zeros(length(L_pl_drug_list),length(template_L.pl_120_L));
aggregate_R=zeros(length(R_pl_drug_list),length(template_R.pl_120_R));
% aggregate
for i=1:length(L_pl_drug_list)
	aggregate_L(i,:)=load(L_pl_drug_list{i}).pl_120_L;
	aggregate_R(i,:)=load(R_pl_drug_list{i}).pl_120_R;
end
% average
avg_L=mean(aggregate_L);
avg_R=mean(aggregate_R);
% plot
Vis_Vertvec(avg_L,avg_R,'~/PL_vs_120_MagDif.png')

%%%%%% pl_bv
% use g_ls to get all pl_bv left
L_pl_drug_list=g_ls(['/scratch/users/apines/data/mdma/*/pl_bv_L_Mag.mat']);
R_pl_drug_list=g_ls(['/scratch/users/apines/data/mdma/*/pl_bv_R_Mag.mat']);
% initialize output vectors, load in a template to do so
template_L=load(L_pl_drug_list{1});
template_R=load(R_pl_drug_list{1});
aggregate_L=zeros(length(L_pl_drug_list),length(template_L.pl_bv_L));
aggregate_R=zeros(length(R_pl_drug_list),length(template_R.pl_bv_R));
% aggregate
for i=1:length(L_pl_drug_list)
	aggregate_L(i,:)=load(L_pl_drug_list{i}).pl_bv_L;
	aggregate_R(i,:)=load(R_pl_drug_list{i}).pl_bv_R;
end
% average
avg_L=mean(aggregate_L);
avg_R=mean(aggregate_R);
% plot
Vis_Vertvec(avg_L,avg_R,'~/PL_vs_BV_MagDif.png')


