% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
structList=g_ls(['/scratch/users/apines/data/mdma/*/*Bold_Vec_Sync_L.mat']);

%%%%%% pl_drug
% use g_ls to get all pl_drug left
L_pl_drug_list=g_ls(['/scratch/users/apines/data/mdma/*/pl_drug_L.mat']);
R_pl_drug_list=g_ls(['/scratch/users/apines/data/mdma/*/pl_drug_R.mat']);
% initialize output vectors, load in a template to do so
template_L=load(L_pl_drug_list{1});
template_R=load(R_pl_drug_list{1});
aggregate_L=zeros(length(L_pl_drug_list),length(template_L.pl_drug_L));
aggregate_R=zeros(length(R_pl_drug_list),length(template_R.pl_drug_R));
% aggregate
for i=1:length(L_pl_drug_list)
	aggregate_L(i,:)=load(L_pl_drug_list{i}).pl_drug_L;
        aggregate_R(i,:)=load(R_pl_drug_list{i}).pl_drug_R;
end
% average
avg_L=mean(aggregate_L);
avg_R=mean(aggregate_R);
% plot
Vis_Vertvec(avg_L,avg_R,'~/PL_vs_Drug_Distance.png')

%%%%% bv_drug
L_bv_drug_list=g_ls(['/scratch/users/apines/data/mdma/*/bv_drug_L.mat']);
R_bv_drug_list=g_ls(['/scratch/users/apines/data/mdma/*/bv_drug_R.mat']);
% initialize output vectors, load in a template to do so
template_L=load(L_bv_drug_list{1});
template_R=load(R_bv_drug_list{1});
aggregate_L=zeros(length(L_bv_drug_list),length(template_L.bv_drug_L));
aggregate_R=zeros(length(R_bv_drug_list),length(template_R.bv_drug_R));
% aggregate
for i=1:length(L_bv_drug_list)
        aggregate_L(i,:)=load(L_bv_drug_list{i}).bv_drug_L;
        aggregate_R(i,:)=load(R_bv_drug_list{i}).bv_drug_R;
end
% average
avg_L=mean(aggregate_L);
avg_R=mean(aggregate_R);
% plot
Vis_Vertvec(avg_L,avg_R,'~/BV_vs_Drug_Distance.png')

%%%%%% sob_drug
L_sob_drug_list=g_ls(['/scratch/users/apines/data/mdma/*/sob_drug_L.mat']);
R_sob_drug_list=g_ls(['/scratch/users/apines/data/mdma/*/sob_drug_R.mat']);
% initialize output vectors, load in a template to do so
template_L=load(L_sob_drug_list{1});
template_R=load(R_sob_drug_list{1});
aggregate_L=zeros(length(L_sob_drug_list),length(template_L.sob_drug_L));
aggregate_R=zeros(length(R_sob_drug_list),length(template_R.sob_drug_R));
% aggregate
for i=1:length(L_sob_drug_list)
        aggregate_L(i,:)=load(L_sob_drug_list{i}).sob_drug_L;
        aggregate_R(i,:)=load(R_sob_drug_list{i}).sob_drug_R;
end
% average
avg_L=mean(aggregate_L);
avg_R=mean(aggregate_R);
% plot
Vis_Vertvec(avg_L,avg_R,'~/Sob_Drug_Distance.png')

%%%%%% pl_bv
L_pl_bv_list=g_ls(['/scratch/users/apines/data/mdma/*/pl_bv_L.mat']);
R_pl_bv_list=g_ls(['/scratch/users/apines/data/mdma/*/pl_bv_R.mat']);
% initialize output vectors, load in a template to do so
template_L=load(L_pl_bv_list{1});
template_R=load(R_pl_bv_list{1});
aggregate_L=zeros(length(L_pl_bv_list),length(template_L.pl_bv_L));
aggregate_R=zeros(length(R_pl_bv_list),length(template_R.pl_bv_R));
% aggregate
for i=1:length(L_pl_bv_list)
        aggregate_L(i,:)=load(L_pl_bv_list{i}).pl_bv_L;
        aggregate_R(i,:)=load(R_pl_bv_list{i}).pl_bv_R;
end
% average
avg_L=mean(aggregate_L);
avg_R=mean(aggregate_R);
% plot
Vis_Vertvec(avg_L,avg_R,'~/PL_BV_Distance.png')

%%%%% formal t-test of pl drug distances vs pl bv distances
L_pl_bv_list=g_ls(['/scratch/users/apines/data/mdma/*/pl_bv_L.mat']);
R_pl_bv_list=g_ls(['/scratch/users/apines/data/mdma/*/pl_bv_R.mat']);
% initialize output vectors, load in a template to do so
template_L=load(L_pl_bv_list{1});
template_R=load(R_pl_bv_list{1});
plbv_aggregate_L=zeros(length(L_pl_bv_list),length(template_L.pl_bv_L));
plbv_aggregate_R=zeros(length(R_pl_bv_list),length(template_R.pl_bv_R));
% aggregate
for i=1:length(L_pl_bv_list)
        plbv_aggregate_L(i,:)=load(L_pl_bv_list{i}).pl_bv_L;
        plbv_aggregate_R(i,:)=load(R_pl_bv_list{i}).pl_bv_R;
end
% use g_ls to get all pl_drug left
L_pl_drug_list=g_ls(['/scratch/users/apines/data/mdma/*/pl_drug_L.mat']);
R_pl_drug_list=g_ls(['/scratch/users/apines/data/mdma/*/pl_drug_R.mat']);
% initialize output vectors, load in a template to do so
template_L=load(L_pl_drug_list{1});
template_R=load(R_pl_drug_list{1});
pld_aggregate_L=zeros(length(L_pl_drug_list),length(template_L.pl_drug_L));
pld_aggregate_R=zeros(length(R_pl_drug_list),length(template_R.pl_drug_R));
% aggregate
for i=1:length(L_pl_drug_list)
        pld_aggregate_L(i,:)=load(L_pl_drug_list{i}).pl_drug_L;
        pld_aggregate_R(i,:)=load(R_pl_drug_list{i}).pl_drug_R;
end
% initialize tvecs
tvec_L=zeros(1,length(template_L.pl_drug_L));
tvec_R=zeros(1,length(template_R.pl_drug_R));
for L=1:length(tvec_L);
	[Lh,Lp,Lci,LStats]=ttest(plbv_aggregate_L(:,L),pld_aggregate_L(:,L));
	tvec_L(L)=LStats.tstat;
end
for R=1:length(tvec_R);
        [Rh,Rp,Rci,RStats]=ttest(plbv_aggregate_R(:,R),pld_aggregate_R(:,R));
        tvec_R(R)=RStats.tstat;
end 




