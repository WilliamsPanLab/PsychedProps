% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot condition maps
L_pl_list=g_ls(['/scratch/users/apines/data/mdma/*/AvgMag_PL_L.mat']);
R_pl_list=g_ls(['/scratch/users/apines/data/mdma/*/AvgMag_PL_R.mat']);
L_bv_list=g_ls(['/scratch/users/apines/data/mdma/*/AvgMag_BV_L.mat']);
R_bv_list=g_ls(['/scratch/users/apines/data/mdma/*/AvgMag_BV_R.mat']);
L_m1_list=g_ls(['/scratch/users/apines/data/mdma/*/AvgMag_M1_L.mat']);
R_m1_list=g_ls(['/scratch/users/apines/data/mdma/*/AvgMag_M1_R.mat']);
L_m2_list=g_ls(['/scratch/users/apines/data/mdma/*/AvgMag_M2_L.mat']);
R_m2_list=g_ls(['/scratch/users/apines/data/mdma/*/AvgMag_M2_R.mat']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot difference maps
%%% FOR FIGURE 2
%%%%% plot avg. placebo vs. avg drug
% initialize output vectors, load in a template to do so
template_L=load(L_pl_list{1});
template_R=load(R_pl_list{1});
aggregate_L=zeros(length(L_pl_list),length(template_L.PL_L));
aggregate_R=zeros(length(R_pl_list),length(template_R.PL_R));
% aggregate
for i=1:length(L_pl_list)
	aggregate_L(i,:)=load(L_pl_list{i}).PL_L;
	aggregate_R(i,:)=load(R_pl_list{i}).PL_R;
end
% average
avg_L=mean(aggregate_L);
avg_R=mean(aggregate_R);
% get average of 80 and 120 for L and R to average and subtract
% initialize output vectors, load in a template to do so
template_L2=load(L_pl_list{1});
template_R2=load(R_pl_list{1});
aggregate_L2=zeros(length(L_m1_list)+length(L_m2_list),length(template_L2.PL_L));
aggregate_R2=zeros(length(R_m1_list)+length(R_m2_list),length(template_R2.PL_R));
% aggregate drug scans for 80 and 120
for i=1:length(L_m1_list)
	aggregate_L2(i,:)=load(L_m1_list{i}).M1_L;
	aggregate_R2(i,:)=load(R_m1_list{i}).M1_R;
end
% 120
for i=1:length(L_m2_list)
	aggregate_L2(i+length(L_m1_list),:)=load(L_m2_list{i}).M2_L;
	aggregate_R2(i+length(R_m1_list),:)=load(R_m2_list{i}).M2_R;
end
% average
avg_L2=mean(aggregate_L2);
avg_R2=mean(aggregate_R2);
% subtract
avg_L_dif=avg_L2-avg_L;
avg_R_dif=avg_R2-avg_R;
% plot avg. no drug vs. avg. drug
Vis_FaceVec(avg_L_dif,avg_R_dif,'~/PL_vs_Drug_BupDif.png')

%%%%% plot avg. no drug vs. avg drug
% initialize output vectors, load in a template to do so
template_L=load(L_pl_list{1});
template_R=load(R_pl_list{1});
aggregate_L=zeros(length(L_pl_list)+length(L_bv_list),length(template_L.PL_L));
aggregate_R=zeros(length(R_pl_list)+length(R_bv_list),length(template_R.PL_R));
% aggregate
for i=1:length(L_pl_list)
        aggregate_L(i,:)=load(L_pl_list{i}).PL_L;
        aggregate_R(i,:)=load(R_pl_list{i}).PL_R;
end
% baseline
for i=1:length(L_bv_list)
        aggregate_L(i+length(L_pl_list),:)=load(L_bv_list{i}).BV_L;
        aggregate_R(i+length(R_pl_list),:)=load(R_bv_list{i}).BV_R;
end
% average
avg_L=mean(aggregate_L);
avg_R=mean(aggregate_R);

% get average of 80 and 120 for L and R to average and subtract
% initialize output vectors, load in a template to do so
template_L2=load(L_pl_list{1});
template_R2=load(R_pl_list{1});
aggregate_L2=zeros(length(L_m1_list)+length(L_m2_list),length(template_L2.PL_L));
aggregate_R2=zeros(length(R_m1_list)+length(R_m2_list),length(template_R2.PL_R));
% aggregate drug scans for 80 and 120
for i=1:length(L_m1_list)
        aggregate_L2(i,:)=load(L_m1_list{i}).M1_L;
        aggregate_R2(i,:)=load(R_m1_list{i}).M1_R;
end
% 120
for i=1:length(L_m2_list)
        aggregate_L2(i+length(L_m1_list),:)=load(L_m2_list{i}).M2_L;
        aggregate_R2(i+length(R_m1_list),:)=load(R_m2_list{i}).M2_R;
end
% average
avg_L2=mean(aggregate_L2);
avg_R2=mean(aggregate_R2);
% subtract
avg_L_dif=avg_L2-avg_L;
avg_R_dif=avg_R2-avg_R;
% plot avg. no drug vs. avg. drug
Vis_FaceVec(avg_L_dif,avg_R_dif,'~/NoDrug_vs_Drug_BupDif.png')
