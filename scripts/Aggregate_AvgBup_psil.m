% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot condition maps
L_bv_list=g_ls(['/scratch/users/apines/data/psil/*/AvgBup_Bf_L.mat']);
R_bv_list=g_ls(['/scratch/users/apines/data/psil/*/AvgBup_Bf_R.mat']);
L_bw_list=g_ls(['/scratch/users/apines/data/psil/*/AvgBup_Bw_L.mat']);
R_bw_list=g_ls(['/scratch/users/apines/data/psil/*/AvgBup_Bw_R.mat']);
L_af_list=g_ls(['/scratch/users/apines/data/psil/*/AvgBup_Af_L.mat']);
R_af_list=g_ls(['/scratch/users/apines/data/psil/*/AvgBup_Af_R.mat']);
L_p_list=g_ls(['/scratch/users/apines/data/psil/*/AvgBup_p_L.mat']);
R_p_list=g_ls(['/scratch/users/apines/data/psil/*/AvgBup_p_R.mat']);
L_m_list=g_ls(['/scratch/users/apines/data/psil/*/AvgBup_m_L.mat']);
R_m_list=g_ls(['/scratch/users/apines/data/psil/*/AvgBup_m_R.mat']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot difference maps
%%% FOR FIGURE 2
%%%%% plot avg. placebo vs. avg drug
% initialize output vectors, load in a template to do so
template_L=load(L_m_list{1});
template_R=load(R_m_list{1});
aggregate_L=zeros(length(L_m_list),length(template_L.avg_L));
aggregate_R=zeros(length(R_m_list),length(template_R.avg_R));
% aggregate
for i=1:length(L_m_list)
	aggregate_L(i,:)=load(L_m_list{i}).avg_L;
	aggregate_R(i,:)=load(R_m_list{i}).avg_R;
end
% average
avg_L=mean(aggregate_L);
avg_R=mean(aggregate_R);
% get average of psil for L and R to average and subtract
% initialize output vectors, load in a template to do so
template_L2=load(L_p_list{1});
template_R2=load(R_p_list{1});
aggregate_L2=zeros(length(L_p_list),length(template_L2.avg_L));
aggregate_R2=zeros(length(R_p_list),length(template_R2.avg_R));
% aggregate drug scans for 80 and 120
for i=1:length(L_p_list)
	aggregate_L2(i,:)=load(L_p_list{i}).avg_L;
	aggregate_R2(i,:)=load(R_p_list{i}).avg_R;
end
% average
avg_L2=mean(aggregate_L2);
avg_R2=mean(aggregate_R2);
% subtract
avg_L_dif=avg_L2-avg_L;
avg_R_dif=avg_R2-avg_R;
% plot avg. no drug vs. avg. drug
Vis_FaceVec(avg_L_dif,avg_R_dif,'~/Methyl_vs_Psilo_BupDif.png')

%%%%% plot avg. no drug vs. avg drug
% initialize output vectors, load in a template to do so
template_L=load(L_bv_list{1});
template_R=load(R_bv_list{1});
aggregate_L=zeros(length(L_bv_list)+length(L_bw_list)+length(L_af_list),length(template_L.avg_L));
aggregate_R=zeros(length(R_bv_list)+length(R_bw_list)+length(R_af_list),length(template_R.avg_R));
% aggregate - BV
for i=1:length(L_bv_list)
        aggregate_L(i,:)=load(L_bv_list{i}).avg_L;
        aggregate_R(i,:)=load(R_bv_list{i}).avg_R;
end
% between
for i=1:length(L_bw_list)
        aggregate_L(i+length(L_bv_list),:)=load(L_bw_list{i}).avg_L;
        aggregate_R(i+length(R_bv_list),:)=load(R_bw_list{i}).avg_R;
end
% after
for i=1:length(L_af_list)
        aggregate_L(i+length(L_bv_list)+length(L_bw_list),:)=load(L_af_list{i}).avg_L;
        aggregate_R(i+length(R_bv_list)+length(L_bw_list),:)=load(R_af_list{i}).avg_R;
end
% average
avg_L=mean(aggregate_L);
avg_R=mean(aggregate_R);
% get average of psil for L and R to average and subtract
% initialize output vectors, load in a template to do so
template_L2=load(L_p_list{1});
template_R2=load(R_p_list{1});
aggregate_L2=zeros(length(L_p_list),length(template_L2.avg_L));
aggregate_R2=zeros(length(R_p_list),length(template_R2.avg_R));
% aggregate drug scans for 80 and 120
for i=1:length(L_p_list)
        aggregate_L2(i,:)=load(L_p_list{i}).avg_L;
        aggregate_R2(i,:)=load(R_p_list{i}).avg_R;
end
% average
avg_L2=mean(aggregate_L2);
avg_R2=mean(aggregate_R2);

% subtract
avg_L_dif=avg_L2-avg_L;
avg_R_dif=avg_R2-avg_R;
% plot avg. no drug vs. avg. drug
Vis_FaceVec(avg_L_dif,avg_R_dif,'~/NoDrug_vs_Drug_BupDif.png')
