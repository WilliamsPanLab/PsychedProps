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

% plot placebo
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
% plot
Vis_Vertvec(avg_L,avg_R,'~/PL_Mag.png')

% plot baseline
template_L=load(L_bv_list{1});
template_R=load(R_bv_list{1});
aggregate_L=zeros(length(L_bv_list),length(template_L.BV_L));
aggregate_R=zeros(length(R_bv_list),length(template_R.BV_R));
% aggregate
for i=1:length(L_pl_list)
        aggregate_L(i,:)=load(L_bv_list{i}).BV_L;
        aggregate_R(i,:)=load(R_bv_list{i}).BV_R;
end
% average
avg_L=mean(aggregate_L);
avg_R=mean(aggregate_R);
% plot
Vis_Vertvec(avg_L,avg_R,'~/BV_Mag.png')

% plot m1
template_L=load(L_m1_list{1});
template_R=load(R_m1_list{1});
aggregate_L=zeros(length(L_m1_list),length(template_L.M1_L));
aggregate_R=zeros(length(R_m1_list),length(template_R.M1_R));
% aggregate
for i=1:length(L_bv_list)
        aggregate_L(i,:)=load(L_m1_list{i}).M1_L;
        aggregate_R(i,:)=load(R_m1_list{i}).M1_R;
end
% average
avg_L=mean(aggregate_L);
avg_R=mean(aggregate_R);
% plot
Vis_Vertvec(avg_L,avg_R,'~/M1_Mag.png')


% plot m2
template_L=load(L_m2_list{1});
template_R=load(R_m2_list{1});
aggregate_L=zeros(length(L_m2_list),length(template_L.M2_L));
aggregate_R=zeros(length(R_m2_list),length(template_R.M2_R));
% aggregate
for i=1:length(L_m2_list)
        aggregate_L(i,:)=load(L_m2_list{i}).M2_L;
        aggregate_R(i,:)=load(R_m2_list{i}).M2_R;
end
% average
avg_L=mean(aggregate_L);
avg_R=mean(aggregate_R);
% plot
Vis_Vertvec(avg_L,avg_R,'~/M2_Mag.png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot difference maps

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
Vis_Vertvec(avg_L_dif,avg_R_dif,'~/PL_vs_Drug_MagDif.png')

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
Vis_Vertvec(avg_L_dif,avg_R_dif,'~/NoDrug_vs_Drug_MagDif.png')
