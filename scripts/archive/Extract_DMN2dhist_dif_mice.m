% extract network props from each session

% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% get subj list
subjList = {'m2000','m7507','m7520','m7522','m7589','m7594'};
% set common fp
commonFP=['/scratch/users/apines/data/mouse/'];
% initialize outdf
outDF=zeros(10,10,6);
% for each mouse
for s=1:length(subjList);
	disp(s)
	% for each run
	for run=1:6
		fp=['/scratch/users/apines/p50_mice/proc/20200228/' subjList{s} '_' num2str(run) '_DMN_2dhist.csv'];
		% if file exists
		if exist(fp,'file')
			% load it in with readmatrix
			angmag=readmatrix(fp);			
			% stack it into the master csv
			outDF(:,:,run)=outDF(:,:,run)+angmag;
		% if it doesnt exist
		else
		end
	end
end

% save out
writematrix(outDF(:,:,1),'/oak/stanford/groups/leanew1/users/apines/data/mice_1_DMN_2dhist_merged.csv')	
writematrix(outDF(:,:,2),'/oak/stanford/groups/leanew1/users/apines/data/mice_2_DMN_2dhist_merged.csv')	
writematrix(outDF(:,:,3),'/oak/stanford/groups/leanew1/users/apines/data/mice_3_DMN_2dhist_merged.csv')	
writematrix(outDF(:,:,4),'/oak/stanford/groups/leanew1/users/apines/data/mice_4_DMN_2dhist_merged.csv')	
writematrix(outDF(:,:,5),'/oak/stanford/groups/leanew1/users/apines/data/mice_5_DMN_2dhist_merged.csv')	
writematrix(outDF(:,:,6),'/oak/stanford/groups/leanew1/users/apines/data/mice_6_DMN_2dhist_merged.csv')	
