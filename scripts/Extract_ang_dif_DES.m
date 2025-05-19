% extract network props from each session

% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% get subj list
dataDir='/oak/stanford/groups/leanew1/SHARED_DATASETS/private/connectome/bids/data/derivatives/fmriprep-20.2.3/fmriprep'
% Use the dir function to get a list of folders starting with "sub-CONN"
dataDir = '/scratch/users/apines/data/mdma/';
% Use the dir function to get a list of items starting with "sub-CONN"
items = dir(fullfile(dataDir, 'sub-CONN*'));
% drop those starting with html (redundant)
subjFolders = items([items.isdir]);
subjList = {subjFolders.name};
% for collecting each task
for task=["rs1" "rs2" "emotion" "gambling" "wm"]
	% initialize output file. Remaining TRs and FD along with td props value for the 4 networks
	outDF=cell(length(subjFolders),7);
	% set common fp
	commonFP=['/scratch/users/apines/data/mdma/'];
	% for each subj except 4 6 10
	for s=1:length(subjFolders);
		disp(s)
		% insert subject ID as first column
		outDF(s,1)=subjList(s);
	        bvFP=[commonFP subjList(s) '/ses-BL/' subjList(s) '_ses-BL_' task '_k4_Prop_Feats_gro.csv'];
	        bvFP=strjoin(bvFP,'');
		% load in baseline csvs
		if exist(bvFP,'file')
			bv=readmatrix(bvFP);
	        	outDF(s,2)=num2cell(bv(1,2));
			outDF(s,3)=num2cell(bv(2,2));
			outDF(s,4)=num2cell(bv(3,2));
			outDF(s,5)=num2cell(bv(4,2));
			% load in remaining frames
		        childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/ses-BL' ];
        	        bvCSIfp=[childfp '/' subjList(s) '_ses-BL_task-' task '_ValidSegments_Trunc.txt'];
        	        bvCSIfp=strjoin(bvCSIfp,'');
        	        CSI = importdata(bvCSIfp);
        	        % trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
        	        numTRsVS=num2cell(sum(CSI(:,2)));
                	outDF(s,6)=numTRsVS;	
			% load in FD
			fdfp=['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' subjList(s) '/ses-BL/func/' subjList(s) '_ses-BL_task-' task '_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_qc.csv'];
			% different fp for resting states
			if task=="rs1"
				fdfp=['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' subjList(s) '/ses-BL/func/' subjList(s) '_ses-BL_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_qc.csv'];
			end
			if task=="rs2"
				fdfp=['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' subjList(s) '/ses-BL/func/' subjList(s) '_ses-BL_task-rs_acq-mb_dir-pe1_run-0_space-fsLR_den-91k_qc.csv'];
			end
			fd=readtable(strjoin(fdfp,''));
			outDF(s,7)=num2cell(fd.meanFD);
		end
	end
	% save out matrix
	writecell(outDF,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/' task '_propsMerged_DES.csv'],''))
end
