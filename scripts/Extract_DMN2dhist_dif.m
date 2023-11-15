% extract network props from each session

% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% read in subj-session-dose correspondence
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');
% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')
% for collecting each task
for task=["rs1" "rs2" "emotion" "gambling" "wm"]
	% initialize output table: difference in alff by each network, for 80 and 120, and in group and individualized atlases
	outDF_bv=zeros(10,10);
	outDF_pl=zeros(10,10);
	outDF_m1=zeros(10,10);
	outDF_m2=zeros(10,10);
	% set common fp
	commonFP=['/scratch/users/apines/data/mdma/'];
	% for each subj except 4 6 10
	for s=[1 2 3 5 7 8 9 11 12 13 14 15 16 17]
		disp(s)
	        % get session info
	        seshInfo=subSeshDose{s,2:5};
		bvFP=[commonFP subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_task-' task '_DMN_2dhist.csv'];
	        bvFP=strjoin(bvFP,'');
		% and grab reverse phase encoding direction nifti
		pFP=[commonFP subjList(s) '/' seshInfo{2} '/' subjList(s) '_' seshInfo{2} '_task-' task '_DMN_2dhist.csv'];
	        pFP=strjoin(pFP,'');
	        m1FP=[commonFP subjList(s) '/' seshInfo{3} '/' subjList(s) '_' seshInfo{3} '_task-' task '_DMN_2dhist.csv'];
	        m1FP=strjoin(m1FP,'');
	        m2FP=[commonFP subjList(s) '/' seshInfo{4} '/' subjList(s) '_' seshInfo{4} '_task-' task '_DMN_2dhist.csv'];
	        m2FP=strjoin(m2FP,'');
		% load in baseline csvs
		if exist(bvFP,'file')
		bv=readmatrix(bvFP);
	        outDF_bv=outDF_bv+bv;
		end
		% load in placebo csvs
		if exist(pFP,'file')
		p=readmatrix(pFP);
		% extract in group consensus atlas: placebo
		outDF_pl=outDF_pl+p;
		end
		% load mdma 1 (80 mg)
		if exist(m1FP,'file')
		m1=readmatrix(m1FP);
		% extact in group consensus atlas: mdma 1
		outDF_m1=outDF_m1+m1;
		end
		% load mdma 2 (120 mg)
		if exist(m2FP,'file')
		m2=readmatrix(m2FP);
		% extract in group consensus atlas: mdma 2
		outDF_m2=outDF_m2+m2;
		end
	end
	% save out matrix
	writematrix(outDF_bv,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/bv_' task '_DMN_2dhist_merged.csv'],''))
	writematrix(outDF_pl,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/pl_' task '_DMN_2dhist_merged.csv'],''))
	writematrix(outDF_m1,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/m1_' task '_DMN_2dhist_merged.csv'],''))
	writematrix(outDF_m2,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/m2_' task '_DMN_2dhist_merged.csv'],''))
end
