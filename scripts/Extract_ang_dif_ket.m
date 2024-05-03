% extract network props from each session

% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% read in subj-session-dose correspondence
subSeshDose=readtable('~/Ketamine_Dose_info_unblinded_first16.csv');
% get subj list
subjPrefix=repmat('sub-P500',14,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14"];
subjList=strcat(subjPrefix,subjSuffix')
% make a cell version for searching
subSeshDoseCell=table2cell(subSeshDose);
isPlacebo = cellfun(@(x) strcmp(x, 'Placebo'), subSeshDoseCell);
isM1 = cellfun(@(x) strcmp(x, '0.05mg/kg'), subSeshDoseCell);
isM2 = cellfun(@(x) strcmp(x, '0.5mg/kg'), subSeshDoseCell);
% for collecting each task
for task=["rs1" "rs2" "emotion" "gambling" "wm"]
	% initialize output table: difference in alff by each network, for 80 and 120, and in group and individualized atlases
	outDF=zeros(14,20);
	% set common fp
	commonFP=['/scratch/users/apines/data/mdma/'];
	% for each subj
	for s=[1 2 4 5 6 7 8 9 10 11 12 13 14]
		disp(s)
	        % get session info
	        seshInfo=subSeshDose{s,2:5};
		% find where bv is 00, where placebo is is 
		bvFP=[commonFP subjList(s) '/ses-00/' subjList(s) '_ses-00_' task '_k4_Prop_Feats_gro.csv'];
	        bvFP=strjoin(bvFP,'');
		% get placebo session (-1 because first row is row number, -1 because session count begins at 0)
		pSesh=find(isPlacebo(s,:)==1)-2;
		m1Sesh=find(isM1(s,:)==1)-2;
		m2Sesh=find(isM2(s,:)==1)-2;
		pFP=[commonFP subjList(s) '/ses-0' num2str(pSesh) '/' subjList(s) '_ses-0' num2str(pSesh) '_' task '_k4_Prop_Feats_gro.csv'];
	        pFP=strjoin(pFP,'');
	        m1FP=[commonFP subjList(s) '/ses-0' num2str(m1Sesh) '/' subjList(s) '_ses-0' num2str(m1Sesh) '_' task '_k4_Prop_Feats_gro.csv'];
	        m1FP=strjoin(m1FP,'');
	        m2FP=[commonFP subjList(s) '/ses-0' num2str(m2Sesh) '/' subjList(s) '_ses-0' num2str(m2Sesh) '_' task '_k4_Prop_Feats_gro.csv'];
	        m2FP=strjoin(m2FP,'');
		% load in baseline csvs
		if exist(bvFP,'file')
		bv=readmatrix(bvFP);
	        outDF(s,1)=bv(1,2);
		outDF(s,2)=bv(2,2);
		outDF(s,3)=bv(3,2);
		outDF(s,4)=bv(4,2);
		% load in remaining frames
	        childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/ses-00'];
                bvCSIfp=[childfp '/' subjList(s) '_ses-00_task-' task '_ValidSegments_Trunc.txt'];
                bvCSIfp=strjoin(bvCSIfp,'');
                CSI = importdata(bvCSIfp);
                % trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
                numTRsVS=sum(CSI(:,2));
                outDF(s,17)=numTRsVS;	
		end
		% load in placebo csvs
		if exist(pFP,'file')
		p=readmatrix(pFP);
		% extract in group consensus atlas: placebo
		outDF(s,5)=p(1,2);
		outDF(s,6)=p(2,2);
		outDF(s,7)=p(3,2);
		outDF(s,8)=p(4,2);
		childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/ses-0' num2str(pSesh) ];
                pCSIfp=[childfp '/' subjList(s) '_ses-0' num2str(m1Sesh) '_task-' task '_ValidSegments_Trunc.txt'];
                pCSIfp=strjoin(pCSIfp,'');
                CSI = importdata(pCSIfp);
                numTRsVS=sum(CSI(:,2));
                outDF(s,18)=numTRsVS;
		end
		% load mdma 1 (80 mg)
		if exist(m1FP,'file')
		m1=readmatrix(m1FP);
		% extact in group consensus atlas: mdma 1
		outDF(s,9)=m1(1,2);
		outDF(s,10)=m1(2,2);
		outDF(s,11)=m1(3,2);
		outDF(s,12)=m1(4,2);
		childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/ses-0' num2str(m1Sesh) ];
                m1CSIfp=[childfp '/' subjList(s) '_ses-0' num2str(m2Sesh) '_task-' task '_ValidSegments_Trunc.txt'];
                m1CSIfp=strjoin(m1CSIfp,'');
                CSI = importdata(m1CSIfp);
                numTRsVS=sum(CSI(:,2));
                outDF(s,19)=numTRsVS;
		end
		% load mdma 2 (120 mg)
		if exist(m2FP,'file')
		m2=readmatrix(m2FP);
		% extract in group consensus atlas: mdma 2
		outDF(s,13)=m2(1,2);
		outDF(s,14)=m2(2,2);
		outDF(s,15)=m2(3,2);
		outDF(s,16)=m2(4,2);
		childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/ses-0' num2str(m2Sesh) ];
                m2CSIfp=[childfp '/' subjList(s) '_ses-0' num2str(m2Sesh) '_task-' task '_ValidSegments_Trunc.txt'];
                m2CSIfp=strjoin(m2CSIfp,'');
                CSI = importdata(m2CSIfp);
                numTRsVS=sum(CSI(:,2));
                outDF(s,20)=numTRsVS;
		end
	end
	% save out matrix
	writematrix(outDF,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/' task '_propsMerged_ket.csv'],''))
end
