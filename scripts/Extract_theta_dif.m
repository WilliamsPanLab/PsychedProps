% extract theta angles from each session and face

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
	% initialize output table: difference by network, for 80 and 120, and in gro and indy atlases
	outDF=zeros(17,20);
	% initialize 1116 faces + rem TRs for left hemi
	thetas_L_bv=zeros(17,1117);
	thetas_L_pcb=zeros(17,1117);
	thetas_L_m1=zeros(17,1117);
	thetas_L_m2=zeros(17,1117);
	% initialize 996 faces + rem TRs for right hemi
	thetas_R_bv=zeros(17,997);
	thetas_R_pcb=zeros(17,997);
	thetas_R_m1=zeros(17,997);
	thetas_R_m2=zeros(17,997);
	% set common fp
	commonFP=['/scratch/users/apines/data/mdma/'];
	% for each subj except 4 6 10
	for s=[1 2 3 5 7 8 9 11 12 13 14 15 16 17]
		disp(s)
	        % get session info
	        seshInfo=subSeshDose{s,2:5};
	        bvFP=[commonFP subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_' task '_k1_Thetas_L.csv'];
	        bvFP=strjoin(bvFP,'');
		% and grab reverse phase encoding direction nifti
		pFP=[commonFP subjList(s) '/' seshInfo{2} '/' subjList(s) '_' seshInfo{2} '_' task '_k1_Thetas_L.csv'];
	        pFP=strjoin(pFP,'');
	        m1FP=[commonFP subjList(s) '/' seshInfo{3} '/' subjList(s) '_' seshInfo{3} '_' task '_k1_Thetas_L.csv'];
	        m1FP=strjoin(m1FP,'');
	        m2FP=[commonFP subjList(s) '/' seshInfo{4} '/' subjList(s) '_' seshInfo{4} '_' task '_k1_Thetas_L.csv'];
	        m2FP=strjoin(m2FP,'');
		% and right hemisphere fps
		bvRFP=[commonFP subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_' task '_k1_Thetas_R.csv'];
	        bvRFP=strjoin(bvRFP,'');
		pRFP=[commonFP subjList(s) '/' seshInfo{2} '/' subjList(s) '_' seshInfo{2} '_' task '_k1_Thetas_R.csv'];
		pRFP=strjoin(pRFP,'');
		m1RFP=[commonFP subjList(s) '/' seshInfo{3} '/' subjList(s) '_' seshInfo{3} '_' task '_k1_Thetas_R.csv'];
		m1RFP=strjoin(m1RFP,'');
		m2RFP=[commonFP subjList(s) '/' seshInfo{4} '/' subjList(s) '_' seshInfo{4} '_' task '_k4_Thetas_R.csv'];
		m2RFP=strjoin(m2RFP,'');
		% load in baseline csvs
		if exist(bvFP,'file')
		bv=readmatrix(bvFP);
		bvr=readmatrix(bvRFP);
		% put thetas into appropriate vector
		% load in remaining frames
	        childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{1} ];
                bvCSIfp=[childfp '/' subjList(s) '_' seshInfo{1} '_task-' task '_ValidSegments_Trunc.txt'];
                bvCSIfp=strjoin(bvCSIfp,'');
                CSI = importdata(bvCSIfp);
                % trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
                numTRsVS=sum(CSI(:,2));
		thetas_L_bv(s,1:1116)=bv(:,2);
		thetas_L_bv(s,1117)=numTRsVS;
		thetas_R_bv(s,1:996)=bvr(:,2);
		thetas_R_bv(s,997)=numTRsVS;
		end
		% load in placebo csvs
		if exist(pFP,'file')
		p=readmatrix(pFP);
		pr=readmatrix(pRFP);
		% put thetas into appropriate vector
		thetas_L_pcb(s,1:1116)=p(:,2);
		thetas_R_pcb(s,1:996)=pr(:,2);
		% load in remaining frames
		childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{2} ];
                pCSIfp=[childfp '/' subjList(s) '_' seshInfo{2} '_task-' task '_ValidSegments_Trunc.txt'];
                pCSIfp=strjoin(pCSIfp,'');
                CSI = importdata(pCSIfp);
                numTRsVS=sum(CSI(:,2));
		thetas_L_pcb(s,1117)=numTRsVS;
		thetas_R_pcb(s,997)=numTRsVS;
		end
		% load mdma 1 (80 mg)
		if exist(m1FP,'file')
		m1=readmatrix(m1FP);
		m1r=readmatrix(m1RFP);
		% mdma 1
		thetas_L_m1(s,1:1116)=m1(:,2);
		thetas_R_m1(s,1:996)=m1r(:,2);
		% load in remaining frames
		childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{3} ];
                m1CSIfp=[childfp '/' subjList(s) '_' seshInfo{3} '_task-' task '_ValidSegments_Trunc.txt'];
                m1CSIfp=strjoin(m1CSIfp,'');
                CSI = importdata(m1CSIfp);
                numTRsVS=sum(CSI(:,2));
		thetas_L_m1(s,1117)=numTRsVS;
		thetas_R_m1(s,997)=numTRsVS;
		end
		% load mdma 2 (120 mg)
		if exist(m2FP,'file')
		m2=readmatrix(m2FP);
		m2r=readmatrix(m2RFP);
		thetas_L_m2(s,1:1116)=m2(:,2);
		thetas_R_m2(s,1:996)=m2r(:,2);
		% extract in group consensus atlas: mdma 2
		childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{4} ];
                m2CSIfp=[childfp '/' subjList(s) '_' seshInfo{4} '_task-' task '_ValidSegments_Trunc.txt'];
                m2CSIfp=strjoin(m2CSIfp,'');
                CSI = importdata(m2CSIfp);
                numTRsVS=sum(CSI(:,2));
		thetas_L_m2(s,1117)=numTRsVS;
		thetas_R_m2(s,997)=numTRsVS;
		end
	end
	% save out matrix
	writematrix(thetas_L_bv,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/' task '_thetas_L_bv.csv'],''))
	writematrix(thetas_L_pcb,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/' task '_thetas_L_pcb.csv'],''))
	writematrix(thetas_L_m1,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/' task '_thetas_L_m1.csv'],''))
	writematrix(thetas_L_m2,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/' task '_thetas_L_m2.csv'],''))
	writematrix(thetas_R_bv,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/' task '_thetas_R_bv.csv'],''))
	writematrix(thetas_R_pcb,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/' task '_thetas_R_pcb.csv'],''))
	writematrix(thetas_R_m1,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/' task '_thetas_R_m1.csv'],''))
	writematrix(thetas_R_m2,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/' task '_thetas_R_m2.csv'],''))
end
