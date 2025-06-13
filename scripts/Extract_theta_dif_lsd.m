% extract network props from each session

%%% resultant csv to look like this:
% Subj | Condition | BUP | Task | Remaining Frames | FD 

% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% get subj list
subjPrefix=repmat('S',20,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20"];
subjList=strcat(subjPrefix,subjSuffix')
tasks=["rs1" "rs2" "mus"];
% for collecting each task
for t =1:3
	task=tasks(t)
	% initialize a vector of 1116 faces + rem TRs + FD for both conditions
	thetas_L_PCB=cell(20,1118);
	thetas_L_LSD=cell(20,1118);
	% right is 996 + remTRs + FD
	thetas_R_PCB=cell(20,998);
	thetas_R_LSD=cell(20,998);
	task=tasks(t);
	% extra columns just for subj name and session
	SubjNameCol=cell(120,2);
	% set common fp
	commonFP=['/scratch/users/apines/LSD_ICL/rest_proc/'];
	% do rec. exclusions in R later
	for s=1:20;
		disp(s)
		% populate subject name column
		SubjNameCol((s+(t*20))-20,1)=cellstr(subjList(s));
		SubjNameCol((s+(t*20))-20,2)=cellstr(num2str(t));
	        % grab PCB
		bvFP=[commonFP subjList(s) '/' subjList(s) '_PCB_' task '_k1_Thetas_L.csv'];
	        bvFP=strjoin(bvFP,'');
		% right hemisphere
		bvRFP=[commonFP subjList(s) '/' subjList(s) '_PCB_' task '_k1_Thetas_R.csv'];
	        bvRFP=strjoin(bvRFP,'');
		% m1 = after
		m1FP=[commonFP subjList(s) '/' subjList(s) '_LSD_' task '_k1_Thetas_L.csv'];
	        m1FP=strjoin(m1FP,'');
		% m1 right hemisphere
		m1RFP=[commonFP subjList(s) '/' subjList(s) '_LSD_' task '_k1_Thetas_R.csv'];
	        m1RFP=strjoin(m1RFP,'');
		% load in baseline csvs
		if exist(bvFP,'file')
			bv=readmatrix(bvFP);
			bvR=readmatrix(bvRFP);
			% s+i*20-20 to insert into subj-sesh row 
			thetas_L_PCB(s,1:1116)={round(bv(:,2),5)};
			thetas_R_PCB(s,1:996)={round(bvR(:,2),5)};
			% load in remaining frames
                	bvCSIfp=[commonFP subjList(s) '/' subjList(s) '_PCB_task-' task '_ValidSegments_Trunc.txt'];
                	bvCSIfp=strjoin(bvCSIfp,'');
                	CSI = importdata(bvCSIfp);
                	% get total numer of TRs
			numTRsVS=sum(CSI(:,2));
			% place into out df structure (needs to be 11 *iteration rows down to keep subj-sesh correspondence)
			thetas_L_PCB(s,1117)={numTRsVS};
			thetas_R_PCB(s,997)={numTRsVS};
			% load in FD
			confFilePath=['/scratch/users/apines/LSD_ICL/rest_proc/' subjList(s) '_PCB_bold_signal_properties.mat'];
			conf1=load(strjoin(confFilePath,''));
			if task=='rs1'
           		     scanInds=1:220;
        		elseif task=='rs2'
        		     scanInds=441:660;
        		elseif task=='mus'
        		     scanInds=221:440;
       			end
        		% extract FD columns
			FD=conf1.FD(scanInds);
			% put mean FD in output
			thetas_L_PCB(s,1118)={mean(FD)};
			thetas_R_PCB(s,998)={mean(FD)};
		else
			disp([bvFP ' not found'])
		end
		% load LSD
		if exist(m1FP,'file')
			m1=readmatrix(m1FP);
			m1R=readmatrix(m1RFP);
                        % s+i*20-20 to insert into subj-sesh row
			thetas_L_LSD(s,1:1116)={round(m1(:,2),5)};
			thetas_R_LSD(s,1:996)={round(m1R(:,2),5)};
                        % load in remaining frames
                        m1CSIfp=[commonFP subjList(s) '/' subjList(s) '_LSD_task-' task '_ValidSegments_Trunc.txt'];
                        m1CSIfp=strjoin(m1CSIfp,'');
                        CSI = importdata(m1CSIfp);
                        % get total numer of TRs
                        numTRsVS=sum(CSI(:,2));
                        % place into out df structure (needs to be 11 *iteration rows down to keep subj-sesh correspondence)
			thetas_L_LSD(s,1117)={numTRsVS};
			thetas_R_LSD(s,997)={numTRsVS};
                        % load in FD
                        confFilePath=['/scratch/users/apines/LSD_ICL/rest_proc/' subjList(s) '_LSD_bold_signal_properties.mat'];
                        conf1=load(strjoin(confFilePath,''));
                        if task=='rs1'
                             scanInds=1:220;
                        elseif task=='rs2'
                             scanInds=441:660;
                        elseif task=='mus'
                             scanInds=221:440;
                        end
                        % extract FD columns
                        FD=conf1.FD(scanInds);
                        % put mean FD in output
			thetas_L_LSD(s,1118)={mean(FD)};
			thetas_R_LSD(s,998)={mean(FD)};
		else
                        disp([m1FP ' not found'])
		end
	end
	% save out tables for this task
	thetas_L_PCB_T=cell2table(thetas_L_PCB);
	thetas_L_LSD_T=cell2table(thetas_L_LSD);
	thetas_R_PCB_T=cell2table(thetas_R_PCB);
	thetas_R_LSD_T=cell2table(thetas_R_LSD);
	writetable(thetas_L_PCB_T,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_' task '_PCB_L.csv'],''));
	writetable(thetas_L_LSD_T,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_' task '_LSD_L.csv'],''));
	writetable(thetas_R_PCB_T,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_' task '_PCB_R.csv'],''));
	writetable(thetas_R_LSD_T,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/lsd_thetas_' task '_LSD_R.csv'],''));
end
% and subject ID column
writetable(table(SubjNameCol),['/oak/stanford/groups/leanew1/users/apines/data/lsd_propsMerged_subjOrder.csv'])
