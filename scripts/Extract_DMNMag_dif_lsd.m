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
% initialize outdf
outDF=cell(120,6);
% for collecting each task
for t =1:3
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
		bvFP=[commonFP subjList(s) '/' subjList(s) '_PCB_' task '_DMNMag.csv'];
	        bvFP=strjoin(bvFP,'');
		% m1 = after
		m1FP=[commonFP subjList(s) '/' subjList(s) '_LSD_' task '_DMNMag.csv'];
	        m1FP=strjoin(m1FP,'');
		% load in baseline csvs
		if exist(bvFP,'file')
			bv=readmatrix(bvFP);
			outDF(((s+(t*20))-20),1)={subjList(s)};
			outDF(((s+(t*20))-20),2)={"PCB"};
			% s+i*20-20 to insert into subj-sesh row 
			outDF(((s+(t*20))-20),3)={bv(1,2)};
			% task
			outDF(((s+(t*20))-20),4)={task};
			% load in remaining frames
                	bvCSIfp=[commonFP subjList(s) '/' subjList(s) '_PCB_task-' task '_ValidSegments_Trunc.txt'];
                	bvCSIfp=strjoin(bvCSIfp,'');
                	CSI = importdata(bvCSIfp);
                	% get total numer of TRs
			numTRsVS=sum(CSI(:,2));
			% place into out df structure (needs to be 11 *iteration rows down to keep subj-sesh correspondence)
                	outDF(((s+(t*20))-20),5)={numTRsVS};
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
			outDF(((s+(t*20))-20),6)={mean(FD)};		
		else
			disp([bvFP ' not found'])
		end
		% load LSD
		if exist(m1FP,'file')
			m1=readmatrix(m1FP);
			outDF(((s+(t*20)+60)-20),1)={subjList(s)};
			outDF(((s+(t*20)+60)-20),2)={"LSD"};
                        % s+i*20-20 to insert into subj-sesh row
                        outDF(((s+(t*20)+60)-20),3)={m1(1,2)};
                        % task
                        outDF(((s+(t*20)+60)-20),4)={task};
                        % load in remaining frames
                        m1CSIfp=[commonFP subjList(s) '/' subjList(s) '_LSD_task-' task '_ValidSegments_Trunc.txt'];
                        m1CSIfp=strjoin(m1CSIfp,'');
                        CSI = importdata(m1CSIfp);
                        % get total numer of TRs
                        numTRsVS=sum(CSI(:,2));
                        % place into out df structure (needs to be 11 *iteration rows down to keep subj-sesh correspondence)
                        outDF(((s+(t*20)+60)-20),5)={numTRsVS};
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
                        outDF(((s+(t*20)+60)-20),6)={mean(FD)};	
		else
                        disp([m1FP ' not found'])
		end
	end
end
% save out matrix
writetable(table(outDF),['/oak/stanford/groups/leanew1/users/apines/data/lsd_DMNMagMerged.csv'])
% and subject ID column
writetable(table(SubjNameCol),['/oak/stanford/groups/leanew1/users/apines/data/lsd_magsMerged_subjOrder.csv'])
