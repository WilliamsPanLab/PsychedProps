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

% INITIALIZE HISTC TOTAL TO PRINT OUT AT END, RUN AGAIN AFTER 25TH 50TH AND 75TH PERCENTILE KNOWN
histCountTot=zeros(1,50);

% SAVE OUT 25TH 50TH AND 75TH PERCENTILE AFTER 

% set common fp
commonFP=['/scratch/users/apines/data/mdma/'];

% for collecting each task
for task=["rs1" "rs2" "emotion" "gambling" "wm"]
	% initialize output table: difference in alff by each network, for 80 and 120, and in group and individualized atlases
	outDF=zeros(17,20);
	% for each subj except 4 6 10
	for s=[1 2 3 5 7 8 9 11 12 13 14 15 16 17]
		disp(s)
	        % get session info NOTE K4 HAS ALL 4 VALUES, K1 has INTENDED VALUE AS WELL
	        seshInfo=subSeshDose{s,2:5};
	        bvFP=[commonFP subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_' task '_BU_counts.mat'];
	        bvFP=strjoin(bvFP,'');
		% and grab reverse phase encoding direction nifti
		pFP=[commonFP subjList(s) '/' seshInfo{2} '/' subjList(s) '_' seshInfo{2} '_' task '_BU_counts.mat'];
	        pFP=strjoin(pFP,'');
	        m1FP=[commonFP subjList(s) '/' seshInfo{3} '/' subjList(s) '_' seshInfo{3} '_' task '_BU_counts.mat'];
	        m1FP=strjoin(m1FP,'');
	        m2FP=[commonFP subjList(s) '/' seshInfo{4} '/' subjList(s) '_' seshInfo{4} '_' task '_BU_counts.mat'];
	        m2FP=strjoin(m2FP,'');
		% load in baseline mats
		if exist(bvFP,'file')
		bv=load(bvFP);
	        outDF(s,1)=bv.BUPcount;
		%%%%
		bvFP2=[commonFP subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_' task '_BU_dur_counts.mat'];
		bvFP2=strjoin(bvFP2,'');
		bv2=load(bvFP2);
		histCountTot=histCountTot+bv2.BU_dur_counts;
		%%%%
		totsum=sum(bv2.BU_dur_counts);
		outDF(s,2)=sum(bv2.BU_dur_counts(4:end))/totsum;
		outDF(s,3)=sum(bv2.BU_dur_counts(6:end))/totsum;
		outDF(s,4)=sum(bv2.BU_dur_counts(10:end))/totsum;
		%%%% MEAN DURATION
		bvFP3=[commonFP subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_' task '_BU_meanDur.mat'];
                bvFP3=strjoin(bvFP3,'');
                bv3=load(bvFP3);
		outDF(s,3)=bv3.meanBUPdur;
		% load in remaining frames
	        childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{1} ];
                bvCSIfp=[childfp '/' subjList(s) '_' seshInfo{1} '_task-' task '_ValidSegments_Trunc.txt'];
                bvCSIfp=strjoin(bvCSIfp,'');
                CSI = importdata(bvCSIfp);
                % trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
                numTRsVS=sum(CSI(:,2));
                outDF(s,17)=numTRsVS;	
		end
		% load in placebo mats
		if exist(pFP,'file')
		p=load(pFP);
		% extract in group consensus atlas: placebo
		outDF(s,5)=p.BUPcount;
		%%%%
                pFP2=[commonFP subjList(s) '/' seshInfo{2} '/' subjList(s) '_' seshInfo{2} '_' task '_BU_dur_counts.mat'];
                pFP2=strjoin(bvFP2,'');
                p2=load(pFP2);
                histCountTot=histCountTot+p2.BU_dur_counts;
                %%%%
		totsum=sum(p2.BU_dur_counts);
		outDF(s,6)=sum(p2.BU_dur_counts(4:end))/totsum;
		outDF(s,7)=sum(p2.BU_dur_counts(6:end))/totsum;
		outDF(s,8)=sum(p2.BU_dur_counts(10:end))/totsum;;
		%%%% MEAN DURATION
                pFP3=[commonFP subjList(s) '/' seshInfo{2} '/' subjList(s) '_' seshInfo{2} '_' task '_BU_meanDur.mat'];
                pFP3=strjoin(pFP3,'');
                p3=load(pFP3);
                outDF(s,7)=p3.meanBUPdur;
		% remaining TRs
		childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{2} ];
                pCSIfp=[childfp '/' subjList(s) '_' seshInfo{2} '_task-' task '_ValidSegments_Trunc.txt'];
                pCSIfp=strjoin(pCSIfp,'');
                CSI = importdata(pCSIfp);
                numTRsVS=sum(CSI(:,2));
                outDF(s,18)=numTRsVS;
		end
		% load mdma 1 (80 mg)
		if exist(m1FP,'file')
		m1=load(m1FP);
		% extact in group consensus atlas: mdma 1
		outDF(s,9)=m1.BUPcount;
		%%%%
                m1FP2=[commonFP subjList(s) '/' seshInfo{3} '/' subjList(s) '_' seshInfo{3} '_' task '_BU_dur_counts.mat'];
                m1FP2=strjoin(m1FP2,'');
                m1_2=load(m1FP2);
                histCountTot=histCountTot+m1_2.BU_dur_counts;
                %%%%
		totsum=sum(m1_2.BU_dur_counts);
		outDF(s,10)=sum(m1_2.BU_dur_counts(4:end))/totsum;
		outDF(s,11)=sum(m1_2.BU_dur_counts(6:end))/totsum;
		outDF(s,12)=sum(m1_2.BU_dur_counts(10:end))/totsum;
		%%%% MEAN DURATION
                m1FP3=[commonFP subjList(s) '/' seshInfo{3} '/' subjList(s) '_' seshInfo{3} '_' task '_BU_meanDur.mat'];
                m1FP3=strjoin(m1FP3,'');
                m13=load(m1FP3);
                outDF(s,11)=m13.meanBUPdur;
		% remaining TRs
		childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{3} ];
                m1CSIfp=[childfp '/' subjList(s) '_' seshInfo{3} '_task-' task '_ValidSegments_Trunc.txt'];
                m1CSIfp=strjoin(m1CSIfp,'');
                CSI = importdata(m1CSIfp);
                numTRsVS=sum(CSI(:,2));
                outDF(s,19)=numTRsVS;
		end
		% load mdma 2 (120 mg)
		if exist(m2FP,'file')
		m2=load(m2FP);
		% extract in group consensus atlas: mdma 2
		outDF(s,13)=m2.BUPcount;
		%%%%
                m2FP2=[commonFP subjList(s) '/' seshInfo{4} '/' subjList(s) '_' seshInfo{4} '_' task '_BU_dur_counts.mat'];
                m2FP2=strjoin(m2FP2,'');
                m2_2=load(m2FP2);
                histCountTot=histCountTot+m2_2.BU_dur_counts;
                %%%%
		totsum=sum(m2_2.BU_dur_counts);
		outDF(s,14)=sum(m2_2.BU_dur_counts(4:end))/totsum;
		outDF(s,15)=sum(m2_2.BU_dur_counts(6:end))/totsum;
		outDF(s,16)=sum(m2_2.BU_dur_counts(10:end))/totsum;
		%%%% MEAN DURATION
                m2FP3=[commonFP subjList(s) '/' seshInfo{4} '/' subjList(s) '_' seshInfo{4} '_' task '_BU_meanDur.mat'];
                m2FP3=strjoin(m2FP3,'');
                m23=load(m2FP3);
                outDF(s,15)=m23.meanBUPdur;
		% remaining TRs
		childfp=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{4} ];
                m2CSIfp=[childfp '/' subjList(s) '_' seshInfo{4} '_task-' task '_ValidSegments_Trunc.txt'];
                m2CSIfp=strjoin(m2CSIfp,'');
                CSI = importdata(m2CSIfp);
                numTRsVS=sum(CSI(:,2));
                outDF(s,20)=numTRsVS;
		end
	end
	% save out matrix
	writematrix(outDF,strjoin(['/oak/stanford/groups/leanew1/users/apines/data/' task '_BUP_freq_merged.csv'],''))
end

%%% for normative dwell time distributions
% percentiles to find
percentiles = [0.25, 0.50, 0.75];
% initialize array to store the bin values corresponding to percentiles
percentileBins = zeros(size(percentiles));
totalCount=sum(histCountTot);
% cumulative sum for percentiles
cumSumHist = cumsum(histCountTot);
% find the bins corresponding to the specified percentiles
for i = 1:length(percentiles)
    % find the index where the cumulative sum first exceeds the percentile threshold
    index = find(cumSumHist >= (percentiles(i) * totalCount), 1);
    % store the bin value corresponding to this index
    percentileBins(i) = index;
end
% looks like it's bins 3 5 and 9: will save out proportion over these percentiles as session-level measurement
