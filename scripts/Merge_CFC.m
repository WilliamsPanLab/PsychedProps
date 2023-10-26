% read in each circular FC matrix, weigh by number of TRs included, aggregate and saveout
% subject list
subjects = {'sub-MDMA001', 'sub-MDMA002', 'sub-MDMA003', 'sub-MDMA005', 'sub-MDMA007', 'sub-MDMA008', 'sub-MDMA009', 'sub-MDMA011', 'sub-MDMA012', 'sub-MDMA013', 'sub-MDMA014', 'sub-MDMA015', 'sub-MDMA016', 'sub-MDMA017'};
% only first session for now
sessions = {'ses-00'};
% task list
tasks = {'rs1', 'rs2','gambling','emotion','wm'};
% initialize aggregate CFC
aggCFC=zeros(5120);
% initialize total of number of TRs
totTRs=0;
% for each subject
% aggregate each subject, session, and task for each streams
for s=[1 2 3 5 7 8 9 11 12 13 14 15 16 17]
    subj=subjList{s}
    % for each task
    for taskcell=tasks
	    task=taskcell{1};
	    % get the file name
	    fn=['/scratch/users/apines/data/mdma/' subj '/' subj '_' sesh '_' task '_CircFC.mat'];
	    %if the file exists
	    if exist(fn, 'file') == 2
	    % load it in
	    CFC=load(fn);
	    % extract from struct
	    CFCL=CFC.adjMats.L;
	    % get number of TRs
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load in continuous segment indices
	    CSIfp=[childfp '/' subj '_' sesh '_task-' task '_ValidSegments_Trunc.txt'];
	    CSI = importdata(CSIfp);
	    % assure that TR count is the same between time series and valid segments txt
	    SegNum=size(CSI);
	    SegNum=SegNum(1);
	    % number of trs in fmri ts
	    of_ts_trs=length(L);
	    % trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
	    numTRsVS=CSI(SegNum,1)+CSI(SegNum,2)-1;
	    % note we lose a "between" volume for each segment
	    if numTRsVS ~= (of_ts_trs+SegNum)
        	disp('TRs from Valid Segments txt and derived vectors do not match. Fix it.')
	    else
	    end
	    % multiply by number of TRs
	    CFC=CFCL*SegNum;
	    % add to aggregate CFCifif
	    aggCFC = aggCFC + CFC;
	    % add to total number of TRs
	    totTRs = totTRs + length(CFC);
	end
    end
end
%   for each task
% if the file exists
% load it in
% multiply by number of TRs
% add to aggregate CFC
% end
% divide by total number of TRs

