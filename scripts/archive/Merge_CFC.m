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
% to match up with medial walls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use native freesurfer command for mw mask indices
surfML = '/oak/stanford/groups/leanew1/users/apines/surf/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/oak/stanford/groups/leanew1/users/apines/surf/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "is medial wall" vector for vertices
mw_L=zeros(1,2562);
mw_L(mwIndVec_l)=1;
mw_R=zeros(1,2562);
mw_R(mwIndVec_r)=1;
% convert to faces
% convert to faces
F_MW_L=sum(mw_L(faces_l),2)./3;
F_MW_R=sum(mw_R(faces_r),2)./3;
% convert "partial" medial wall to medial wall
F_MW_L=ceil(F_MW_L);
F_MW_R=ceil(F_MW_R);
% face mask indices
fmwIndVec_l=find(F_MW_L);
fmwIndVec_r=find(F_MW_R);
% make medial wall vector
g_noMW_combined_L=setdiff([1:5120],fmwIndVec_l);
g_noMW_combined_R=setdiff([1:5120],fmwIndVec_r);

% for each subject
% aggregate each subject, session, and task for each streams
for s=1:14
    subj=subjects{s}
    % for this session
    sesh=sessions{1};
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
	    	childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
		CSIfp=[childfp '/' subj '_' sesh '_task-' task '_ValidSegments_Trunc.txt'];
	    	CSI = importdata(CSIfp);
	    	% assure that TR count is the same between time series and valid segments txt
	    	SegNum=size(CSI);
	    	SegNum=SegNum(1);
	    	% trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
	    	numTRsVS=CSI(SegNum,1)+CSI(SegNum,2)-1;
	    	% multiply by number of TRs
	    	CFC=CFCL*numTRsVS;
	    	% add to aggregate CFC
	    	aggCFC(g_noMW_combined_L,g_noMW_combined_L) = aggCFC(g_noMW_combined_L,g_noMW_combined_L) + CFC;
		%aggCFC = aggCFC + CFC;
		% get sum of segments over 14
		CSIOver=CSI(CSI(:,2)>14,2);
		% add to total number of TRs
	    	totTRs = totTRs + sum(CSIOver)-1;
	end
    end
end
% divide by total number of TRs
aggCFC=aggCFC./totTRs;
% saveout
save('/scratch/users/apines/data/mdma/aggregateCFC.mat','aggCFC')
