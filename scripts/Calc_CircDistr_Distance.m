% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% get subjseshdose correspondence
% read in subj-session-dose correspondence
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');
% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')
% for each subject
for s=[1 2 3 5 7 8 9 11 12 13 14 15 16 17];
	% display s
	disp(s)
	% get session info
	seshInfo=subSeshDose{s,2:5};
	
	% write each session out instead of for loop for sesh's
		% load in disributions
		bvfpl=[parentFP '/' subj '_' seshInfo{1} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
		bvfpr=[parentFP '/' subj '_' seshInfo{1} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
		bvopFl_L=load(bvfpl);
		bvopFl_R=load(bvfpr);
		% placebo
		pfpl=[parentFP '/' subj '_' seshInfo{2} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
                pfpr=[parentFP '/' subj '_' seshInfo{2} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
                popFl_L=load(pfpl);
                popFl_R=load(pfpr);
		% 80 mg
		m1fpl=[parentFP '/' subj '_' seshInfo{3} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
                m1fpr=[parentFP '/' subj '_' seshInfo{3} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
                m1pFl_L=load(m1fpl);
                m1pFl_R=load(m1fpr);
                % 120 mg
		m2fpl=[parentFP '/' subj '_' seshInfo{4} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
                m2fpr=[parentFP '/' subj '_' seshInfo{4} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
                m2pFl_L=load(m2fpl);
                m2pFl_R=load(m2fpr);



		% combine distributions across resting states
		
		% store in arrays
		% end for each sesh
	% kuiper test each vertex
	% pl vs. 80 
	% pl vs 120
	% pl. vs drug
	% bv vs. drug
	% bv vs. plac
	% vert surface print .pngs
% end for each subject
% average k stat over each vertex
% print comparisons
% pl vs. 80 
% pl vs 120
% pl. vs drug
% bv vs. drug
% bv vs. plac

