function Mats_to_VertWisecsvs_Left
%restoredefaultpath

% read in subj-session-dose correspondence
%subSeshDose=readtable('~/subjSeshDoseCorresp.csv');

% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')

% for each vertex (MAINTAIN OG VERTEX NUMBER 1:2562 FOR CLARITY)
for v=1:1116
	disp(v)
	% initialize vertex table - left
	vertTableL = [];
	% initialize a vector for subject
	subjsVec={};
	% sesssion
	seshsVec={};
	% task
	tasksVec={};
	% initialize vertex output, Left and Right
	vertValsL=[];
	% for each task
	for tasks=["rs1" "rs2" "wm" "gambling"]
		task=char(tasks);
		% for each subject (not 4, 6, or 10)
		% CHANGING TO MEAN 6/28/24
		for s=[1 2 3 5 7 8 9 11 12 13 14 15 16 17]
			subj=subjList{s};
			inFP=['/scratch/users/apines/data/mdma/' subj '/'];
		 	ses01fp=[inFP '/ses-01/' subj '_ses-01_' task '_k1_AngDistMat.mat'];
		 	if exist(ses01fp,'file')
				subjVals01=load(ses01fp);
				% extract this vertex for each timepoint - left
				valuesL01=mean(subjVals01.AngDist.Left(v,:));
				% append vertex table with every timepoint value for scan1, with repmat of subj sesh and task as 1st 3 columns	
				vertTableL = [vertTableL; repmat({subj, 'ses-01', task}, size(valuesL01, 2), 1), num2cell(valuesL01)'];
			else
			end
			% session 2
			ses02fp=[inFP '/ses-02/' subj '_ses-02_' task '_k1_AngDistMat.mat'];
			if exist(ses02fp,'file')
				subjVals02=load(ses02fp);
		 		% extract this vertex for each timepoint - left
                                valuesL02=mean(subjVals02.AngDist.Left(v,:));
                                % append vertex table with every timepoint value for scan1, with repmat of subj sesh and task as 1st 3 columns
                                vertTableL = [vertTableL; repmat({subj, 'ses-02', task}, size(valuesL02, 2), 1), num2cell(valuesL02)'];
			else
			end
			% session 3
			ses03fp=[inFP '/ses-03/' subj '_ses-03_' task '_k1_AngDistMat.mat'];
			if exist(ses03fp,'file')
				subjVals03=load(ses03fp);
				% extract this vertex for each timepoint - left
                                valuesL03=mean(subjVals03.AngDist.Left(v,:));
                                % append vertex table with every timepoint value for scan1, with repmat of subj sesh and task as 1st 3 columns
                                vertTableL = [vertTableL; repmat({subj, 'ses-03', task}, size(valuesL03, 2), 1), num2cell(valuesL03)'];
                        else
                        end
		% end for each subject
		end
	% end for each task
	end
	% saveout long df for this vertex, left
	vertTableL = cell2table(vertTableL, 'VariableNames', {'Subject', 'Session', 'Task', 'Value'});
	filename=['/scratch/users/apines/taskVerts/v' num2str(v) '_L_DMN.csv'];
	writetable(vertTableL,filename);
% end for each vertex
end
% THE GOAL IS A LONG DATAFRAME FOR EACH VERTEX FOR EACH HEMISPHERE. 
% COLUMNS AS SUBJECT, FRAMEWISE ANGULAR DISTANCE FROM INTONETWORKANGLE, TASK, MEAN FD, REMAINING TRS, DRUG VS. NODRUG
