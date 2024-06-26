%restoredefaultpath

% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% get subj list
subjList=readtable('~/connectomeSubjs.txt','ReadVariableNames',false)
numSubjs=size(subjList,1);
% load in FD,remTRs from an R export?
% and a binary include from main analyses
% no, will do this in the rscript (vertexwise stat testing

% add TSNR mask, includes medial wall
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
mw_L=zeros(1,2562);
mw_L(mwAndTSNR_L>0)=1;
mw_R=zeros(1,2562);
mw_R(mwAndTSNR_R>0)=1;

% for each vertex (MAINTAIN OG VERTEX NUMBER 1:2562 FOR CLARITY)
for v=1:2562
	disp(v)
	% initialize vertex table - left
	vertTableL = [];
	% initialize vertex table - right
	vertTableR = [];
	% initialize a vector for subject
	subjsVec={};
	% sesssion
	seshsVec={};
	% task
	tasksVec={};
	% initialize vertex output, Left and Right
	vertValsL=[];
	vertValsR=[];
	% for each task
	for tasks=["rs1" "rs2" "wm"]
		task=char(tasks);
		% for each DES subject
		for s=1:numSubjs
			subj=char(subjList{s,1});
			inFP=['/scratch/users/apines/data/mdma/' subj '/'];
		 	ses01fp=[inFP '/ses-BL/' subj '_ses-BL_' task '_k1_AngDistMat_task.mat'];
		 	if exist(ses01fp,'file')
				subjVals01=load(ses01fp);
				% extract this vertex for each timepoint - left
				valuesL01=subjVals01.AngDist.Left(v,:);
				% extract this vertex for each timepoint - right
				valuesR01=subjVals01.AngDist.Right(v,:);
				% append vertex table with every timepoint value for scan1, with repmat of subj sesh and task as 1st 3 columns	
				vertTableL = [vertTableL; repmat({subj, 'ses-BL', task}, size(valuesL01, 2), 1), num2cell(valuesL01)'];
				vertTableR = [vertTableR; repmat({subj, 'ses-BL', task}, size(valuesR01, 2), 1), num2cell(valuesR01)'];
			else
			end
		% end for each subject
		end
	% end for each task
	end
	% if this vertex isnt in noMW_orlowTSNR_L
	if mw_L(v)==0;
		% saveout long df for this vertex, left
		vertTableL = cell2table(vertTableL, 'VariableNames', {'Subject', 'Session', 'Task', 'Value'});
		filename=['/scratch/users/apines/taskVerts/v' num2str(v) '_DES_L.csv'];
		writetable(vertTableL,filename);
	else
	end
	% if this vertex isnt in noMW_orlowTSNR_R
	if mw_R(v)==0;
		% saveout long df for this vertex, left
                vertTableR = cell2table(vertTableR, 'VariableNames', {'Subject', 'Session', 'Task', 'Value'});
		filename=['/scratch/users/apines/taskVerts/v' num2str(v) '_DES_R.csv'];
                writetable(vertTableR,filename);
        else
        end
% end for each vertex
end
% THE GOAL IS A LONG DATAFRAME FOR EACH VERTEX FOR EACH HEMISPHERE. 
% COLUMNS AS SUBJECT, FRAMEWISE ANGULAR DISTANCE FROM INTONETWORKANGLE, TASK, MEAN FD, REMAINING TRS, DRUG VS. NODRUG
