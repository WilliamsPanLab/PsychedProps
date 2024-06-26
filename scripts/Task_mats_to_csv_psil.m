function Task_mats_to_csv(subj,sesh,task)
%restoredefaultpath

% read in subj-session-dose correspondence
%subSeshDose=readtable('~/subjSeshDoseCorresp.csv');

% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% get subj list
% get subj list
subjPrefix=repmat('PS',11,1);
subjSuffix=["02","03","16","18","19","21","24","93","96","98","99"];
subjList=strcat(subjPrefix,subjSuffix')

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

	% for each "task"
	for tasks=["rs1" "rs2" "rs3" "rs4" "rs5" "rs6"]
		task=char(tasks);
		% for each subj except 2
        	for s=[2 3 4 5 6 7 8 9 10 11];
                	disp(s)
			subj=cellstr(subjList(s));	
			inFP=['/scratch/users/apines/data/psil/' subj '/'];
			% for up to 8 iterations
			for i=1:8
				subj=subjList(s);	
		 		% baseline
				bvFP=[inFP '/Baseline' num2str(i) '/' subj '_Baseline' num2str(i) '_' task '_k1_AngDistMat_task.mat'];
				bvFP=strjoin(bvFP,'');
				if exist(bvFP,'file')
					subjValsBV=load(bvFP);
                                	% extract this vertex for each timepoint - left
                                	valuesLBV=subjValsBV.AngDist.Left(v,:);
                                	% extract this vertex for each timepoint - right
                                	valuesRBV=subjValsBV.AngDist.Right(v,:);
                                	% append vertex table with every timepoint value for scan1, with repmat of subj sesh and task as 1st 3 columns
                                	vertTableL = [vertTableL; repmat({subj, 'Baseline', task}, size(valuesLBV, 2), 1), num2cell(valuesLBV)'];
                                	vertTableR = [vertTableR; repmat({subj, 'Baseline', task}, size(valuesRBV, 2), 1), num2cell(valuesRBV)'];
				else
				end
				% between
                                bwFP=[inFP '/Between' num2str(i) '/' subj '_Between' num2str(i) '_' task '_k1_AngDistMat_task.mat'];
                                bwFP=strjoin(bvFP,'');
                                if exist(bwFP,'file')
                                        subjValsBW=load(bwFP);
                                        % extract this vertex for each timepoint - left
                                        valuesLBW=subjValsBW.AngDist.Left(v,:);
                                        % extract this vertex for each timepoint - right
                                        valuesRBW=subjValsBW.AngDist.Right(v,:);
                                        % append vertex table with every timepoint value for scan1, with repmat of subj sesh and task as 1st 3 columns
                                        vertTableL = [vertTableL; repmat({subj, 'Between', task}, size(valuesLBW, 2), 1), num2cell(valuesLBW)'];
                                        vertTableR = [vertTableR; repmat({subj, 'Between', task}, size(valuesRBW, 2), 1), num2cell(valuesRBW)'];
                                else
                                end
				% after
                                aFP=[inFP '/After' num2str(i) '/' subj '_After' num2str(i) '_' task '_k1_AngDistMat_task.mat'];
                                aFP=strjoin(aFP,'');
                                if exist(aFP,'file')
                                        subjValsa=load(aFP);
                                        % extract this vertex for each timepoint - left
                                        valuesLa=subjValsa.AngDist.Left(v,:);
                                        % extract this vertex for each timepoint - right
                                        valuesRa=subjValsa.AngDist.Right(v,:);
                                        % append vertex table with every timepoint value for scan1, with repmat of subj sesh and task as 1st 3 columns
                                        vertTableL = [vertTableL; repmat({subj, 'After', task}, size(valuesLa, 2), 1), num2cell(valuesLa)'];
                                        vertTableR = [vertTableR; repmat({subj, 'After', task}, size(valuesRa, 2), 1), num2cell(valuesRa)'];
                                else
                                end
				% drug
                                dFP=[inFP '/Drug' num2str(i) '/' subj '_Drug' task '_k1_AngDistMat_task.mat'];
                                dFP=strjoin(dFP,'');
                                if exist(dFP,'file')
                                        subjValsd=load(dFP);
                                        % extract this vertex for each timepoint - left
                                        valuesLd=subjValsd.AngDist.Left(v,:);
                                        % extract this vertex for each timepoint - right
                                        valuesRd=subjValsd.AngDist.Right(v,:);
					% NEED TO RECORD DRUG 1 OR 2 FOR DECODING LATER WITH R OUTPUT
					seshstring=['Drug' num2str(i)];
					seshstring=strjoin(seshstring,'');
                                        % append vertex table with every timepoint value for scan1, with repmat of subj sesh and task as 1st 3 columns
                                        vertTableL = [vertTableL; repmat({subj, seshstring, task}, size(valuesLd, 2), 1), num2cell(valuesLd)'];
                                        vertTableR = [vertTableR; repmat({subj, seshstring, task}, size(valuesRd, 2), 1), num2cell(valuesRd)'];
                                else
                                end
			% end for each iteration
			end
		% end for each subject
		end
	% end for each task
	end
	% if this vertex isnt in noMW_orlowTSNR_L
	if mw_L(v)==0;
		% saveout long df for this vertex, left
		vertTableL = cell2table(vertTableL, 'VariableNames', {'Subject', 'Session', 'Task', 'Value'});
		filename=['/scratch/users/apines/taskVerts/v' num2str(v) '_psil_L.csv'];
		writetable(vertTableL,filename);
	else
	end
	% if this vertex isnt in noMW_orlowTSNR_R
	if mw_R(v)==0;
		% saveout long df for this vertex, left
                vertTableR = cell2table(vertTableR, 'VariableNames', {'Subject', 'Session', 'Task', 'Value'});
		filename=['/scratch/users/apines/taskVerts/v' num2str(v) '_psil_R.csv'];
                writetable(vertTableR,filename);
        else
        end
% end for each vertex
end
% THE GOAL IS A LONG DATAFRAME FOR EACH VERTEX FOR EACH HEMISPHERE. 
% COLUMNS AS SUBJECT, FRAMEWISE ANGULAR DISTANCE FROM INTONETWORKANGLE, TASK, MEAN FD, REMAINING TRS, DRUG VS. NODRUG
