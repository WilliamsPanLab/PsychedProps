function AmygFC
%establish normative FC between Nac Shell,Lat Amyg, MD thal, and DMN in this dataset
restoredefaultpath

% loop over each subject
% read in subject dosage correspondence, has to be before addpath for some silly reason
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');
tasks = {'rs1', 'rs2','gambling','emotion','wm'};
% combine angular distance from reference streams across subjects
subjects = {'sub-MDMA001', 'sub-MDMA002', 'sub-MDMA003', 'sub-MDMA005', 'sub-MDMA007', 'sub-MDMA008', 'sub-MDMA009', 'sub-MDMA011', 'sub-MDMA012', 'sub-MDMA013', 'sub-MDMA014', 'sub-MDMA015', 'sub-MDMA016', 'sub-MDMA017'};

subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')
% start iterator
iterator=0
corrMats=zeros(8,8,length([1 2 3 5 7 9 11 12 13 14 15 16 17]));
% initialize master table to save out
mastTab=cell(1,5);
mastTabIterator=0;
% aggregate each subject, session, and task for each streams
for s=[1 2 3 5 7 9 11 12 13 14 15 16 17]
        % add 1 to iterator
        iterator=iterator+1
        % get subject name
        subj=subjList{s}
	% get session info
        seshInfo=subSeshDose{s,2:5};
        % for placebo, 80mg, and 120mg
        seshArray={seshInfo{2} seshInfo{3} seshInfo{4}};
	% make an iterator for session
        seshIterator=0;
	for sessioncell=seshArray
		sesh=sessioncell{1}
		% update session iterator
		seshIterator=seshIterator+1;		
		% set cifti path
		subjdir = ['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' subj  '/' sesh  '/func/'];
		% get list of files
		subjfiles=dir(fullfile(subjdir));
		% get list of fmriprep output files
		fmriprepdir = ['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subj '/' sesh '/func/'];
		fmriprepfiles=dir(fullfile(fmriprepdir));
		% for each task
		for taskcell=tasks
			task=taskcell{1};
			% props data filepath
			parFP=['/scratch/users/apines/data/mdma/' subj '/' sesh];
			propsFP=[parFP '/' subj '_' sesh '_' task '_Prop_Feats_groNew.csv'];
			% if file exists
			if isfile(propsFP)==1
				props=readtable(propsFP);
				DMNProps=table2array(props(1,2));
				mastTabIterator=mastTabIterator+1		
				% get FD for this run 
				confFile={};
				% sep. naming convention for rs
				if string(task)=="rs1"
					for i = 1:length(fmriprepfiles)
						if contains(fmriprepfiles(i).name, 'rs_acq-mb_dir-pe0') && contains(fmriprepfiles(i).name, '_desc-confounds_timeseries.tsv')
                                                	confFile = [confFile, fmriprepfiles(i).name];
                                        	end
					end
				elseif string(task)=="rs2"
					for i = 1:length(fmriprepfiles)
						if contains(fmriprepfiles(i).name, 'rs_acq-mb_dir-pe1') && contains(fmriprepfiles(i).name, '_desc-confounds_timeseries.tsv')
                                        	        confFile = [confFile, fmriprepfiles(i).name];
                                       		end
					end
				else
					for i = 1:length(fmriprepfiles)
						if contains(fmriprepfiles(i).name, task) && contains(fmriprepfiles(i).name, '_desc-confounds_timeseries.tsv')
							confFile = [confFile, fmriprepfiles(i).name];
						end
					end
				end
				confpath = [fmriprepdir confFile{1}];
				conf=readtable(confpath,"FileType","text",'Delimiter', '\t');
				FD=table2array(conf(:,'framewise_displacement'));
				% populate unitary DMN props
				mastTab{mastTabIterator,1}=subj;
				mastTab{mastTabIterator,2}=seshIterator;		
				mastTab{mastTabIterator,3}=task;
				mastTab{mastTabIterator,4}=mean(FD(~isnan(FD)));
				mastTab{mastTabIterator,5}=DMNProps
			% if file doesnt exist
			else
			end
		end
	end
end
% create column names
mastTabColNames = 'subj,sesh,task,FD,DMNProps';

% Split the string into a cell array using the ',' delimiter
columnNames = strsplit(mastTabColNames, ',');

outTable=cell2table(mastTab,'VariableNames',columnNames)
writetable(outTable,'~/DMNProps.csv')
