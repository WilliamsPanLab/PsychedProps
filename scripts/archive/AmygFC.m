function AmygFC
%establish normative FC between Nac Shell,Lat Amyg, MD thal, and DMN in this dataset
scale3fp='/oak/stanford/groups/leanew1/users/apines/maps/Tian_Subcortex_S3_3T_32k.dlabel.nii';
% load in concatenated cifti, calculate FC for this subject, print out onto two maps (L and R)
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(Paths{1}))

% load in atlas
scale_3=read_cifti(scale3fp);

% load in scale 3 data
% 9 for right DM thal
DMT_rInd=find(scale_3.cdata==9);
% 34 for left DM thal
DMT_lInd=find(scale_3.cdata==34);
% 22 for nac shell right hemi
NAC_rInd=find(scale_3.cdata==22);
%% 47 for nac shell left hemi
NAC_lInd=find(scale_3.cdata==47);
% rh lateral amyg - 19
AmyLat_rInd=find(scale_3.cdata==19);
% lh lateral amyg 44
AmyLat_lInd=find(scale_3.cdata==44);

% get default mode indices from NMF
dmnFile='/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/group_Nets_fullRes.mat';
dmn=load(dmnFile).nets;
% DMN is second column
dmnL=dmn.Lnets(:,2);
dmnR=dmn.Rnets(:,2);
%%% get indices that correspond to left hemisphere out of the 32492
% this is a multi-step procedure, god knows why. First get true 32492 vetors
L32k=zeros(1,32492);
R32k=zeros(1,32492);
% add values of 1 where not medial wall using vertlist (+1 because vertlist starts at 0)
L32k(scale_3.diminfo{1}.models{1}.vertlist+1)=1;
R32k(scale_3.diminfo{1}.models{2}.vertlist+1)=1;
% get medial-wall masked version of DMN values for 59k total cortical vertices (for matching to 91k cifti)
L59k_dmn=dmnL(logical(L32k));
R59k_dmn=dmnR(logical(R32k));
% get dmn indices
dmnIndL=find(L59k_dmn>.4);
dmnIndR=find(R59k_dmn>.4);
% add starting value to DMNindR
dmnIndR=dmnIndR+scale_3.diminfo{1}.models{2}.start-1;

% test printout cifti
templateCifti=read_cifti('/oak/stanford/groups/leanew1/users/apines/NeuroSynthMaps/dmn_smooth.dscalar.nii')
templateCifti.cdata=zeros(91282,1);
templateCifti.cdata(dmnIndL)=1;
templateCifti.cdata(dmnIndR)=1;
templateCifti.cdata(DMT_lInd)=2;
templateCifti.cdata(DMT_rInd)=2;
templateCifti.cdata(NAC_lInd)=3;
templateCifti.cdata(NAC_rInd)=3;
templateCifti.cdata(AmyLat_lInd)=4;
templateCifti.cdata(AmyLat_rInd)=4;
write_cifti(templateCifti,'/oak/stanford/groups/leanew1/users/apines/maps/testAtlas.dscalar.nii')
% test checks out

% loop over each subject
% read in subject dosage correspondence, has to be before addpath for some silly reason
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');
tasks = {'rs1', 'rs2','gambling','emotion','wm','conscious','gonogo','nonconscious'};
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
			% get proper dtseries
			matchingFiles = {};
			% sep naming convention for rs
			if string(task)=="rs1"
				for i = 1:length(subjfiles)
                                	if contains(subjfiles(i).name, 'rs_acq-mb_dir-pe0') && contains(subjfiles(i).name, 'desc-denoisedSmoothed_bold.dtseries.nii')
                                	        matchingFiles = [matchingFiles, subjfiles(i).name];
                                	end
                   		end
			elseif string(task)=="rs2"
                                for i = 1:length(subjfiles)
                                	if contains(subjfiles(i).name, 'rs_acq-mb_dir-pe1') && contains(subjfiles(i).name, 'desc-denoisedSmoothed_bold.dtseries.nii')
                                        	matchingFiles = [matchingFiles, subjfiles(i).name];
                                	end
				end
                        else	
				for i = 1:length(subjfiles)
					if contains(subjfiles(i).name, task) && contains(subjfiles(i).name, 'desc-denoisedSmoothed_bold.dtseries.nii')
                				matchingFiles = [matchingFiles, subjfiles(i).name];
            				end
        			end
			end	
			% if file exists
			if isempty(matchingFiles)==0
				mastTabIterator=mastTabIterator+1		
				% reconstruct cifti path
				ciftipath=[subjdir matchingFiles{1}];	
				% load in cifti
				concatData=read_cifti(ciftipath);
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
				% extract DM thalamic TS - left and right
				% left
				DMT_ts_L=concatData.cdata(DMT_lInd,:);
				DMT_mts_L=mean(DMT_ts_L);
				% right
				DMT_ts_R=concatData.cdata(DMT_rInd,:);
				DMT_mts_R=mean(DMT_ts_R);

				% extract nac shell TS - left and right
				% left
				NAC_ts_L=concatData.cdata(NAC_lInd,:);
				NAC_mts_L=mean(NAC_ts_L);
				% right
				NAC_ts_R=concatData.cdata(NAC_rInd,:);
				NAC_mts_R=mean(NAC_ts_R);

				% extract amygdalar TS - left and right
				% left
				A_ts_L=concatData.cdata(AmyLat_lInd,:);
				A_mts_L=mean(A_ts_L);
				% right
				A_ts_R=concatData.cdata(AmyLat_rInd,:);
				A_mts_R=mean(A_ts_R);
	
				% DMN TS - left and right
				% left
				DMN_ts_L=concatData.cdata(dmnIndL,:);
				DMN_mts_L=mean(DMN_ts_L);
				% right
				DMN_ts_R=concatData.cdata(dmnIndR,:);
				DMN_mts_R=mean(DMN_ts_R);
	
				% get big correlation matrix between all derived mean-time series
				% initialize matrix
				corrMat=zeros(8,8);
				% fill in
				corrMat(1,2)=corr(DMT_mts_L',DMT_mts_R');
				corrMat(1,3)=corr(DMT_mts_L',NAC_mts_L');
				corrMat(1,4)=corr(DMT_mts_L',NAC_mts_R');
				corrMat(1,5)=corr(DMT_mts_L',A_mts_L');
				corrMat(1,6)=corr(DMT_mts_L',A_mts_R');
				corrMat(1,7)=corr(DMT_mts_L',DMN_mts_L');
				corrMat(1,8)=corr(DMT_mts_L',DMN_mts_R');
				corrMat(2,3)=corr(DMT_mts_R',NAC_mts_L');
				corrMat(2,4)=corr(DMT_mts_R',NAC_mts_R');
				corrMat(2,5)=corr(DMT_mts_R',A_mts_L');
				corrMat(2,6)=corr(DMT_mts_R',A_mts_R');
				corrMat(2,7)=corr(DMT_mts_R',DMN_mts_L');
				corrMat(2,8)=corr(DMT_mts_R',DMN_mts_R');
				corrMat(3,4)=corr(NAC_mts_L',NAC_mts_R');
				corrMat(3,5)=corr(NAC_mts_L',A_mts_L');
				corrMat(3,6)=corr(NAC_mts_L',A_mts_R');
				corrMat(3,7)=corr(NAC_mts_L',DMN_mts_L');
				corrMat(3,8)=corr(NAC_mts_L',DMN_mts_R');
				corrMat(4,5)=corr(NAC_mts_R',A_mts_L');
				corrMat(4,6)=corr(NAC_mts_R',A_mts_R');
				corrMat(4,7)=corr(NAC_mts_R',DMN_mts_L');
				corrMat(4,8)=corr(NAC_mts_R',DMN_mts_R');
				corrMat(5,6)=corr(A_mts_L',A_mts_R');
				corrMat(5,7)=corr(A_mts_L',DMN_mts_L');
				corrMat(5,8)=corr(A_mts_L',DMN_mts_R');
				corrMat(6,7)=corr(A_mts_R',DMN_mts_L');
				corrMat(6,8)=corr(A_mts_R',DMN_mts_R');
				corrMat(7,8)=corr(DMN_mts_L',DMN_mts_R');
				% fill in lower triangle
				corrMat=corrMat+corrMat'-diag(diag(corrMat));
			        % extract unitary AmFug FC
				UniAmFugFC=mean([corrMat(1,2) corrMat(1,3) corrMat(1,4) corrMat(1,5) corrMat(1,6) corrMat(2,3) corrMat(2,4) corrMat(2,5) corrMat(2,6) corrMat(3,4) corrMat(3,5) corrMat(3,6) corrMat(4,5) corrMat(4,6) corrMat(5,6)]);
				% populate unitary AmFug FC
				mastTab{mastTabIterator,1}=subj;
				mastTab{mastTabIterator,2}=seshIterator;		
				mastTab{mastTabIterator,3}=task;
				mastTab{mastTabIterator,4}=mean(FD(~isnan(FD)));
				mastTab{mastTabIterator,5}=UniAmFugFC	
			% if file doesnt exist
			else
			end
		end
	end
end
% create column names
mastTabColNames = 'subj,sesh,task,FD,AmFugFC';

% Split the string into a cell array using the ',' delimiter
columnNames = strsplit(mastTabColNames, ',');

outTable=cell2table(mastTab,'VariableNames',columnNames)
writetable(outTable,'~/AmFugFC.csv')
% save out a long-format subj, dose, task dataframe from uniary AmFug FC structure
% get average correlation matrix
meanCor=mean(corrMats,3);
