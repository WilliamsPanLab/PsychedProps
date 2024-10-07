function Extract_AmygFC(subj,sesh)
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
dmnIndL=find(L59k_dmn>.3);
dmnIndR=find(R59k_dmn>.3);
% add starting value to DMNindR
dmnIndR=dmnIndR+scale_3.diminfo{1}.models{2}.start-1;
% add tasks
tasks = {'rs1', 'rs2','gambling','emotion','wm','conscious','gonogo','nonconscious'};
% initialize master table to save out
mastTab=cell(1,5);
mastTabIterator=0;
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
		% annoying way to deal with readtable malfunction on set path #1 
		addpath(genpath(Paths{1}))
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
		% annoying way to deal with readtable malfunction on set path #2 
		restoredefaultpath
		conf=readtable(confpath,"FileType","text",'Delimiter', '\t');
		FD=table2array(conf(:,'framewise_displacement'));

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
	
		% get average correlation between DMN regions and amygdalae (DMNL - AL, DMNL - AR, DMNR - AL, DMNR - AR)
		DMNL_AL=corr(A_ts_L',DMN_mts_L');
		DMNL_AR=corr(A_ts_R',DMN_mts_L');
		DMNR_AL=corr(A_ts_L',DMN_mts_R');
	        DMNR_AR=corr(A_ts_R',DMN_mts_R');	
		% save out DMN-amyg FC (note FD not saved out in this iter, stored in other csv)
		T=table(mean([DMNL_AL' DMNL_AR' DMNR_AL' DMNR_AR']),'RowNames',"Row1");
		outFP=['/scratch/users/apines/data/mdma/' subj '/' sesh];
		writetable(T,[outFP '/' subj '_' sesh '_' task '_DMN_AmygFC_gro.csv'],'WriteRowNames',true)
	% if file doesnt exist
	else
	end
end
