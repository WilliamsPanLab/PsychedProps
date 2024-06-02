% extract network props from each session

% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% list current subjects
subjects={'PS03','PS16','PS18','PS19','PS21','PS24'};

% Pre data
pre_data = {
    {'Baseline1', 'Baseline2', 'Baseline3', 'Baseline4', 'Baseline5'},
    {'Baseline1', 'Baseline2', 'Baseline3', 'Baseline4', 'Baseline5', 'Baseline6', 'Baseline8'},
    {'Baseline1', 'Baseline2', 'Baseline3', 'Baseline4', 'Baseline5'},
    {'Baseline1', 'Baseline3', 'Baseline4', 'Baseline5'},
    {'Baseline1', 'Baseline2', 'Baseline3', 'Baseline4'},
    {'Baseline1', 'Baseline2', 'Baseline3'}
};

% Psilocybin data
psil_data = {
    {'Drug2'},
    {'Drug2'},
    {'Drug2'},
    {'Drug1'},
    {'Drug2'},
    {'Drug1'}
};

% Methylphenidate data
meth_data = {
    {'Drug1'},
    {'Drug1'},
    {'Drug1'},
    {'Drug2'},
    {'Drug1'},
    {'Drug2'}
};

% including follow up psil scans for these folks
additional_data = {
    {'PS93'}, % For ps03
    {'PS96'}, % For ps16
    {'PS98'}, % For ps18
    {}, % For ps19
    {}, % For ps21
    {} % For ps24
};

% Initialize the struct
preStruct = struct();

% Populate the struct
for i = 1:length(subjects)
    subject = subjects{i};
    preStruct.(subject).pre = pre_data{i};
    preStruct.(subject).psil = psil_data{i};
    preStruct.(subject).meth = meth_data{i};
    if ~isempty(additional_data{i})
        preStruct.(subject).additional = additional_data{i};
    end
end

% initialize group-level dfs 
outDF_pre=zeros(10,10,6);
outDF_psil=zeros(10,10,6);
outDF_meth=zeros(10,10,6);

% set base fp for loading
basefp='/scratch/users/apines/data/psil/';

% and specific rs scans for looping
tasks={'rs1' 'rs2' 'rs3' 'rs4' 'rs5' 'rs6'};
% for collecting each task (up to resting state 6: will need if exist thingy)
for t=1:length(tasks);
	task=tasks{t};
	for s=1:length(subjects)
		subject=subjects{s};
		% return this subjects specific scans
		preScans=preStruct.(subject).pre;
		psilScans=preStruct.(subject).psil;
		methScans=preStruct.(subject).meth;
		% for each "pre" scan
		for p=1:length(preScans);
			preScan=preScans{p};
			% load in da hist
			fp=[basefp subject '/' preScan '/' subject '_' preScan '_task-' task '_DMN_2dhist.csv'];
			if exist(fp,'file')
				preScan2dhist=readmatrix(fp);
				outDF_pre(:,:,s)=outDF_pre(:,:,s)+preScan2dhist;
			end
		end
		% for each psil scan
                for p=1:length(psilScans);
                        psilScan=psilScans{p};
                        % load in da hist
                        fp=[basefp subject '/' psilScan '/' subject '_' psilScan '_task-' task '_DMN_2dhist.csv'];
                        if exist(fp,'file')
                                psilScan2dhist=readmatrix(fp);
                                outDF_psil(:,:,s)=outDF_psil(:,:,s)+psilScan2dhist;
                        end
                end
		% for each meth scan
                for p=1:length(methScans);
                        methScan=methScans{p};
                        % load in da hist
                        fp=[basefp subject '/' methScan '/' subject '_' methScan '_task-' task '_DMN_2dhist.csv'];
                        if exist(fp,'file')
                                methScan2dhist=readmatrix(fp);
                                outDF_meth(:,:,s)=outDF_meth(:,:,s)+methScan2dhist;
                        end
                end
		% if subject has extra psilscans
		if isfield(preStruct.(subject), 'additional_data')
			surrogateSubjName=preStruct.(subject).additional_data;
			% load in drug1 to add to psil sessions
			fp=[basefp surrogateSubjName '/Drug1/' surrogateSubjName '_Drug1_task-' task '_DMN_2dhist.csv'];
			if exist(fp,'file')
                                psilScan2dhist=readmatrix(fp);
                                outDF_psil(:,:,s)=outDF_psil(:,:,s)+psilScan2dhist;
                        end
		end
	end
end
% save out matrices
save('/oak/stanford/groups/leanew1/users/apines/data/pre_psil_DMN_2dhist_merged.mat','outDF_pre');
save('/oak/stanford/groups/leanew1/users/apines/data/during_psil_DMN_2dhist_merged.mat','outDF_psil');
save('/oak/stanford/groups/leanew1/users/apines/data/during_meth_DMN_2dhist_merged.mat','outDF_meth');
