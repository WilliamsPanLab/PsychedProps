function Calc_AvgBup_psil(subj)

% just need whether it's drug 1 or 2 that corresponds to psil
restoredefaultpath
subSeshDose=readtable('~/subjSeshDoseCorresp_psilo.csv');

% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/')
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% get spherical coordinates of each vertex
% Load in surface data
surfL = ['/oak/stanford/groups/leanew1/users/apines/surf/lh.sphere'];
surfR = ['/oak/stanford/groups/leanew1/users/apines/surf/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% faces_L
F_L=faces_l;
% vertices V
V_L=vx_l;
% faces_R
F_R=faces_r;
% vertices V
V_R=vx_r;

% load in medial wall + SNR so we don't have to loop over every single vertex and then mask out later
% add TSNR mask, includes medial wall
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
mw_L=zeros(1,2562);
mw_L(mwAndTSNR_L==1)=1;
mw_R=zeros(1,2562);
mw_R(mwAndTSNR_R==1)=1;
% convert to faces
F_MW_L=sum(mw_L(faces_l),2)./3;
F_MW_R=sum(mw_R(faces_r),2)./3;
% convert "partial" medial wall to medial wall
F_MW_L=ceil(F_MW_L);
F_MW_R=ceil(F_MW_R);
% face mask indices
fmwIndVec_l=find(F_MW_L);
fmwIndVec_r=find(F_MW_R);
% load in DMN to make more thorough mask
networks=load(['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/gro_Nets_fs4.mat']);
%% k = 1 to select DMN.
Dnet_LH=networks.nets.Lnets(:,1);
Dnet_RH=networks.nets.Rnets(:,1);
nets_LH=networks.nets.Lnets(:,1);
nets_RH=networks.nets.Rnets(:,1);
% create face-wise network mask
DMN_bool_L=sum(nets_LH(faces_l),2)./3;
DMN_bool_R=sum(nets_RH(faces_r),2)./3;
DMN_bool_L(DMN_bool_L>.3)=1;
DMN_bool_R(DMN_bool_R>.3)=1;
DMN_bool_L(DMN_bool_L<.3)=0;
DMN_bool_R(DMN_bool_R<.3)=0;
DMN_bool_L=logical(DMN_bool_L);
DMN_bool_R=logical(DMN_bool_R);
% combine with medial wall mask
MasterMask_L=DMN_bool_L;
MasterMask_R=DMN_bool_R;
MasterMask_L(fmwIndVec_l)=0;
MasterMask_R(fmwIndVec_r)=0;
% should be 1116 faces for left, 996 for right
mw_L=MasterMask_L;
mw_R=MasterMask_R;

% get subj list
subjPrefix=repmat('PS',10,1);
subjSuffix=["03","16","18","19","21","24","93","96","98","99"];
subjList=strcat(subjPrefix,subjSuffix')

% set common filepath
commonFP=['/scratch/users/apines/data/psil/'];

% to inherit same code as group-level scripting
s=find(strcmp(subj,subjList));
% now set subj to acutal subject name for filepaths
% should be redundant with input
subj=subjList(s)
% get session-condition correspondence
seshInfo=subSeshDose{1,s};
% initialize opfl structure to hold each condition/hemisphere
OpFl_rs_bv = struct();

%%%% now loop over each condition to load in each and concatenate resting-state angular time series
% get all baseline scans
baselineStrSearch_L=strjoin(['/scratch/users/apines/data/psil/' subj '/Baseline*/*_k1_Prop_TS_dmn_L.csv'],'');
baselineStrSearch_R=strjoin(['/scratch/users/apines/data/psil/' subj '/Baseline*/*_k1_Prop_TS_dmn_R.csv'],'');
bvscans_L=g_ls(char(baselineStrSearch_L));
bvscans_R=g_ls(char(baselineStrSearch_R));
% loop over all baseline %%
for c=1:length(bvscans_L);
	% extract Baseline # and task (rs1-rs6) from filename
	info=strsplit(bvscans_L{c},'/');
	sesh=info{8};
	filename=info{9};
	fileinfo=strsplit(filename,'_');
	task=fileinfo{3};
	% use extracted info to pull in number of remaining TRs
	survivingTrsFP=strjoin([commonFP '/' subj '/' sesh '/' subj '_' sesh '_task-' task '_ValidSegments_Trunc.txt'],'');
	survivingTrs=load(survivingTrsFP);
	% if remaining TRs > 250, concatenate it onto struct
	if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;	
		OpFl_rs_bv.(['f' num2str(c)])=struct();
		% LOAD IN ANGULAR TIME SERIES instead (and facewise)
		fpl=bvscans_L{c};
		fpr=bvscans_R{c};
		OpFl_rs_bv.(['f' num2str(c)]).L=dlmread(fpl);
		OpFl_rs_bv.(['f' num2str(c)]).R=dlmread(fpr);
	end
end

% initialize opfl structure to hold each condition/hemisphere
OpFl_rs_bw = struct();
% get all between scans
betweenStrSearch_L=strjoin(['/scratch/users/apines/data/psil/' subj '/Between*/*_k1_Prop_TS_dmn_L.csv'],'');
betweenStrSearch_R=strjoin(['/scratch/users/apines/data/psil/' subj '/Between*/*_k1_Prop_TS_dmn_R.csv'],'');
bwscans_L=g_ls(char(betweenStrSearch_L));
bwscans_R=g_ls(char(betweenStrSearch_R));
% loop over all between %%
for c=1:length(bwscans_L);
        % extract between # and task (rs1-rs6) from filename
        info=strsplit(bwscans_L{c},'/');
        sesh=info{8};
        filename=info{9};
        fileinfo=strsplit(filename,'_');
        task=fileinfo{3};
        % use extracted info to pull in number of remaining TRs
        survivingTrsFP=strjoin([commonFP '/' subj '/' sesh '/' subj '_' sesh '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        % if remaining TRs > 250, concatenate it onto struct
        if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250; 
                OpFl_rs_bw.(['f' num2str(c)])=struct();
                % LOAD IN ANGULAR TIME SERIES instead (and facewise)
                fpl=bwscans_L{c};
                fpr=bwscans_R{c};
                OpFl_rs_bw.(['f' num2str(c)]).L=dlmread(fpl);
                OpFl_rs_bw.(['f' num2str(c)]).R=dlmread(fpr);
        end
end

% loop over all after
% initialize opfl structure to hold each condition/hemisphere
OpFl_rs_af = struct();
% get all after scans
afterStrSearch_L=strjoin(['/scratch/users/apines/data/psil/' subj '/After*/*_k1_Prop_TS_dmn_L.csv'],'');
afterStrSearch_R=strjoin(['/scratch/users/apines/data/psil/' subj '/After*/*_k1_Prop_TS_dmn_R.csv'],'');
afscans_L=g_ls(char(afterStrSearch_L));
afscans_R=g_ls(char(afterStrSearch_R));
% loop over all after %%
for c=1:length(afscans_L);
        % extract after # and task (rs1-rs6) from filename
        info=strsplit(afscans_L{c},'/');
        sesh=info{8};
        filename=info{9};
        fileinfo=strsplit(filename,'_');
        task=fileinfo{3};
        % use extracted info to pull in number of remaining TRs
        survivingTrsFP=strjoin([commonFP '/' subj '/' sesh '/' subj '_' sesh '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        % if remaining TRs > 250, concatenate it onto struct
        if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;
                OpFl_rs_af.(['f' num2str(c)])=struct();
                % LOAD IN ANGULAR TIME SERIES instead (and facewise)
                fpl=afscans_L{c};
                fpr=afscans_R{c};
                OpFl_rs_af.(['f' num2str(c)]).L=dlmread(fpl);
                OpFl_rs_af.(['f' num2str(c)]).R=dlmread(fpr);
        end
end

% loop over all psil
% initialize opfl structure to hold each condition/hemisphere
OpFl_rs_p = struct();
% figure out which is psil and which is methyl for this PT
psilNum=subSeshDose{1,s};
% get all psil scans
psilStrSearch_L=strjoin(['/scratch/users/apines/data/psil/' subj '/Drug' num2str(psilNum) '/*_k1_Prop_TS_dmn_L.csv'],'');
psilStrSearch_R=strjoin(['/scratch/users/apines/data/psil/' subj '/Drug' num2str(psilNum) '/*_k1_Prop_TS_dmn_R.csv'],'');
pscans_L=g_ls(char(psilStrSearch_L));
pscans_R=g_ls(char(psilStrSearch_R));
% loop over all between %%
for c=1:length(pscans_L);
        % extract drug # and task (rs1-rs6) from filename
        info=strsplit(pscans_L{c},'/');
        sesh=info{8};
        filename=info{9};
        fileinfo=strsplit(filename,'_');
        task=fileinfo{3};
        % use extracted info to pull in number of remaining TRs
        survivingTrsFP=strjoin([commonFP '/' subj '/' sesh '/' subj '_' sesh '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        % if remaining TRs > 250, concatenate it onto struct
        if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;
                OpFl_rs_p.(['f' num2str(c)])=struct();
                % LOAD IN ANGULAR TIME SERIES instead (and facewise)
                fpl=pscans_L{c};
                fpr=pscans_R{c};
                OpFl_rs_p.(['f' num2str(c)]).L=dlmread(fpl);
                OpFl_rs_p.(['f' num2str(c)]).R=dlmread(fpr);
        end
end

% loop over all methyl
OpFl_rs_m = struct();
% figure out which one is methyl (drug)
methyNum=setdiff([1 2],psilNum);

% get all methyl scans
methStrSearch_L=strjoin(['/scratch/users/apines/data/psil/' subj '/Drug' num2str(methyNum) '/*_k1_Prop_TS_dmn_L.csv'],'');
methStrSearch_R=strjoin(['/scratch/users/apines/data/psil/' subj '/Drug' num2str(methyNum) '/*_k1_Prop_TS_dmn_R.csv'],'');
mscans_L=g_ls(char(methStrSearch_L));
mscans_R=g_ls(char(methStrSearch_R));
% loop over all between %%
for c=1:length(mscans_L);
        % extract drug # and task (rs1-rs6) from filename
        info=strsplit(mscans_L{c},'/');
        sesh=info{8};
        filename=info{9};
        fileinfo=strsplit(filename,'_');
        task=fileinfo{3};
        % use extracted info to pull in number of remaining TRs
        survivingTrsFP=strjoin([commonFP '/' subj '/' sesh '/' subj '_' sesh '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        % if remaining TRs > 250, concatenate it onto struct
        if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;
                OpFl_rs_m.(['f' num2str(c)])=struct();
                % LOAD IN ANGULAR TIME SERIES instead (and facewise)
                fpl=mscans_L{c};
                fpr=mscans_R{c};
                OpFl_rs_m.(['f' num2str(c)]).L=dlmread(fpl);
                OpFl_rs_m.(['f' num2str(c)]).R=dlmread(fpr);
        end
end

% saveout filepath
outFP=['/scratch/users/apines/data/psil/' subj];

% aggregate before scans
% Get the list of fields in the struct
fields = fieldnames(OpFl_rs_bv);
% Initialize matrices to accumulate the sum
sum_L = 0;
sum_R = 0;
% Loop over each field
for i = 1:length(fields)
    fieldName = fields{i};
    % Extract the L and R matrices
    bv_Angles_L = (sum(OpFl_rs_bv.(fieldName).L < 90,2) / size(OpFl_rs_bv.(fieldName).L,2))*100;
    bv_Angles_R = (sum(OpFl_rs_bv.(fieldName).R < 90,2) / size(OpFl_rs_bv.(fieldName).R,2))*100;
    % Accumulate the sum across fields
    sum_L = sum_L + bv_Angles_L;
    sum_R = sum_R + bv_Angles_R;
end
% average across each scan
avg_L = sum_L / length(fields);
avg_R = sum_R / length(fields);
% some conditions will be unpopulated
if (~isnan(avg_L))
	save(strjoin([outFP '/AvgBup_Bf_L.mat'],""),'avg_L');
	save(strjoin([outFP '/AvgBup_Bf_R.mat'],""),'avg_R');
end

% aggregate between scans
% Get the list of fields in the struct
fields = fieldnames(OpFl_rs_bw);
% Initialize matrices to accumulate the sum
sum_L = 0;
sum_R = 0;
% Loop over each field
for i = 1:length(fields)
    fieldName = fields{i};
    % Extract the L and R matrices
    bw_Angles_L = (sum(OpFl_rs_bw.(fieldName).L < 90,2) / size(OpFl_rs_bw.(fieldName).L,2))*100;
    bw_Angles_R = (sum(OpFl_rs_bw.(fieldName).R < 90,2) / size(OpFl_rs_bw.(fieldName).R,2))*100;
    % Accumulate the sum across fields
    sum_L = sum_L + bw_Angles_L;
    sum_R = sum_R + bw_Angles_R;
end
% average across each scan
avg_L = sum_L / length(fields);
avg_R = sum_R / length(fields);
if (~isnan(avg_L))
	save(strjoin([outFP '/AvgBup_Bw_L.mat'],""),'avg_L');
	save(strjoin([outFP '/AvgBup_Bw_R.mat'],""),'avg_R');
end

% aggregate after scans
% Get the list of fields in the struct
fields = fieldnames(OpFl_rs_af);
% Initialize matrices to accumulate the sum
sum_L = 0;
sum_R = 0;
% Loop over each field
for i = 1:length(fields)
    fieldName = fields{i};
    % Extract the L and R matrices
    af_Angles_L = (sum(OpFl_rs_af.(fieldName).L < 90,2) / size(OpFl_rs_af.(fieldName).L,2))*100;
    af_Angles_R = (sum(OpFl_rs_af.(fieldName).R < 90,2) / size(OpFl_rs_af.(fieldName).R,2))*100;
    % Accumulate the sum across fields
    sum_L = sum_L + af_Angles_L;
    sum_R = sum_R + af_Angles_R;
end
% average across each scan
avg_L = sum_L / length(fields);
avg_R = sum_R / length(fields);
if (~isnan(avg_L))
	save(strjoin([outFP '/AvgBup_Af_L.mat'],""),'avg_L');
	save(strjoin([outFP '/AvgBup_Af_R.mat'],""),'avg_R');
end
% aggregate psil scans
% Get the list of fields in the struct
fields = fieldnames(OpFl_rs_p);
% Initialize matrices to accumulate the sum
sum_L = 0;
sum_R = 0;
% Loop over each field
for i = 1:length(fields)
    fieldName = fields{i};
    % Extract the L and R matrices
    p_Angles_L = (sum(OpFl_rs_p.(fieldName).L < 90,2) / size(OpFl_rs_p.(fieldName).L,2))*100;
    p_Angles_R = (sum(OpFl_rs_p.(fieldName).R < 90,2) / size(OpFl_rs_p.(fieldName).R,2))*100;
    % Accumulate the sum across fields
    sum_L = sum_L + p_Angles_L;
    sum_R = sum_R + p_Angles_R;
end
% average across each scan
avg_L = sum_L / length(fields);
avg_R = sum_R / length(fields);
if (~isnan(avg_L))
	save(strjoin([outFP '/AvgBup_p_L.mat'],""),'avg_L');
	save(strjoin([outFP '/AvgBup_p_R.mat'],""),'avg_R');
end

% aggregate methyl scans
% Get the list of fields in the struct
fields = fieldnames(OpFl_rs_m);
% Initialize matrices to accumulate the sum
sum_L = 0;
sum_R = 0;
% Loop over each field
for i = 1:length(fields)
    fieldName = fields{i};
    % Extract the L and R matrices
    m_Angles_L = (sum(OpFl_rs_m.(fieldName).L < 90,2) / size(OpFl_rs_m.(fieldName).L,2))*100;
    m_Angles_R = (sum(OpFl_rs_m.(fieldName).R < 90,2) / size(OpFl_rs_m.(fieldName).R,2))*100;
    % Accumulate the sum across fields
    sum_L = sum_L + m_Angles_L;
    sum_R = sum_R + m_Angles_R;
end
% average across each scan
avg_L = sum_L / length(fields);
avg_R = sum_R / length(fields);
if (~isnan(avg_L))
	save(strjoin([outFP '/AvgBup_m_L.mat'],""),'avg_L');
	save(strjoin([outFP '/AvgBup_m_R.mat'],""),'avg_R');
end
