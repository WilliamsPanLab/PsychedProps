function Calc_AvgMagnitude_psil(subj)

% just need whether it's drug 1 or 2 that corresponds to psil
restoredefaultpath
subSeshDose=readtable('~/subjSeshDoseCorresp_psilo.csv');

% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/')
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder/PANDA_1.3.1_64')

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
% get mesh triangles
TR_L = TriRep(F_L,V_L);
P_L = TR_L;
TR_R = TriRep(F_R,V_R);
P_R = TR_R;
% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));
% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);


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
baselineStrSearch=strjoin([commonFP '/' subj '/Baseline*/*_OpFl.mat'],'');
bvscans=g_ls(char(baselineStrSearch));
% loop over all baseline %%
for c=1:length(bvscans);
	% extract Baseline # and task (rs1-rs6) from filename
	info=strsplit(bvscans{c},'/');
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
		% load in optical flow output
		OpFl_rs_bv.(['f' num2str(c)]).L=load(bvscans{c}).us.vf_left;
		OpFl_rs_bv.(['f' num2str(c)]).R=load(bvscans{c}).us.vf_right;
	end
end

% initialize opfl structure to hold each condition/hemisphere
OpFl_rs_bw = struct();
% get all between scans
betweenStrSearch=strjoin([commonFP '/' subj '/Between*/*_OpFl.mat'],'');
bwscans=g_ls(char(betweenStrSearch));
% loop over all between %%
for c=1:length(bwscans);
        % extract between # and task (rs1-rs6) from filename
        info=strsplit(bwscans{c},'/');
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
                % load in optical flow output
                OpFl_rs_bw.(['f' num2str(c)]).L=load(bwscans{c}).us.vf_left;
                OpFl_rs_bw.(['f' num2str(c)]).R=load(bwscans{c}).us.vf_right;
        end
end

% loop over all after
% initialize opfl structure to hold each condition/hemisphere
OpFl_rs_af = struct();
% get all after scans
afterStrSearch=strjoin([commonFP '/' subj '/After*/*_OpFl.mat'],'');
afscans=g_ls(char(afterStrSearch));
% loop over all after %%
for c=1:length(afscans);
        % extract after # and task (rs1-rs6) from filename
        info=strsplit(afscans{c},'/');
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
                % load in optical flow output
                OpFl_rs_af.(['f' num2str(c)]).L=load(afscans{c}).us.vf_left;
                OpFl_rs_af.(['f' num2str(c)]).R=load(afscans{c}).us.vf_right;
        end
end

% loop over all psil
% initialize opfl structure to hold each condition/hemisphere
OpFl_rs_p = struct();
% figure out which is psil and which is methyl for this PT
psilNum=subSeshDose{1,s};
% get all psil scans
psilStrSearch=strjoin([commonFP '/' subj '/Drug' num2str(psilNum) '/*_OpFl.mat'],'');
pscans=g_ls(char(psilStrSearch));
% loop over all psil %%
for c=1:length(pscans);
        % extract drug # and task (rs1-rs6) from filename
        info=strsplit(pscans{c},'/');
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
                % load in optical flow output
                OpFl_rs_p.(['f' num2str(c)]).L=load(pscans{c}).us.vf_left;
                OpFl_rs_p.(['f' num2str(c)]).R=load(pscans{c}).us.vf_right;
        end
end

% loop over all methyl
OpFl_rs_m = struct();
% figure out which one is methyl (drug)
methyNum=setdiff([1 2],psilNum);

% get all methyl scans
methStrSearch=strjoin([commonFP '/' subj '/Drug' num2str(methyNum) '/*_OpFl.mat'],'');
mscans=g_ls(char(methStrSearch));
% loop over all meth %%
for c=1:length(mscans);
        % extract drug # and task (rs1-rs6) from filename
        info=strsplit(mscans{c},'/');
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
		% load in optical flow output
                OpFl_rs_m.(['f' num2str(c)]).L=load(mscans{c}).us.vf_left;
                OpFl_rs_m.(['f' num2str(c)]).R=load(mscans{c}).us.vf_right;                
        end
end

% saveout filepath
outFP=['/scratch/users/apines/data/psil/' subj];

% initialize 2d vector arrays
bv_L=zeros(5120,1);
bv_R=zeros(5120,1);
bw_L=zeros(5120,1);
bw_R=zeros(5120,1);
af_L=zeros(5120,1);
af_R=zeros(5120,1);
m_L=zeros(5120,1);
m_R=zeros(5120,1);
p_L=zeros(5120,1);
p_R=zeros(5120,1);

% remove cell structure for more straightforward array
% Loop through each field in OpFl_rs_bv
for field = fieldnames(OpFl_rs_bv)'
    f = field{1}; 
    % Get the number of cells
    num_cells = numel(OpFl_rs_bv.(f).L);
    % Preallocate arrays for L and R with the correct size
    L_array = zeros(5120, num_cells, 3);
    R_array = zeros(5120, num_cells, 3);
    % Populate the arrays by preserving the original correspondence
    for j = 1:num_cells
        L_array(:, j, :) = OpFl_rs_bv.(f).L{j};
        R_array(:, j, :) = OpFl_rs_bv.(f).R{j};
    end
    % Update the struct with the new arrays
    OpFl_rs_bv.(f).L = L_array;
    OpFl_rs_bv.(f).R = R_array;
end
% Loop through each field in OpFl_rs_bw
for field = fieldnames(OpFl_rs_bw)'
    f = field{1};
    % Get the number of cells
    num_cells = numel(OpFl_rs_bw.(f).L);
    % Preallocate arrays for L and R with the correct size
    L_array = zeros(5120, num_cells, 3);
    R_array = zeros(5120, num_cells, 3);
    % Populate the arrays by preserving the original correspondence
    for j = 1:num_cells
        L_array(:, j, :) = OpFl_rs_bw.(f).L{j};
        R_array(:, j, :) = OpFl_rs_bw.(f).R{j};
    end
    % Update the struct with the new arrays
    OpFl_rs_bw.(f).L = L_array;
    OpFl_rs_bw.(f).R = R_array;
end
% Loop through each field in OpFl_rs_af
for field = fieldnames(OpFl_rs_af)'
    f = field{1};
    % Get the number of cells
    num_cells = numel(OpFl_rs_af.(f).L);
    % Preallocate arrays for L and R with the correct size
    L_array = zeros(5120, num_cells, 3);
    R_array = zeros(5120, num_cells, 3);
    % Populate the arrays by preserving the original correspondence
    for j = 1:num_cells
        L_array(:, j, :) = OpFl_rs_af.(f).L{j};
        R_array(:, j, :) = OpFl_rs_af.(f).R{j};
    end
    % Update the struct with the new arrays
    OpFl_rs_af.(f).L = L_array;
    OpFl_rs_af.(f).R = R_array;
end
% Loop through each field in OpFl_rs_p
for field = fieldnames(OpFl_rs_p)'
    f = field{1};
    % Get the number of cells
    num_cells = numel(OpFl_rs_p.(f).L);
    % Preallocate arrays for L and R with the correct size
    L_array = zeros(5120, num_cells, 3);
    R_array = zeros(5120, num_cells, 3);
    % Populate the arrays by preserving the original correspondence
    for j = 1:num_cells
        L_array(:, j, :) = OpFl_rs_p.(f).L{j};
        R_array(:, j, :) = OpFl_rs_p.(f).R{j};
    end
    % Update the struct with the new arrays
    OpFl_rs_p.(f).L = L_array;
    OpFl_rs_p.(f).R = R_array;
end
% Loop through each field in OpFl_rs_m
for field = fieldnames(OpFl_rs_m)'
    f = field{1};
    % Get the number of cells
    num_cells = numel(OpFl_rs_m.(f).L);
    % Preallocate arrays for L and R with the correct size
    L_array = zeros(5120, num_cells, 3);
    R_array = zeros(5120, num_cells, 3);
    % Populate the arrays by preserving the original correspondence
    for j = 1:num_cells
        L_array(:, j, :) = OpFl_rs_m.(f).L{j};
        R_array(:, j, :) = OpFl_rs_m.(f).R{j};
    end
    % Update the struct with the new arrays
    OpFl_rs_m.(f).L = L_array;
    OpFl_rs_m.(f).R = R_array;
end


% aggregate before scans
% Get the list of fields in the struct
fields = fieldnames(OpFl_rs_bv);
% Loop over each field
for i = 1:length(fields)
   fieldName = fields{i};
   % initialize facewise vector for this field
   fwise_vec_L=zeros(5120,1);
   fwise_vec_R=zeros(5120,1);
   % loop over each vertex: left hemisphere
   for f=1:5120
   	% Extract the L and R matrices
   	bv_vecs_L = OpFl_rs_bv.(fieldName).L;
   	% get the average magnitude per vertex
   	x_comps_L = bv_vecs_L(f,:,1);
	y_comps_L = bv_vecs_L(f,:,2);
	z_comps_L = bv_vecs_L(f,:,3);
	% for each timepoint
	for tp=1:size(bv_vecs_L,2);
		x_comp_L=x_comps_L(tp);
		y_comp_L=y_comps_L(tp);
		z_comp_L=z_comps_L(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_L;y_comp_L;z_comp_L]),azd_L(f),eld_L(f));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		fwise_vec_L(f)=fwise_vec_L(f)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    fwise_vec_L=fwise_vec_L./(size(bv_vecs_L,2));
    % right hemisphere
    for f=1:5120
	% Extract the L and R matrices
	bv_vecs_R = OpFl_rs_bv.(fieldName).R;
	% get the average magnitude per vertex
	x_comps_R = bv_vecs_R(f,:,1);
	y_comps_R = bv_vecs_R(f,:,2);
	z_comps_R = bv_vecs_R(f,:,3);
	% for each timepoint
	for tp=1:size(bv_vecs_R,2);
		x_comp_R=x_comps_R(tp);
		y_comp_R=y_comps_R(tp);
		z_comp_R=z_comps_R(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_R;y_comp_R;z_comp_R]),azd_R(f),eld_R(f));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		fwise_vec_R(f)=fwise_vec_R(f)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    fwise_vec_R=fwise_vec_R./(size(bv_vecs_R,2));
    % add in vwise vecs to bv_L and bv_R
    bv_L=bv_L+fwise_vec_L;
    bv_R=bv_R+fwise_vec_R;
end
% average across each scan
avg_L = bv_L / length(fields);
avg_R = bv_R / length(fields);

% some conditions will be unpopulated
if (~isnan(avg_L))
	save(strjoin([outFP '/AvgMag_Bf_L.mat'],""),'avg_L');
	save(strjoin([outFP '/AvgMag_Bf_R.mat'],""),'avg_R');
end

% aggregate between scans
% Get the list of fields in the struct
fields = fieldnames(OpFl_rs_bw);
% Loop over each field
for i = 1:length(fields)
    fieldName = fields{i};
    % initialize vertexwise vector for this field
    fwise_vec_L=zeros(5120,1);
    fwise_vec_R=zeros(5120,1);
    % loop over each vertex: left hemisphere
    for f=1:5120
	% Extract the L and R matrices
	bw_vecs_L = OpFl_rs_bw.(fieldName).L;
	% get the average magnitude per vertex
	x_comps_L = bw_vecs_L(f,:,1);
	y_comps_L = bw_vecs_L(f,:,2);
	z_comps_L = bw_vecs_L(f,:,3);
	% for each timepoint
	for tp=1:size(bw_vecs_L,2);
		x_comp_L=x_comps_L(tp);
		y_comp_L=y_comps_L(tp);
		z_comp_L=z_comps_L(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_L;y_comp_L;z_comp_L]),azd_L(f),eld_L(f));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		fwise_vec_L(f)=fwise_vec_L(f)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    fwise_vec_L=fwise_vec_L./(size(bw_vecs_L,2));
    % right hemisphere
    for f=1:5120
	% Extract the L and R matrices
	bw_vecs_R = OpFl_rs_bw.(fieldName).R;
	% get the average magnitude per vertex
	x_comps_R = bw_vecs_R(f,:,1);
	y_comps_R = bw_vecs_R(f,:,2);
	z_comps_R = bw_vecs_R(f,:,3);
	% for each timepoint
	for tp=1:size(bw_vecs_R,2);
		x_comp_R=x_comps_R(tp);
		y_comp_R=y_comps_R(tp);
		z_comp_R=z_comps_R(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_R;y_comp_R;z_comp_R]),azd_R(f),eld_R(f));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		fwise_vec_R(f)=fwise_vec_R(f)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    fwise_vec_R=fwise_vec_R./(size(bw_vecs_R,2));
    % add in vwise vecs to bw_L and bw_R
    bw_L=bw_L+fwise_vec_L;
    bw_R=bw_R+fwise_vec_R;
end
% average across each scan
avg_L = bw_L / length(fields);
avg_R = bw_R / length(fields);

% saveout
if (~isnan(avg_L))
	save(strjoin([outFP '/AvgMag_Bw_L.mat'],""),'avg_L');
	save(strjoin([outFP '/AvgMag_Bw_R.mat'],""),'avg_R');
end

% aggregate after scans
% Get the list of fields in the struct
fields = fieldnames(OpFl_rs_af);
% Loop over each field
for i = 1:length(fields)
    fieldName = fields{i};
    % initialize vertexwise vector for this field
    fwise_vec_L=zeros(5120,1);
    fwise_vec_R=zeros(5120,1);
    % loop over each vertex: left hemisphere
    for f=1:5120
	% Extract the L and R matrices
	af_vecs_L = OpFl_rs_af.(fieldName).L;
	% get the average magnitude per vertex
	x_comps_L = af_vecs_L(f,:,1);
	y_comps_L = af_vecs_L(f,:,2);
	z_comps_L = af_vecs_L(f,:,3);
	% for each timepoint
	for tp=1:size(af_vecs_L,2);
		x_comp_L=x_comps_L(tp);
		y_comp_L=y_comps_L(tp);
		z_comp_L=z_comps_L(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_L;y_comp_L;z_comp_L]),azd_L(f),eld_L(f));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		fwise_vec_L(f)=fwise_vec_L(f)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    fwise_vec_L=fwise_vec_L./(size(af_vecs_L,2));
    % right hemisphere
    for f=1:5120
	% Extract the L and R matrices
	af_vecs_R = OpFl_rs_af.(fieldName).R;
	% get the average magnitude per vertex
	x_comps_R = af_vecs_R(f,:,1);
	y_comps_R = af_vecs_R(f,:,2);
	z_comps_R = af_vecs_R(f,:,3);
	% for each timepoint
	for tp=1:size(af_vecs_R,2);
		x_comp_R=x_comps_R(tp);
		y_comp_R=y_comps_R(tp);
		z_comp_R=z_comps_R(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_R;y_comp_R;z_comp_R]),azd_R(f),eld_R(f));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		fwise_vec_R(f)=fwise_vec_R(f)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    fwise_vec_R=fwise_vec_R./(size(af_vecs_R,2));
    % add in vwise vecs to af_L and af_R
    af_L=af_L+fwise_vec_L;
    af_R=af_R+fwise_vec_R;
end

% average across each scan
avg_L = af_L / length(fields);
avg_R = af_R / length(fields);
% saveout
if (~isnan(avg_L))
	save(strjoin([outFP '/AvgMag_Af_L.mat'],""),'avg_L');
	save(strjoin([outFP '/AvgMag_Af_R.mat'],""),'avg_R');
end
% aggregate psil scans
% Get the list of fields in the struct
fields = fieldnames(OpFl_rs_p);
% Loop over each field
for i = 1:length(fields)
    fieldName = fields{i};
    % initialize vertexwise vector for this field
    fwise_vec_L=zeros(5120,1);
    fwise_vec_R=zeros(5120,1);
    % loop over each vertex: left hemisphere
    for f=1:5120
	% Extract the L and R matrices
	p_vecs_L = OpFl_rs_p.(fieldName).L;
	% get the average magnitude per vertex
	x_comps_L = p_vecs_L(f,:,1);
	y_comps_L = p_vecs_L(f,:,2);
	z_comps_L = p_vecs_L(f,:,3);
	% for each timepoint
	for tp=1:size(p_vecs_L,2);
		x_comp_L=x_comps_L(tp);
		y_comp_L=y_comps_L(tp);
		z_comp_L=z_comps_L(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_L;y_comp_L;z_comp_L]),azd_L(f),eld_L(f));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		fwise_vec_L(f)=fwise_vec_L(f)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    fwise_vec_L=fwise_vec_L./(size(p_vecs_L,2));
    % right hemisphere
    for f=1:5120
	% Extract the L and R matrices
	p_vecs_R = OpFl_rs_p.(fieldName).R;
	% get the average magnitude per vertex
	x_comps_R = p_vecs_R(f,:,1);
	y_comps_R = p_vecs_R(f,:,2);
	z_comps_R = p_vecs_R(f,:,3);
	% for each timepoint
	for tp=1:size(p_vecs_R,2);
		x_comp_R=x_comps_R(tp);
		y_comp_R=y_comps_R(tp);
		z_comp_R=z_comps_R(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_R;y_comp_R;z_comp_R]),azd_R(f),eld_R(f));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		fwise_vec_R(f)=fwise_vec_R(f)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    fwise_vec_R=fwise_vec_R./(size(p_vecs_R,2));
    % add in vwise vecs to p_L and p_R
    p_L=p_L+fwise_vec_L;
    p_R=p_R+fwise_vec_R;
end
% average across each scan
avg_L = p_L / length(fields);
avg_R = p_R / length(fields);
% saveout
if (~isnan(avg_L))
	save(strjoin([outFP '/AvgMag_p_L.mat'],""),'avg_L');
	save(strjoin([outFP '/AvgMag_p_R.mat'],""),'avg_R');
end

% aggregate methyl scans
% Get the list of fields in the struct
fields = fieldnames(OpFl_rs_m);
% Loop over each field
for i = 1:length(fields)
    fieldName = fields{i};
    % initialize vertexwise vector for this field
    fwise_vec_L=zeros(5120,1);
    fwise_vec_R=zeros(5120,1);
    % loop over each vertex: left hemisphere
    for f=1:5120
	% Extract the L and R matrices
	m_vecs_L = OpFl_rs_m.(fieldName).L;
	% get the average magnitude per vertex
	x_comps_L = m_vecs_L(f,:,1);
	y_comps_L = m_vecs_L(f,:,2);
	z_comps_L = m_vecs_L(f,:,3);
	% for each timepoint
	for tp=1:size(m_vecs_L,2);
		x_comp_L=x_comps_L(tp);
		y_comp_L=y_comps_L(tp);
		z_comp_L=z_comps_L(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_L;y_comp_L;z_comp_L]),azd_L(f),eld_L(f));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		fwise_vec_L(f)=fwise_vec_L(f)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    fwise_vec_L=fwise_vec_L./(size(m_vecs_L,2));
    % right hemisphere
    for f=1:5120
	% Extract the L and R matrices
	m_vecs_R = OpFl_rs_m.(fieldName).R;
	% get the average magnitude per vertex
	x_comps_R = m_vecs_R(f,:,1);
	y_comps_R = m_vecs_R(f,:,2);
	z_comps_R = m_vecs_R(f,:,3);
	% for each timepoint
	for tp=1:size(m_vecs_R,2);
		x_comp_R=x_comps_R(tp);
		y_comp_R=y_comps_R(tp);
		z_comp_R=z_comps_R(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_R;y_comp_R;z_comp_R]),azd_R(f),eld_R(f));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		fwise_vec_R(f)=fwise_vec_R(f)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    fwise_vec_R=fwise_vec_R./(size(m_vecs_R,2));
    % add in vwise vecs to m_L and m_R
    m_L=m_L+fwise_vec_L;
    m_R=m_R+fwise_vec_R;
end
% average across each scan
avg_L = m_L / length(fields);
avg_R = m_R / length(fields);
% saveout
if (~isnan(avg_L))
	save(strjoin([outFP '/AvgMag_m_L.mat'],""),'avg_L');
	save(strjoin([outFP '/AvgMag_m_R.mat'],""),'avg_R');
end
