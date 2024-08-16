function Calc_AvgMagnitude_psil(subj)

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
%%% need this to get more surface geometry information
% get incenters of triangles
TR_L = TriRep(F_L,V_L);
% note X returns vertices rather than incenters of faces
P_L = TR_L.X;
TR_R = TriRep(F_R,V_R);
P_R = TR_R.X;
% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));
% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);

% load in medial wall + SNR so we don't have to loop over every single vertex and then mask out later
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
% note this is indexing for VALID vertices as opposed to some other scripts with 1 at INVALID vertices
mw_L=ones(1,2562);
mw_L(mwAndTSNR_L>0)=0;
mw_R=ones(1,2562);
mw_R(mwAndTSNR_R>0)=0;
mw_L=logical(mw_L);
mw_R=logical(mw_R);

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
baselineStrSearch_L=strjoin([commonFP '/' subj '/Baseline*/*_p2mm_masked_Vert_Angles_L.mat'],'');
baselineStrSearch_R=strjoin([commonFP '/' subj '/Baseline*/*_p2mm_masked_Vert_Angles_R.mat'],'');
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
	survivingTrsFP=strjoin([commonFP '/' subj '/' sesh '/' subj '_' sesh '_' task '_ValidSegments_Trunc.txt'],'');
	survivingTrs=load(survivingTrsFP);
	% if remaining TRs > 250, concatenate it onto struct
	if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;	
		OpFl_rs_bv.(['f' num2str(c)])=struct();
		% LOAD IN ANGULAR TIME SERIES instead (and facewise)
		fpl=bvscans_L{c};
		fpr=bvscans_R{c};
		OpFl_rs_bv.(['f' num2str(c)]).L=load(fpl).vertWise_Vecs_l(mw_L,:,:);
		OpFl_rs_bv.(['f' num2str(c)]).R=load(fpr).vertWise_Vecs_r(mw_R,:,:);
	end
end

% initialize opfl structure to hold each condition/hemisphere
OpFl_rs_bw = struct();
% get all between scans
betweenStrSearch_L=strjoin([commonFP '/' subj '/Between*/*_p2mm_masked_Vert_Angles_L.mat'],'');
betweenStrSearch_R=strjoin([commonFP '/' subj '/Between*/*_p2mm_masked_Vert_Angles_R.mat'],'');
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
        survivingTrsFP=strjoin([commonFP '/' subj '/' sesh '/' subj '_' sesh '_' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        % if remaining TRs > 250, concatenate it onto struct
        if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250; 
                OpFl_rs_bw.(['f' num2str(c)])=struct();
                % LOAD IN ANGULAR TIME SERIES instead (and facewise)
                fpl=bwscans_L{c};
                fpr=bwscans_R{c};
                OpFl_rs_bw.(['f' num2str(c)]).L=load(fpl).vertWise_Vecs_l(mw_L,:,:);
                OpFl_rs_bw.(['f' num2str(c)]).R=load(fpr).vertWise_Vecs_r(mw_R,:,:);
        end
end

% loop over all after
% initialize opfl structure to hold each condition/hemisphere
OpFl_rs_af = struct();
% get all after scans
afterStrSearch_L=strjoin([commonFP '/' subj '/After*/*_p2mm_masked_Vert_Angles_L.mat'],'');
afterStrSearch_R=strjoin([commonFP '/' subj '/After*/*_p2mm_masked_Vert_Angles_R.mat'],'');
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
        survivingTrsFP=strjoin([commonFP '/' subj '/' sesh '/' subj '_' sesh '_' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        % if remaining TRs > 250, concatenate it onto struct
        if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;
                OpFl_rs_af.(['f' num2str(c)])=struct();
                % LOAD IN ANGULAR TIME SERIES instead (and facewise)
                fpl=afscans_L{c};
                fpr=afscans_R{c};
                OpFl_rs_af.(['f' num2str(c)]).L=load(fpl).vertWise_Vecs_l(mw_L,:,:);
                OpFl_rs_af.(['f' num2str(c)]).R=load(fpr).vertWise_Vecs_r(mw_R,:,:);
        end
end

% loop over all psil
% initialize opfl structure to hold each condition/hemisphere
OpFl_rs_p = struct();
% figure out which is psil and which is methyl for this PT
psilNum=subSeshDose{1,s};
% get all psil scans
psilStrSearch_L=strjoin([commonFP '/' subj '/Drug' num2str(psilNum) '/*_p2mm_masked_Vert_Angles_L.mat'],'');
psilStrSearch_R=strjoin([commonFP '/' subj '/Drug' num2str(psilNum) '/*_p2mm_masked_Vert_Angles_R.mat'],'');
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
        survivingTrsFP=strjoin([commonFP '/' subj '/' sesh '/' subj '_' sesh '_' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        % if remaining TRs > 250, concatenate it onto struct
        if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;
                OpFl_rs_p.(['f' num2str(c)])=struct();
                % LOAD IN ANGULAR TIME SERIES instead (and facewise)
                fpl=pscans_L{c};
                fpr=pscans_R{c};
                OpFl_rs_p.(['f' num2str(c)]).L=load(fpl).vertWise_Vecs_l(mw_L,:,:);
                OpFl_rs_p.(['f' num2str(c)]).R=load(fpr).vertWise_Vecs_r(mw_R,:,:);
        end
end

% loop over all methyl
OpFl_rs_m = struct();
% figure out which one is methyl (drug)
methyNum=setdiff([1 2],psilNum);

% get all methyl scans
methStrSearch_L=strjoin([commonFP '/' subj '/Drug' num2str(methyNum) '/*_p2mm_masked_Vert_Angles_L.mat'],'');
methStrSearch_R=strjoin([commonFP '/' subj '/Drug' num2str(methyNum) '/*_p2mm_masked_Vert_Angles_R.mat'],'');
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
        survivingTrsFP=strjoin([commonFP '/' subj '/' sesh '/' subj '_' sesh '_' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        % if remaining TRs > 250, concatenate it onto struct
        if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;
                OpFl_rs_m.(['f' num2str(c)])=struct();
                % LOAD IN ANGULAR TIME SERIES instead (and facewise)
                fpl=mscans_L{c};
                fpr=mscans_R{c};
                OpFl_rs_m.(['f' num2str(c)]).L=load(fpl).vertWise_Vecs_l(mw_L,:,:);
                OpFl_rs_m.(['f' num2str(c)]).R=load(fpr).vertWise_Vecs_r(mw_R,:,:);
        end
end

% saveout filepath
outFP=['/scratch/users/apines/data/psil/' subj];

% create valid verts of surface geometry
valid_verts_L=find(mw_L);
valid_verts_R=find(mw_R);

% initialize 2d vector arrays
bv_L=zeros(sum(mw_L),1);
bv_R=zeros(sum(mw_R),1);
bw_L=zeros(sum(mw_L),1);
bw_R=zeros(sum(mw_R),1);
af_L=zeros(sum(mw_L),1);
af_R=zeros(sum(mw_R),1);
m_L=zeros(sum(mw_L),1);
m_R=zeros(sum(mw_R),1);
p_L=zeros(sum(mw_L),1);
p_R=zeros(sum(mw_R),1);

% aggregate before scans
% Get the list of fields in the struct
fields = fieldnames(OpFl_rs_bv);
% Loop over each field
for i = 1:length(fields)
   fieldName = fields{i};
   % initialize vertexwise vector for this field
   vwise_vec_L=zeros(sum(mw_L),1);
   vwise_vec_R=zeros(sum(mw_R),1);
   % loop over each vertex: left hemisphere
   for v=1:sum(mw_L)
   	% Extract the L and R matrices
   	bv_vecs_L = OpFl_rs_bv.(fieldName).L;
   	% get the average magnitude per vertex
   	x_comps_L = bv_vecs_L(v,:,1);
	y_comps_L = bv_vecs_L(v,:,2);
	z_comps_L = bv_vecs_L(v,:,3);
	% for each timepoint
	for tp=1:size(bv_vecs_L,2);
		x_comp_L=x_comps_L(tp);
		y_comp_L=y_comps_L(tp);
		z_comp_L=z_comps_L(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_L;y_comp_L;z_comp_L]),azd_L(valid_verts_L(v)),eld_L(valid_verts_L(v)));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		vwise_vec_L(v)=vwise_vec_L(v)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    vwise_vec_L=vwise_vec_L./(size(bv_vecs_L,2));
    % right hemisphere
    for v=1:sum(mw_R)
	% Extract the L and R matrices
	bv_vecs_R = OpFl_rs_bv.(fieldName).R;
	% get the average magnitude per vertex
	x_comps_R = bv_vecs_R(v,:,1);
	y_comps_R = bv_vecs_R(v,:,2);
	z_comps_R = bv_vecs_R(v,:,3);
	% for each timepoint
	for tp=1:size(bv_vecs_R,2);
		x_comp_R=x_comps_R(tp);
		y_comp_R=y_comps_R(tp);
		z_comp_R=z_comps_R(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_R;y_comp_R;z_comp_R]),azd_R(valid_verts_R(v)),eld_R(valid_verts_R(v)));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		vwise_vec_R(v)=vwise_vec_R(v)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    vwise_vec_R=vwise_vec_R./(size(bv_vecs_R,2));
    % add in vwise vecs to bv_L and bv_R
    bv_L=bv_L+vwise_vec_L;
    bv_R=bv_R+vwise_vec_R;
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
    vwise_vec_L=zeros(sum(mw_L),1);
    vwise_vec_R=zeros(sum(mw_R),1);
    % loop over each vertex: left hemisphere
    for v=1:sum(mw_L)
	% Extract the L and R matrices
	bw_vecs_L = OpFl_rs_bw.(fieldName).L;
	% get the average magnitude per vertex
	x_comps_L = bw_vecs_L(v,:,1);
	y_comps_L = bw_vecs_L(v,:,2);
	z_comps_L = bw_vecs_L(v,:,3);
	% for each timepoint
	for tp=1:size(bw_vecs_L,2);
		x_comp_L=x_comps_L(tp);
		y_comp_L=y_comps_L(tp);
		z_comp_L=z_comps_L(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_L;y_comp_L;z_comp_L]),azd_L(valid_verts_L(v)),eld_L(valid_verts_L(v)));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		vwise_vec_L(v)=vwise_vec_L(v)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    vwise_vec_L=vwise_vec_L./(size(bw_vecs_L,2));
    % right hemisphere
    for v=1:sum(mw_R)
	% Extract the L and R matrices
	bw_vecs_R = OpFl_rs_bw.(fieldName).R;
	% get the average magnitude per vertex
	x_comps_R = bw_vecs_R(v,:,1);
	y_comps_R = bw_vecs_R(v,:,2);
	z_comps_R = bw_vecs_R(v,:,3);
	% for each timepoint
	for tp=1:size(bw_vecs_R,2);
		x_comp_R=x_comps_R(tp);
		y_comp_R=y_comps_R(tp);
		z_comp_R=z_comps_R(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_R;y_comp_R;z_comp_R]),azd_R(valid_verts_R(v)),eld_R(valid_verts_R(v)));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		vwise_vec_R(v)=vwise_vec_R(v)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    vwise_vec_R=vwise_vec_R./(size(bw_vecs_R,2));
    % add in vwise vecs to bw_L and bw_R
    bw_L=bw_L+vwise_vec_L;
    bw_R=bw_R+vwise_vec_R;
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
    vwise_vec_L=zeros(sum(mw_L),1);
    vwise_vec_R=zeros(sum(mw_R),1);
    % loop over each vertex: left hemisphere
    for v=1:sum(mw_L)
	% Extract the L and R matrices
	af_vecs_L = OpFl_rs_af.(fieldName).L;
	% get the average magnitude per vertex
	x_comps_L = af_vecs_L(v,:,1);
	y_comps_L = af_vecs_L(v,:,2);
	z_comps_L = af_vecs_L(v,:,3);
	% for each timepoint
	for tp=1:size(af_vecs_L,2);
		x_comp_L=x_comps_L(tp);
		y_comp_L=y_comps_L(tp);
		z_comp_L=z_comps_L(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_L;y_comp_L;z_comp_L]),azd_L(valid_verts_L(v)),eld_L(valid_verts_L(v)));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		vwise_vec_L(v)=vwise_vec_L(v)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    vwise_vec_L=vwise_vec_L./(size(af_vecs_L,2));
    % right hemisphere
    for v=1:sum(mw_R)
	% Extract the L and R matrices
	af_vecs_R = OpFl_rs_af.(fieldName).R;
	% get the average magnitude per vertex
	x_comps_R = af_vecs_R(v,:,1);
	y_comps_R = af_vecs_R(v,:,2);
	z_comps_R = af_vecs_R(v,:,3);
	% for each timepoint
	for tp=1:size(af_vecs_R,2);
		x_comp_R=x_comps_R(tp);
		y_comp_R=y_comps_R(tp);
		z_comp_R=z_comps_R(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_R;y_comp_R;z_comp_R]),azd_R(valid_verts_R(v)),eld_R(valid_verts_R(v)));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		vwise_vec_R(v)=vwise_vec_R(v)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    vwise_vec_R=vwise_vec_R./(size(af_vecs_R,2));
    % add in vwise vecs to af_L and af_R
    af_L=af_L+vwise_vec_L;
    af_R=af_R+vwise_vec_R;
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
    vwise_vec_L=zeros(sum(mw_L),1);
    vwise_vec_R=zeros(sum(mw_R),1);
    % loop over each vertex: left hemisphere
    for v=1:sum(mw_L)
	% Extract the L and R matrices
	p_vecs_L = OpFl_rs_p.(fieldName).L;
	% get the average magnitude per vertex
	x_comps_L = p_vecs_L(v,:,1);
	y_comps_L = p_vecs_L(v,:,2);
	z_comps_L = p_vecs_L(v,:,3);
	% for each timepoint
	for tp=1:size(p_vecs_L,2);
		x_comp_L=x_comps_L(tp);
		y_comp_L=y_comps_L(tp);
		z_comp_L=z_comps_L(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_L;y_comp_L;z_comp_L]),azd_L(valid_verts_L(v)),eld_L(valid_verts_L(v)));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		vwise_vec_L(v)=vwise_vec_L(v)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    vwise_vec_L=vwise_vec_L./(size(p_vecs_L,2));
    % right hemisphere
    for v=1:sum(mw_R)
	% Extract the L and R matrices
	p_vecs_R = OpFl_rs_p.(fieldName).R;
	% get the average magnitude per vertex
	x_comps_R = p_vecs_R(v,:,1);
	y_comps_R = p_vecs_R(v,:,2);
	z_comps_R = p_vecs_R(v,:,3);
	% for each timepoint
	for tp=1:size(p_vecs_R,2);
		x_comp_R=x_comps_R(tp);
		y_comp_R=y_comps_R(tp);
		z_comp_R=z_comps_R(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_R;y_comp_R;z_comp_R]),azd_R(valid_verts_R(v)),eld_R(valid_verts_R(v)));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		vwise_vec_R(v)=vwise_vec_R(v)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    vwise_vec_R=vwise_vec_R./(size(p_vecs_R,2));
    % add in vwise vecs to p_L and p_R
    p_L=p_L+vwise_vec_L;
    p_R=p_R+vwise_vec_R;
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
    vwise_vec_L=zeros(sum(mw_L),1);
    vwise_vec_R=zeros(sum(mw_R),1);
    % loop over each vertex: left hemisphere
    for v=1:sum(mw_L)
	% Extract the L and R matrices
	m_vecs_L = OpFl_rs_m.(fieldName).L;
	% get the average magnitude per vertex
	x_comps_L = m_vecs_L(v,:,1);
	y_comps_L = m_vecs_L(v,:,2);
	z_comps_L = m_vecs_L(v,:,3);
	% for each timepoint
	for tp=1:size(m_vecs_L,2);
		x_comp_L=x_comps_L(tp);
		y_comp_L=y_comps_L(tp);
		z_comp_L=z_comps_L(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_L;y_comp_L;z_comp_L]),azd_L(valid_verts_L(v)),eld_L(valid_verts_L(v)));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		vwise_vec_L(v)=vwise_vec_L(v)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    vwise_vec_L=vwise_vec_L./(size(m_vecs_L,2));
    % right hemisphere
    for v=1:sum(mw_R)
	% Extract the L and R matrices
	m_vecs_R = OpFl_rs_m.(fieldName).R;
	% get the average magnitude per vertex
	x_comps_R = m_vecs_R(v,:,1);
	y_comps_R = m_vecs_R(v,:,2);
	z_comps_R = m_vecs_R(v,:,3);
	% for each timepoint
	for tp=1:size(m_vecs_R,2);
		x_comp_R=x_comps_R(tp);
		y_comp_R=y_comps_R(tp);
		z_comp_R=z_comps_R(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_comp_R;y_comp_R;z_comp_R]),azd_R(valid_verts_R(v)),eld_R(valid_verts_R(v)));
		% convert now x-y equivalent vectors to magnitude
		xy=azelrho(1:2);
		magnitude = sqrt(sum(xy .^ 2));
		vwise_vec_R(v)=vwise_vec_R(v)+magnitude;
	end
    end
    % correct magnitude by length of time series to get average
    vwise_vec_R=vwise_vec_R./(size(m_vecs_R,2));
    % add in vwise vecs to m_L and m_R
    m_L=m_L+vwise_vec_L;
    m_R=m_R+vwise_vec_R;
end
% average across each scan
avg_L = m_L / length(fields);
avg_R = m_R / length(fields);
% saveout
if (~isnan(avg_L))
	save(strjoin([outFP '/AvgMag_m_L.mat'],""),'avg_L');
	save(strjoin([outFP '/AvgMag_m_R.mat'],""),'avg_R');
end
