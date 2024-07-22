function Extract_DMNMag(subj,sesh,task)
% define some paths 
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(Paths{1}))

% set common filepath
commonFP=['/scratch/users/apines/data/psil/'];

% load in optical flow data
fpl=[commonFP '/' subj '/' sesh '/' subj '_' sesh '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
fpr=[commonFP '/' subj '/' sesh '/' subj '_' sesh '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
OpFl_L=load(fpl);
OpFl_R=load(fpr);
% get length of time series
lengthOpFl=size(OpFl_L.vertWise_Vecs_l,2);

% load in Networks
networks=load(['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/gro_Nets_fs4.mat']);
%% k = 1 to select DMN.
Dnet_LH=networks.nets.Lnets(:,1);
Dnet_RH=networks.nets.Rnets(:,1);

% load in medial wall + SNR so we don't have to loop over every single vertex and then mask out later
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
% note this is indexing for VALID vertices as opposed to some other scripts with 1 at INVALID vertices
mw_L=ones(1,2562);
mw_L(mwAndTSNR_L==0)=0;
mw_R=ones(1,2562);
mw_R(mwAndTSNR_R==0)=0;
mw_L=logical(mw_L);
mw_R=logical(mw_R);

% get DMN indices (>.3)
DMN_bool_L=zeros(1,2562);
DMN_bool_R=zeros(1,2562);
DMN_bool_L(Dnet_LH>.3)=1;
DMN_bool_R(Dnet_RH>.3)=1;

% ensure MW and SNR not represented in DMN mask
DMN_bool_L(mw_L==1)=0;
DMN_bool_R(mw_R==1)=0;

% for indexing within the forloops
DMN_bool_L_inds=find(DMN_bool_L);
DMN_bool_R_inds=find(DMN_bool_R);

%%%% andddd we'll need template surface stuff for calculating angle transformations
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


% initialize output vectors: note these are in TANGENT PLANE TO SPHERICAL VERTEX not in original xyz space
%xy_comps_L=zeros(sum(DMN_bool_L),lengthOpFl,2);
%xy_comps_R=zeros(sum(DMN_bool_R),lengthOpFl,2);
% initialize vertexwise output magnitudes
mags_L=zeros(sum(DMN_bool_L),lengthOpFl);
mags_R=zeros(sum(DMN_bool_R),lengthOpFl);

% for every vertex, get angles in 2d instead of 3d (3 dimensions is redudnant for spherical surface)
% LEFT
for v=1:sum(DMN_bool_L);
	% get vertex index
	vertInd=DMN_bool_L_inds(v);
	% for each timepoint
	for tp=1:lengthOpFl
		x_compOG=OpFl_L.vertWise_Vecs_l(vertInd,tp,1);
		y_compOG=OpFl_L.vertWise_Vecs_l(vertInd,tp,2);
		z_compOG=OpFl_L.vertWise_Vecs_l(vertInd,tp,3);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([x_compOG;y_compOG;z_compOG]),azd_L(DMN_bool_L_inds(v)),eld_L(DMN_bool_L_inds(v)));
		xy=azelrho(1:2);
		% calculate magnitude!
		mags_L(v,tp)=sqrt(sum(xy .^ 2));
	end
end
% RIGHT
for v=1:sum(DMN_bool_R);
        % get vertex index
        vertInd=DMN_bool_R_inds(v);
        % for each timepoint
        for tp=1:lengthOpFl
                x_compOG=OpFl_R.vertWise_Vecs_r(vertInd,tp,1);
                y_compOG=OpFl_R.vertWise_Vecs_r(vertInd,tp,2);
                z_compOG=OpFl_R.vertWise_Vecs_r(vertInd,tp,3);
                % convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
                azelrho=cart2sphvec(double([x_compOG;y_compOG;z_compOG]),azd_R(DMN_bool_R_inds(v)),eld_R(DMN_bool_R_inds(v)));
                xy=azelrho(1:2);
                % calculate magnitude!
                mags_R(v,tp)=sqrt(sum(xy .^ 2));
        end     
end

% read in angular distance from extract_relangleverts script
inFP=['/scratch/users/apines/data/psil/' subj '/' sesh];
AngDistFP=[inFP '/' subj '_' sesh '_' task '_k1_AngDistVertsMat.mat'];
AngDist=load(AngDistFP).AngDist;
% make BUP mask for both hemis
BUP_L=AngDist.Left<90;
BUP_R=AngDist.Right<90;
% make TD mask for both hemis
TD_L=AngDist.Left>90;
TD_R=AngDist.Right>90;

% get average vector magnitude for each vertex
% LEFT
leftAvgMag=mean(mags_L,2);
% RIGHT
rightAvgMag=mean(mags_R,2);
% combined for averaging with weight-by-vertex rather than by-hemisphere
Mags=[leftAvgMag; rightAvgMag];
% get average across all DMN vertices
avgMag=mean(Mags);
%%% now for BUP
leftAvgMag=mean(mags_L(BUP_L),2);
% RIGHT
rightAvgMag=mean(mags_R(BUP_R),2);
% combined for averaging with weight-by-vertex rather than by-hemisphere
Mags=[leftAvgMag; rightAvgMag];
% get average across all DMN vertices
avgMag_bup=mean(Mags);
%%% now for TD
leftAvgMag=mean(mags_L(TD_L),2);
% RIGHT
rightAvgMag=mean(mags_R(TD_R),2);
% combined for averaging with weight-by-vertex rather than by-hemisphere
Mags=[leftAvgMag; rightAvgMag];
% get average across all DMN vertices 
avgMag_td=mean(Mags);

% save out
T=table(avgMag,'RowNames',"Row1");
outFP=['/scratch/users/apines/data/psil/' subj '/' sesh];
writetable(T,[outFP '/' subj '_' sesh '_' task '_DMNMag.csv'],'WriteRowNames',true)

% save out
T=table(avgMag_bup,'RowNames',"Row1");
outFP=['/scratch/users/apines/data/psil/' subj '/' sesh];
writetable(T,[outFP '/' subj '_' sesh '_' task '_DMNMag_BUP.csv'],'WriteRowNames',true)

% save out
T=table(avgMag_td,'RowNames',"Row1");
outFP=['/scratch/users/apines/data/psil/' subj '/' sesh];
writetable(T,[outFP '/' subj '_' sesh '_' task '_DMNMag_TD.csv'],'WriteRowNames',true)
