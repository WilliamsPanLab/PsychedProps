function Extract_DMNMag(subj,sesh,task)
% define some paths 
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(Paths{1}))

% set common filepath
commonFP=['/scratch/users/apines/data/mdma/'];
% load in optical flow data
fp=[commonFP '/' subj '/' sesh '/' subj '_' sesh '_' task '_OpFl.mat'];
OpFl=load(fp);
OpFl_L=OpFl.us.vf_left;
OpFl_R=OpFl.us.vf_right;

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
P_L = TR_L.incenters;
TR_R = TriRep(F_R,V_R);
P_R = TR_R.incenters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% make medial wall vector
g_noMW_combined_L=setdiff([1:5120],fmwIndVec_l);
g_noMW_combined_R=setdiff([1:5120],fmwIndVec_r);

% get length of time series
lengthOpFl=size(OpFl_L,2);

% load in Networks
networks=load(['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/gro_Nets_fs4.mat']);
%% k = 4 to select FPN.
Dnet_LH=networks.nets.Lnets(:,4);
Dnet_RH=networks.nets.Rnets(:,4);
% create face-wise network mask
DMN_bool_L=sum(Dnet_LH(faces_l),2)./3;
DMN_bool_R=sum(Dnet_RH(faces_r),2)./3;
DMN_bool_L(DMN_bool_L>.3)=1;
DMN_bool_R(DMN_bool_R>.3)=1;
DMN_bool_L(DMN_bool_L<.3)=0;
DMN_bool_R(DMN_bool_R<.3)=0;
DMN_bool_L=logical(DMN_bool_L);
DMN_bool_R=logical(DMN_bool_R);

% ensure MW and SNR not represented in DMN mask
DMN_bool_L(F_MW_L==1)=0;
DMN_bool_R(F_MW_R==1)=0;

% for indexing within the forloops
DMN_bool_L_inds=find(DMN_bool_L);
DMN_bool_R_inds=find(DMN_bool_R);

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

% for every face, get angles in 2d instead of 3d (3 dimensions is redudnant for spherical surface)
% LEFT
for f=1:sum(DMN_bool_L);
	% get face index
	FInd=DMN_bool_L_inds(f);
	% for each timepoint
	for tp=1:lengthOpFl
		relVf_L=OpFl_L{tp};
		% xyz components
                xComp_L=relVf_L(f,1);
                yComp_L=relVf_L(f,2);
                zComp_L=relVf_L(f,3);	
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(DMN_bool_L_inds(f)),eld_L(DMN_bool_L_inds(f)));
		xy=azelrho(1:2);
		% calculate magnitude!
		mags_L(f,tp)=sqrt(sum(xy .^ 2));
	end
end
% RIGHT
for f=1:sum(DMN_bool_R);
        % get face index
        FInd=DMN_bool_R_inds(f);
        % for each timepoint
        for tp=1:lengthOpFl
		relVf_R=OpFl_R{tp};
                % xyz components
                xComp_R=relVf_R(f,1);
                yComp_R=relVf_R(f,2);
                zComp_R=relVf_R(f,3); 
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([xComp_R;yComp_R;zComp_R]),azd_R(DMN_bool_R_inds(f)),eld_R(DMN_bool_R_inds(f)));
                xy=azelrho(1:2);
                % calculate magnitude!
                mags_R(f,tp)=sqrt(sum(xy .^ 2));
        end     
end

% get average vector magnitude for each vertex
% LEFT
leftAvgMag=mean(mags_L,2);
% RIGHT
rightAvgMag=mean(mags_R,2);
% combined for averaging with weight-by-vertex rather than by-hemisphere
Mags=[leftAvgMag; rightAvgMag];
% get average across all DMN vertices
avgMag=mean(Mags);
% keep full facewise vectors as tables
stringVecAng_L=compose("Face%d", 1:size(leftAvgMag,1));
stringVecAng_R=compose("Face%d", 1:size(rightAvgMag,1));
leftAvgMag=table(leftAvgMag,'RowNames',stringVecAng_L);
rightAvgMag=table(rightAvgMag,'RowNames',stringVecAng_R);
% save out 
T=table(avgMag,'RowNames',"Row1");
outFP=['/scratch/users/apines/data/mdma/' subj '/' sesh];
writetable(T,[outFP '/' subj '_' sesh '_' task '_FPNMag.csv'],'WriteRowNames',true)
% save out facewise
writetable(leftAvgMag,[outFP '/' subj '_' sesh '_' task '_FPNMags_L.csv'],'WriteRowNames',true)
writetable(rightAvgMag,[outFP '/' subj '_' sesh '_' task '_FPNMags_R.csv'],'WriteRowNames',true)
