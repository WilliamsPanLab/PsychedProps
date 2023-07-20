% read in subject dosage correspondence, has to be before addpath for some silly reason
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');

% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfL = [SubjectsFolder '/lh.sphere'];
surfR = [SubjectsFolder '/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% get incenters of triangles
TR_L = TriRep(faces_l,vx_l);
P_L = TR_L.incenters;
TR_R = TriRep(faces_r,vx_r);
P_R = TR_R.incenters;

% use native freesurfer command for mw mask indices
surfML = [SubjectsFolder '/lh.Medial_wall.label'];
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = [SubjectsFolder '/rh.Medial_wall.label'];
mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "is medial wall" vector for vertices
mw_L=zeros(1,2562);
mw_L(mwIndVec_l)=1;
mw_R=zeros(1,2562);
mw_R(mwIndVec_r)=1;
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

% normalize verts to unit sphere
% left
numV=length(vx_l);
vx_l(numV+1:end, :) = VecNormalize(vx_l(numV+1:end, :));
% right
numV=length(vx_r);
vx_r(numV+1:end, :) = VecNormalize(vx_r(numV+1:end, :));

% calculate normative directionality streamlines for this subject and session
n = size(faces_l, 1);
T = TriRep(faces_l, vx_l);
P = T.incenters./100;

% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));

% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);


% initialize baseline, placebo, 80 mg, and 120 mg vectors
bv_g_angles=zeros(17,5120);
p_g_angles=zeros(17,5120);
m1_g_angles=zeros(17,5120);
m2_g_angles=zeros(17,5120);

% load each subject and get mean value within and across subjects at each vertex
% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')
for s=[1 2 3 5 6 7 8 9 10 11 12 13 14 15 16 17]
	s
	% get session info
        seshInfo=subSeshDose{s,2:5};
        %% ADD IN REST OF FILEPATHS
        bvFP=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_OpFl_rs.mat'];
        bvFP=strjoin(bvFP,'');
	pFP=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{2} '/' subjList(s) '_' seshInfo{2} '_OpFl_rs.mat'];
	pFP=strjoin(pFP,'');
	m1FP=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{3} '/' subjList(s) '_' seshInfo{3} '_OpFl_rs.mat'];
        m1FP=strjoin(m1FP,'');
	m2FP=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{4} '/' subjList(s) '_' seshInfo{4} '_OpFl_rs.mat'];
        m2FP=strjoin(m2FP,'');
	% load in optical flow output
	bv=load(bvFP);
	p=load(pFP);
	m1=load(m1FP);
	m2=load(m2FP);

	% calculate normative directionality streamlines for this subject and session
	n = size(faces_l, 1);
	T = TriRep(faces_l, vx_l);
	P = T.incenters./100;

	%%% bv as left Hemi
	bv=bv.us.vf_left;
	NumTRs=size(bv);
	NumTRs=NumTRs(2);
	lenOpFl=NumTRs;
	% initialize vector field to be extracted
	plotVF=zeros(n,3);
	% note, by not saving rho (just theta), we are discarding magnitude information at this point
	Thetas_L=zeros(1,5120);
	% loop over each face
	for F=g_noMW_combined_L
		% loop over each frame to extract from unfortunate cell structure
		for fr=1:lenOpFl
			% current vector field
	        	relVf_L=bv{fr};
			% xyz components
	        	xComp_L=relVf_L(F,1);
	        	yComp_L=relVf_L(F,2);
	        	zComp_L=relVf_L(F,3);
			% convert to spherical coord system
	        	vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(F),eld_L(F));
	       		% store in output vector (r is redundant across all vecs, only using az and el)
			[Thetas_L(F),rho]=cart2pol(vs_L(1),vs_L(2));
		end
	end
	% extract mean angle at each face
	SubjThetas=circ_mean(Thetas_L);
	% populate baseline dataframe
	bv_g_angles(s,:)=SubjThetas;

	%%% p as left Hemi
	p=p.us.vf_left;
	NumTRs=size(p);
	NumTRs=NumTRs(2);
	lenOpFl=NumTRs;
	% initialize vector field to be extracted
	plotVF=zeros(n,3);
	% note, by not saving rho (just theta), we are discarding magnitude information at this point
	Thetas_L=zeros(1,5120);
	% loop over each face
	for F=g_noMW_combined_L
		% loop over each frame to extract from unfortunate cell structure
		for fr=1:lenOpFl
			% current vector field
	        	relVf_L=p{fr};
			% xyz components
	        	xComp_L=relVf_L(F,1);
	        	yComp_L=relVf_L(F,2);
	        	zComp_L=relVf_L(F,3);
			% convert to spherical coord system
	        	vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(F),eld_L(F));
	       		% store in output vector (r is redundant across all vecs, only using az and el)
			[Thetas_L(F),rho]=cart2pol(vs_L(1),vs_L(2));
		end
	end
	% extract mean angle at each face
	SubjThetas=circ_mean(Thetas_L);
	% populate placebo dataframe
	p_g_angles(s,:)=SubjThetas;

	%%% m1 as left Hemi
	m1=m1.us.vf_left;
	NumTRs=size(m1);
	NumTRs=NumTRs(2);
	lenOpFl=NumTRs;
	% initialize vector field to be extracted
	plotVF=zeros(n,3);
	% note, by not saving rho (just theta), we are discarding magnitude information at this point
	Thetas_L=zeros(1,5120);
	% loop over each face
	for F=g_noMW_combined_L
		% loop over each frame to extract from unfortunate cell structure
		for fr=1:lenOpFl
			% current vector field
	        	relVf_L=m1{fr};
			% xyz components
	        	xComp_L=relVf_L(F,1);
	        	yComp_L=relVf_L(F,2);
	        	zComp_L=relVf_L(F,3);
			% convert to spherical coord system
	        	vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(F),eld_L(F));
	       		% store in output vector (r is redundant across all vecs, only using az and el)
			[Thetas_L(F),rho]=cart2pol(vs_L(1),vs_L(2));
		end
	end
	% extract mean angle at each face
	SubjThetas=circ_mean(Thetas_L);
	% populate 80mg dataframe
	m1_g_angles(s,:)=SubjThetas;

	%%% m2 as left Hemi
	m2=m2.us.vf_left;
	NumTRs=size(m2);
	NumTRs=NumTRs(2);
	lenOpFl=NumTRs;
	% initialize vector field to be extracted
	plotVF=zeros(n,3);
	% note, by not saving rho (just theta), we are discarding magnitude information at this point
	Thetas_L=zeros(1,5120);
	% loop over each face
	for F=g_noMW_combined_L
		% loop over each frame to extract from unfortunate cell structure
		for fr=1:lenOpFl
			% current vector field
	        	relVf_L=m2{fr};
			% xyz components
	        	xComp_L=relVf_L(F,1);
	        	yComp_L=relVf_L(F,2);
	        	zComp_L=relVf_L(F,3);
			% convert to spherical coord system
	        	vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(F),eld_L(F));
	       		% store in output vector (r is redundant across all vecs, only using az and el)
			[Thetas_L(F),rho]=cart2pol(vs_L(1),vs_L(2));
		end
	end
	% extract mean angle at each face
	SubjThetas=circ_mean(Thetas_L);
	% populate 120mg dataframe
	m2_g_angles(s,:)=SubjThetas;
end
% remove 4th row
bv_g_angles(4,:)=[];
p_g_angles(4,:)=[];
m1_g_angles(4,:)=[];
m2_g_angles(4,:)=[];
% get circular means across angles for each participant
bv_g_angles_mean=circ_mean(bv_g_angles);
p_g_angles_mean=circ_mean(p_g_angles);
m1_g_angles_mean=circ_mean(m1_g_angles);
m2_g_angles_mean=circ_mean(m2_g_angles);

% read in DMN as background
networks=load(['/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/group_Nets_fs4.mat']);
nets_LH=networks.nets.Lnets(:,2);

% set rhos to 1 to match unit sphere
rhos=ones(1,5120);
% convert thetas back to cartersian coordinates with pol2cart
[x,y]=pol2cart(bv_g_angles_mean,1);
% make x y rhos three rows in a matrix
vecs_l=zeros(3,5120);
vecs_l(1,:)=x;
vecs_l(2,:)=y;
vecs_l(3,:)=rhos;
vectors_cart=zeros(3,5120);
% for each face
for F=g_noMW_combined_L
	% get vectors in cartesian space
	vectors_cart(:,F)=sph2cartvec(vecs_l(:,F),azd_L(F),eld_L(F));
end
% Streamlines for first component: transpose to meet format
v = vectors_cart';
% scale cortical mantle
vx_l=vx_l./101;
P_L=P_L./101;
% move z coord of cortical mantle backwards
% vx_l(:,3)=vx_l(:,3)-1;


% Set parameters.
nmax = max(sqrt(sum((v).^2, 2)));
h = 0.1/nmax;
maxit = 50;
lw = .7;
% Define your data
data = nets_LH;
% Define the colormap
cmap = colormap('autumn');
% Set values below 0.4 to gray
grayValue = [0.3 0.3 0.3]; % RGB values for gray color (adjust as needed)
data(data < 0.4) = 0;
cmap(1,:) = grayValue;

% set figure
figure
% Plot the data using trisurf
aplot = trisurf(faces_l, vx_l(:,1), vx_l(:,2), vx_l(:,3), data);

% Adjust the color axis limits to match the 'autumn' colormap
caxis([0.4, max(nets_LH)]);

% add cortical mantle %
aplot = trisurf(faces_l, vx_l(:,1), vx_l(:,2), vx_l(:,3),data)
colormap(cmap);
freezeColors;
hold on


% Plot the data using trisurf
aplot = trisurf(faces_l, vx_l(:,1), vx_l(:,2), vx_l(:,3), data);
% Adjust the color axis limits to match the 'autumn' colormap
caxis([0.4, max(nets_LH)]);
% Add cortical mantle
hold on;
% Add quiver plot for the vector field
quiver3(P_L(g_noMW_combined_L, 1), P_L(g_noMW_combined_L, 2), P_L(g_noMW_combined_L, 3), v(g_noMW_combined_L, 1), v(g_noMW_combined_L, 2), v(g_noMW_combined_L, 3), 'Color', 'g');
% Set colormap
colormap(cmap);
freezeColors;
view(230, 185);
daspect([1 1 1]);

%h=get(gca,'Children');
%set(gca,'Children',[h(76009) h(1:76008)]);
%adjustFigure;
% one rotation for insula
%savefigure(F, fullfile(childfp, filename), '-png', '-r600');
print(['~/streams/groupBV.png'],'-dpng')
%%%%%%%%%%%%%%

