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
bv_g_anglesx=zeros(17,5120);
p_g_anglesx=zeros(17,5120);
m1_g_anglesx=zeros(17,5120);
m2_g_anglesx=zeros(17,5120);
bv_g_anglesy=zeros(17,5120);
p_g_anglesy=zeros(17,5120);
m1_g_anglesy=zeros(17,5120);
m2_g_anglesy=zeros(17,5120);

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
	% save thetas a rho
	Thetas_L=zeros(1,5120);
	Rhos_L=zeros(1,5120);
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
			[Thetas_L(F),Rhos_L(F)]=cart2pol(vs_L(1),vs_L(2));
		end
	end
	% populate baseline dataframe
	[bv_g_anglesx(s,:),bv_g_anglesy(s,:)]=pol2cart(Thetas_L,Rhos_L);

	%%% p as left Hemi
	p=p.us.vf_left;
	NumTRs=size(p);
	NumTRs=NumTRs(2);
	lenOpFl=NumTRs;
	% initialize vector field to be extracted
	plotVF=zeros(n,3);
	% rhos and thetas
	Thetas_L=zeros(1,5120);
	Rhos_L=zeros(1,5120);
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
			[Thetas_L(F),Rhos_L(F)]=cart2pol(vs_L(1),vs_L(2));
		end
	end
	% populate placebo dataframe
        [p_g_anglesx(s,:),p_g_anglesy(s,:)]=pol2cart(Thetas_L,Rhos_L);

	%%% m1 as left Hemi
	m1=m1.us.vf_left;
	NumTRs=size(m1);
	NumTRs=NumTRs(2);
	lenOpFl=NumTRs;
	% initialize vector field to be extracted
	plotVF=zeros(n,3);
	% rhos and thetas
	Thetas_L=zeros(1,5120);
	Rhos_L=zeros(1,5120);
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
			[Thetas_L(F),Rhos_L(F)]=cart2pol(vs_L(1),vs_L(2));
		end
	end
	% populate m1 dataframe
        [m1_g_anglesx(s,:),m1_g_anglesy(s,:)]=pol2cart(Thetas_L,Rhos_L);

	%%% m2 as left Hemi
	m2=m2.us.vf_left;
	NumTRs=size(m2);
	NumTRs=NumTRs(2);
	lenOpFl=NumTRs;
	% initialize vector field to be extracted
	plotVF=zeros(n,3);
	% rhos and thetas
	Thetas_L=zeros(1,5120);
	Rhos_L=zeros(1,5120);
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
			[Thetas_L(F),Rhos_L(F)]=cart2pol(vs_L(1),vs_L(2));
		end
	end
	% populate m2 dataframe
        [m2_g_anglesx(s,:),m2_g_anglesy(s,:)]=pol2cart(Thetas_L,Rhos_L);
end
% remove 4th row
bv_g_anglesx(4,:)=[];
p_g_anglesx(4,:)=[];
m1_g_anglesx(4,:)=[];
m2_g_anglesx(4,:)=[];
bv_g_anglesy(4,:)=[];
p_g_anglesy(4,:)=[];
m1_g_anglesy(4,:)=[];
m2_g_anglesy(4,:)=[];
% get means across angles across each participant
bv_g_angles_meanx=mean(bv_g_anglesx);
p_g_angles_meanx=mean(p_g_anglesx);
m1_g_angles_meanx=mean(m1_g_anglesx);
m2_g_angles_meanx=mean(m2_g_anglesx);
% y components
bv_g_angles_meany=mean(bv_g_anglesy);
p_g_angles_meany=mean(p_g_anglesy);
m1_g_angles_meany=mean(m1_g_anglesy);
m2_g_angles_meany=mean(m2_g_anglesy);
% read in DMN as background
networks=load(['/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/group_Nets_fs4.mat']);
nets_LH=networks.nets.Lnets(:,2);
% convert to faces for equivalence to vector fields
F_nets_LH=sum(nets_LH(faces_l),2)./3;
vx_l=vx_l./101;
P_L=P_L./101;

% set rhos to 0: prior to spherical transform it represents vector coming out from sphere
rhos=zeros(1,5120);
% make x y rhos three rows in a matrix
vectors_s=zeros(3,5120);
vectors_s(1,:)=bv_g_angles_meanx;
vectors_s(2,:)=bv_g_angles_meany;
vectors_s(3,:)=rhos;
vectors_cart=zeros(3,5120);
% Streamlines for first component: transpose to meet format
for F=g_noMW_combined_L
        % get vectors in cartesian space
        vectors_cart(:,F)=sph2cartvec(vectors_s(:,F),azd_L(F),eld_L(F));
end
v=vectors_cart';
% Define your data
data = F_nets_LH;
% zero mw
data(fmwIndVec_l)=0;
% Define the colormap
cmap = colormap('autumn');
% Set values below 0.4 to gray
grayValue = [0.7 0.7 0.7]; % RGB values for gray color (adjust as needed)
data(data < 0.4) = 0;
cmap(1,:) = grayValue;

% set figure
figure('units','pixels','position',[0 0 2500 2500])
scalingfactor=5
% Adjust the color axis limits to match the 'autumn' colormap
caxis([0.4, max(data)]);
% add cortical mantle %
aplot = trisurf(faces_l, vx_l(:,1)./scalingfactor, vx_l(:,2)./scalingfactor, vx_l(:,3)./scalingfactor,data)
colormap(cmap);
hold on
% convert dmn to faces
set(aplot,'FaceColor','flat','FaceVertexCData',data,'CDataMapping','scaled');
freezeColors
scalingfactor=5;
ret=v;
% Adjust the color axis limits to match the 'autumn' colormap
caxis([0.4, max(nets_LH)]);
hold on;
% color vector!
% Calculate the magnitudes of the vectors in ret
magnitudes = sqrt(sum(ret.^2, 2));
% Normalize the magnitudes to the range [0, 1], .029 is max in placebo so we will use throughout for consistent coloration
normalized_magnitudes = (magnitudes - 0) / (.029 - min(magnitudes));
% Define the color range from [.1 .1 .1] to [.9 .9 .9]
color_range_low = [.9 .9 .9];
color_range_high = [.1 .1 .1];
% Interpolate the colors based on the normalized magnitudes
color_vector = color_range_low + (color_range_high - color_range_low) .* normalized_magnitudes;
% add quivers
quiver3D(P_L(g_noMW_combined_L,1)./scalingfactor,P_L(g_noMW_combined_L,2)./scalingfactor,P_L(g_noMW_combined_L,3)./scalingfactor,ret(g_noMW_combined_L,1), ret(g_noMW_combined_L,2), ret(g_noMW_combined_L,3),color_vector(g_noMW_combined_L),0)
% Set colormap
colormap('gray');
view(230, 185);
daspect([1 1 1]);

print(['~/streams/groupBV.png'],'-dpng')
%%%%%%%%%%%%%%
% placebo
% make x y rhos three rows in a matrix
vectors_s=zeros(3,5120);
vectors_s(1,:)=p_g_angles_meanx;
vectors_s(2,:)=p_g_angles_meany;
vectors_s(3,:)=rhos;
vectors_cart=zeros(3,5120);
% Streamlines for first component: transpose to meet format
for F=g_noMW_combined_L
        % get vectors in cartesian space
        vectors_cart(:,F)=sph2cartvec(vectors_s(:,F),azd_L(F),eld_L(F));
end
v=vectors_cart';
% set figure
figure('units','pixels','position',[0 0 2500 2500])
scalingfactor=5
% Adjust the color axis limits to match the 'autumn' colormap
caxis([0.4, max(data)]);

% add cortical mantle %
aplot = trisurf(faces_l, vx_l(:,1)./scalingfactor, vx_l(:,2)./scalingfactor, vx_l(:,3)./scalingfactor,data)
colormap(cmap);
hold on
% convert dmn to faces
set(aplot,'FaceColor','flat','FaceVertexCData',data,'CDataMapping','scaled');
scalingfactor=5;
freezeColors;
ret=v;
% Adjust the color axis limits to match the 'autumn' colormap
caxis([0.4, max(nets_LH)]);
hold on;
% color vector!
% Calculate the magnitudes of the vectors in ret
magnitudes = sqrt(sum(ret.^2, 2));
% Normalize the magnitudes to the range [0, 1], .029 is max in placebo so we will use throughout for consistent coloration
normalized_magnitudes = (magnitudes - 0) / (.029 - min(magnitudes));
% Interpolate the colors based on the normalized magnitudes
color_vector = color_range_low + (color_range_high - color_range_low) .* normalized_magnitudes;
% add quivers
quiver3D(P_L(g_noMW_combined_L,1)./scalingfactor,P_L(g_noMW_combined_L,2)./scalingfactor,P_L(g_noMW_combined_L,3)./scalingfactor,ret(g_noMW_combined_L,1), ret(g_noMW_combined_L,2), ret(g_noMW_combined_L,3),color_vector(g_noMW_combined_L),0)
% Set colormap
colormap('gray');
view(230, 185);
daspect([1 1 1]);

print(['~/streams/groupP.png'],'-dpng')
%%%%%%%%%%%%%%
% 80 mg
%%%%%%%%%%%%%%
% make x y rhos three rows in a matrix
vectors_s=zeros(3,5120);
vectors_s(1,:)=m1_g_angles_meanx;
vectors_s(2,:)=m1_g_angles_meany;
vectors_s(3,:)=rhos;
vectors_cart=zeros(3,5120);
% Streamlines for first component: transpose to meet format
for F=g_noMW_combined_L
        % get vectors in cartesian space
        vectors_cart(:,F)=sph2cartvec(vectors_s(:,F),azd_L(F),eld_L(F));
end
v=vectors_cart';
% set figure
figure('units','pixels','position',[0 0 2500 2500])
scalingfactor=5
% Adjust the color axis limits to match the 'autumn' colormap
caxis([0.4, max(data)]);
% add cortical mantle %
aplot = trisurf(faces_l, vx_l(:,1)./scalingfactor, vx_l(:,2)./scalingfactor, vx_l(:,3)./scalingfactor,data)
colormap(cmap);
hold on
% convert dmn to faces
set(aplot,'FaceColor','flat','FaceVertexCData',data,'CDataMapping','scaled');
freezeColors;
ret=v;
% Adjust the color axis limits to match the 'autumn' colormap
caxis([0.4, max(nets_LH)]);
hold on;
% color vector!
% Calculate the magnitudes of the vectors in ret
magnitudes = sqrt(sum(ret.^2, 2));
% Normalize the magnitudes to the range [0, 1], .029 is max here so we will use throughout for consistent coloration
normalized_magnitudes = (magnitudes - 0) / (.029 - min(magnitudes));
% Interpolate the colors based on the normalized magnitudes
color_vector = color_range_low + (color_range_high - color_range_low) .* normalized_magnitudes;
% add quivers
quiver3D(P_L(g_noMW_combined_L,1)./scalingfactor,P_L(g_noMW_combined_L,2)./scalingfactor,P_L(g_noMW_combined_L,3)./scalingfactor,ret(g_noMW_combined_L,1), ret(g_noMW_combined_L,2), ret(g_noMW_combined_L,3),color_vector(g_noMW_combined_L),0)
% Set colormap
colormap('gray');
view(230, 185);
daspect([1 1 1]);

print(['~/streams/groupm1.png'],'-dpng')
%%%%%%%%%%%%%%
% 120 mg
%%%%%%%%%%%%%%
% make x y rhos three rows in a matrix
vectors_s=zeros(3,5120);
vectors_s(1,:)=m2_g_angles_meanx;
vectors_s(2,:)=m2_g_angles_meany;
vectors_s(3,:)=rhos;
vectors_cart=zeros(3,5120);
% Streamlines for first component: transpose to meet format
for F=g_noMW_combined_L
        % get vectors in cartesian space
        vectors_cart(:,F)=sph2cartvec(vectors_s(:,F),azd_L(F),eld_L(F));
end
v=vectors_cart';
% set figure
figure('units','pixels','position',[0 0 2500 2500])
scalingfactor=5
% Adjust the color axis limits to match the 'autumn' colormap
caxis([0.4, max(data)]);

% add cortical mantle %
aplot = trisurf(faces_l, vx_l(:,1)./scalingfactor, vx_l(:,2)./scalingfactor, vx_l(:,3)./scalingfactor,data)
colormap(cmap);
hold on
% convert dmn to faces
set(aplot,'FaceColor','flat','FaceVertexCData',data,'CDataMapping','scaled');
ret=v;
hold on;
freezeColors;
% color vector!
% Calculate the magnitudes of the vectors in ret
magnitudes = sqrt(sum(ret.^2, 2));
% Normalize the magnitudes to the range [0, 1], .029 is max in placebo so we will use throughout for consistent coloration
normalized_magnitudes = (magnitudes - 0) / (.029 - min(magnitudes));
% Interpolate the colors based on the normalized magnitudes
color_vector = color_range_low + (color_range_high - color_range_low) .* normalized_magnitudes;
% add quivers
quiver3D(P_L(g_noMW_combined_L,1)./scalingfactor,P_L(g_noMW_combined_L,2)./scalingfactor,P_L(g_noMW_combined_L,3)./scalingfactor,ret(g_noMW_combined_L,1), ret(g_noMW_combined_L,2), ret(g_noMW_combined_L,3),color_vector(g_noMW_combined_L),0)
% Set colormap
colormap('gray');
view(230, 185);
daspect([1 1 1]);

print(['~/streams/groupm2.png'],'-dpng')

colorbar

print(['~/streams/groupm2cb.png'],'-dpng')
% make a drug no drug difference image 
% merge placebo and baseline vectors
vectors_nondrug=zeros(3,5120);
vectors_nondrug(1,:)=(bv_g_angles_meanx+p_g_angles_meanx)./2;
vectors_nondrug(2,:)=(bv_g_angles_meany+p_g_angles_meany)./2;
vectors_nondrug(3,:)=rhos;
% merge m1 and m2 vectors (80 and 120mg)
vectors_d=zeros(3,5120);
vectors_d(1,:)=(m1_g_angles_meanx+m2_g_angles_meanx)./2;
vectors_d(2,:)=(m1_g_angles_meany+m2_g_angles_meany)./2;
vectors_d(3,:)=rhos;
% get difference between merged nondrug and merged drug vectors 
vectors_diff=vectors_nondrug-vectors_d;
% Streamlines for first component: transpose to meet format
for F=g_noMW_combined_L
	% get vectors in cartesian space
	vectors_cart(:,F)=sph2cartvec(vectors_diff(:,F),azd_L(F),eld_L(F));
end
v=vectors_cart';
% set figure
figure('units','pixels','position',[0 0 2500 2500])
% Adjust the color axis limits to match the 'autumn' colormap
caxis([0.4, max(data)]);
% add cortical mantle %
aplot = trisurf(faces_l, vx_l(:,1)./scalingfactor, vx_l(:,2)./scalingfactor, vx_l(:,3)./scalingfactor,data)
colormap(cmap);
freezeColors;
hold on
% convert dmn to faces
set(aplot,'FaceColor','flat','FaceVertexCData',data,'CDataMapping','scaled');
scalingfactor=5;
% convert vectors to unit vectors for plotting independent of magnitude
% thank you https://stackoverflow.com/questions/40778629/matlab-convert-vector-to-unit-vector
ret=v;
% Adjust the color axis limits to match the 'autumn' colormap
caxis([0.4, max(nets_LH)]);
hold on;
% add quivers
quiver3D(P_L(g_noMW_combined_L,1)./scalingfactor,P_L(g_noMW_combined_L,2)./scalingfactor,P_L(g_noMW_combined_L,3)./scalingfactor,ret(g_noMW_combined_L,1), ret(g_noMW_combined_L,2), ret(g_noMW_combined_L,3),[.2 .2 .2],0)
% Set colormap
colormap(cmap);
freezeColors;
view(230, 185);
daspect([1 1 1]);
print(['~/streams/nond_mind.png'],'-dpng')





% make a sep. plot showing scale of arrows as legend?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

