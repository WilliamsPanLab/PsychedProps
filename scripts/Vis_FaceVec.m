function Vis_FaceVec(FaceVecL,FaceVecR,Fn) 

addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))

%%% Load in surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfL = [SubjectsFolder '/lh.sphere'];
surfR = [SubjectsFolder '/rh.sphere'];
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
DMN_masked_L=sum(nets_LH(faces_l),2)./3;
DMN_masked_L(DMN_bool_L==0)=0;
DMN_masked_R=sum(nets_RH(faces_r),2)./3;
DMN_masked_R(DMN_bool_R==0)=0;
% combine with medial wall mask
MasterMask_L=DMN_bool_L;
MasterMask_R=DMN_bool_R;
MasterMask_L(fmwIndVec_l)=0;
MasterMask_R(fmwIndVec_r)=0;
% should be 1116 faces for left, 996 for right
mw_L=MasterMask_L;
mw_R=MasterMask_R;

% get gradient of DMN
%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/lukaslang-ofd-614a2ffc50d6'))
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/lh.inflated']);
F_L=faces;
V_L=vertices;
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/rh.inflated']);
F_R=faces;
V_R=vertices;
% get incenters of triangles
TR_L = TriRep(F_L,V_L);
P_L = TR_L.incenters;
TR_R = TriRep(F_R,V_R);
P_R = TR_R.incenters;
% calculate network gradients on sphere
ng_L = grad(F_L, V_L, nets_LH);
ng_R = grad(F_R, V_R, nets_RH);
% same masking procedure
ng_L(MasterMask_L)=0;
ng_R(MasterMask_R)=0;

% LEFT HEMISPHERE
% ORTHOGONALIZE TO SURFACE OF INFLATED
ret_L=ng_L;
% for each vector, subtract weighted surface-orthogonal component from original vector
for f=1:length(F_L)
        % retrieve original vector
        OGvec=ret_L(f,:);
        % find the three vertices involved in this face 
	verticesIdx = F_L(f, :);
        % get euclidean coords of involved vertices
	v1 = V_L(verticesIdx(1), :);
    	v2 = V_L(verticesIdx(2), :);
    	v3 = V_L(verticesIdx(3), :);
	% get normal vectors of each involved face
	normalVector = cross(v2 - v1, v3 - v1);
	% and norm it
	normalVector = normalVector / norm(normalVector);
        % get dot product of orthogonal vector and original vector
        OGvecOrthogonal = dot(OGvec, normalVector) * normalVector;
        % subtract orthogonal component of original vector from original vector
        modVec = OGvec - OGvecOrthogonal;;
        % add modified vector to initialized matrix
        ret_L(f,:)=modVec;
end
ret_L=VecNormalize(ret_L);
ret_L(~DMN_bool_L,:)=0;
data=zeros(1,5120);
%data(g_noMW_combined_L)=FaceVecL;
data(mw_L)=FaceVecL;
%%%%%%% fixed colorscale varities

% match r fig color scheme
color1 = [9 65 107];     % #09416b
color2 = [126 126 126];  % gray midpoint
color3 = [193 33 57];    % #c12139

% Normalize to [0,1]
color1 = color1 / 255;
color2 = color2 / 255;
color3 = color3 / 255;

% Interpolate across 255 steps using 3 anchor points
custommap = interp1([0 0.5 1], [color1; color2; color3], linspace(0, 1, 255));

figure
surfL = [SubjectsFolder '/lh.inflated'];
surfR = [SubjectsFolder '/rh.inflated'];
% +1 the faces: begins indexing at 0
[vertices, faces] = freesurfer_read_surf(surfL);
asub = subaxis(2,2,1, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0,'Holdaxis',1);

aplot = trisurf(F_L, V_L(:,1), V_L(:,2), V_L(:,3));
bplot=quiver3D(P_L(:,1),P_L(:,2),P_L(:,3),ret_L(:,1), ret_L(:,2), ret_L(:,3),data',1)
view([90 0]);
colormap(custommap)
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
shading flat;
camlight;
	alpha(1)

% reset mincol here
mincol=-5.1;
maxcol=5.1;
set(gca,'CLim',[mincol,maxcol]);
%set(aplot,'FaceColor','flat','FaceVertexCData',DMN_masked_L,'CDataMapping','scaled');
%set(aplot, 'FaceColor', [0.5, 0.5, 0.5], 'FaceVertexCData', []);
set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');
set(bplot, 'FaceColor', [.8, .8, .8]);

asub = subaxis(2,2,4, 'sh', 0.00, 'sv', 0.00, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
%bplot=quiver3D(P_L(:,1),P_L(:,2),P_L(:,3),ret_L(:,1), ret_L(:,2), ret_L(:,3),data',1)
bplot=quiver3D(P_L(:,1),P_L(:,2),P_L(:,3),ret_L(:,1), ret_L(:,2), ret_L(:,3),data',1)
view([90 0]);
rotate(aplot, [0 0 1], 180)
rotate(bplot, [0 0 1], 180)
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
material metal %shiny %metal;
shading flat;
camlight;
alpha(1)
 pos = get(asub, 'Position');
 posnew = pos; posnew(2) = posnew(2) + 0.13; posnew(1) = posnew(1) -.11; set(asub, 'Position', posnew);
set(gcf,'Color','w')

set(gca,'CLim',[mincol,maxcol]);
set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');
%set(aplot, 'FaceColor', [0.5, 0.5, 0.5], 'FaceVertexCData', []);
set(bplot, 'FaceColor', [.8, .8, .8]);

% RIGHT HEMISPHERE
% ORTHOGONALIZE TO SURFACE OF INFLATED
ret_R=ng_R;
% for each vector, subtract weighted surface-orthogonal component from original vector
for f=1:length(F_R)
        % retrieve original vector
        OGvec=ret_R(f,:);
        % find the three vertices involved in this face 
        verticesIdx = F_R(f, :);
        % get euclidean coords of involved vertices
        v1 = V_R(verticesIdx(1), :);
        v2 = V_R(verticesIdx(2), :);
        v3 = V_R(verticesIdx(3), :);
        % get normal vectors of each involved face
        normalVector = cross(v2 - v1, v3 - v1);
        % and norm it
        normalVector = normalVector / norm(normalVector);
        % get dot product of orthogonal vector and original vector
        OGvecOrthogonal = dot(OGvec, normalVector) * normalVector;
        % subtract orthogonal component of original vector from original vector
        modVec = OGvec - OGvecOrthogonal;;
        % add modified vector to initialized matrix
        ret_R(f,:)=modVec;
end
ret_R=VecNormalize(ret_R);
ret_R(~DMN_bool_R,:)=0;
data=zeros(1,5120);
data(mw_R)=FaceVecR;

[vertices, faces] = freesurfer_read_surf(surfR);

asub = subaxis(2,2,2, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
bplot=quiver3D(P_R(:,1),P_R(:,2),P_R(:,3),ret_R(:,1), ret_R(:,2), ret_R(:,3),data',1)
view([90 0]);
rotate(aplot, [0 0 1], 180)
rotate(bplot, [0 0 1], 180)
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none
material metal %shiny %metal;%shading flat;
shading flat;
camlight;
 pos = get(asub, 'Position');
 posnew = pos; posnew(1) = posnew(1) - 0.11; set(asub, 'Position', posnew);
alpha(1)


set(gca,'CLim',[mincol,maxcol]);
set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');
%set(aplot, 'FaceColor', [0.5, 0.5, 0.5], 'FaceVertexCData', []);
set(bplot, 'FaceColor', [.8, .8, .8]);

asub = subaxis(2,2,3, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
bplot=quiver3D(P_R(:,1),P_R(:,2),P_R(:,3),ret_R(:,1), ret_R(:,2), ret_R(:,3),data',1)
view([90 0]);
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
material metal %shiny %metal;
shading flat;
camlight;
alpha(1)
 pos = get(asub, 'Position');
 posnew = pos; posnew(2) = posnew(2) + 0.13; set(asub, 'Position', posnew);
set(gcf,'Color','w')


set(gca,'CLim',[mincol,maxcol]);
set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');
set(bplot, 'FaceColor', [.8, .8, .8]);
%set(aplot, 'FaceColor', [0.5, 0.5, 0.5], 'FaceVertexCData', []);

%c=colorbar;
%c=colorbar('XTickLabel',{'-5', '0', '5'},'XTick', -5:5:5)
%c.Location='southoutside'

print(Fn,'-dpng','-r2000')
