function Vis_MWtsnr_mask(Fn) 

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
mw_L=fmwIndVec_l;
mw_R=fmwIndVec_r;

% black and white custommap colormap
% Define a custom colormap: white everywhere, black at value = 1
n = 256; % Number of colors in the colormap
custommap = ones(n, 3); % Initialize with white (R=1, G=1, B=1)

% Make the last row (corresponding to value=1) black
%custommap(1,:) = [.5, .5, .5];
custommap(end, :) = [0, 0, 0]; % Black (R=0, G=0, B=0)


figure
surfL = [SubjectsFolder '/lh.inflated'];
surfR = [SubjectsFolder '/rh.inflated'];
% +1 the faces: begins indexing at 0
[vertices, faces] = freesurfer_read_surf(surfL);
asub = subaxis(2,2,1, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0,'Holdaxis',1);

data=zeros(1,5120);
data(mw_L)=1;

aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
view([90 0]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
shading flat;
camlight;
	alpha(1)

% reset mincol here
mincol=0;
maxcol=1;
set(gca,'CLim',[mincol,maxcol]);
%set(aplot,'FaceColor','flat','FaceVertexCData',DMN_masked_L,'CDataMapping','scaled');
%set(aplot, 'FaceColor', [0.5, 0.5, 0.5], 'FaceVertexCData', []);
set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');

asub = subaxis(2,2,4, 'sh', 0.00, 'sv', 0.00, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
view([90 0]);
rotate(aplot, [0 0 1], 180)
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
data=zeros(1,5120);
data(mw_R)=1;

[vertices, faces] = freesurfer_read_surf(surfR);

asub = subaxis(2,2,2, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
view([90 0]);
rotate(aplot, [0 0 1], 180)
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

asub = subaxis(2,2,3, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
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
%set(aplot, 'FaceColor', [0.5, 0.5, 0.5], 'FaceVertexCData', []);

%c=colorbar;
%c=colorbar('XTickLabel',{'-5', '0', '5'},'XTick', -5:5:5)
%c.Location='southoutside'

print(Fn,'-dpng','-r2000')
