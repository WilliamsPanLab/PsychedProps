function Vis_Vertvec(VertVecL,VertVecR,Fn) 

addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))

%%% Load in surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfL = [SubjectsFolder '/lh.sphere'];
surfR = [SubjectsFolder '/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
%faces_l = faces_l + 1;
%faces_r = faces_r + 1;
% faces_L
%F_L=faces_l;
% vertices V
%V_L=vx_l;
% faces_R
%F_R=faces_r;
% vertices V
%V_R=vx_r;

% use native freesurfer command for mw mask indices
%surfML = [SubjectsFolder '/lh.Medial_wall.label'];
%mwIndVec_l = read_medial_wall_label(surfML);
%surfMR = [SubjectsFolder '/rh.Medial_wall.label'];
%mwIndVec_r = read_medial_wall_label(surfMR);
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
mw_L=ones(1,2562);
mw_L(mwAndTSNR_L>0)=0;
mw_R=ones(1,2562);
mw_R(mwAndTSNR_R>0)=0;
mwIndVec_l=find(mw_L);
mwIndVec_r=find(mw_R);


%%%%%%%%%%%%%%%%%%%%%%%%
data=zeros(1,2562);
data(mwIndVec_l)=VertVecL;
%data(mwIndVec_l)=0;
% for medial wall plotting
%data=VertVecL;

%%%%%%% fixed colorscale varities

% blue-orange color scheme
BO_cm=inferno(9);
BO_cm(1,:)=[49 197 244];
BO_cm(2,:)=[71 141 203];
BO_cm(3,:)=[61 90 168];
BO_cm(4,:)=[64 104 178];
BO_cm(5,:)=[126 126 126];
BO_cm(6,:)=[240 74 35];
BO_cm(7,:)=[243 108 33];
BO_cm(8,:)=[252 177 11];
BO_cm(9,:)=[247 236 31];
% scale to 1
BO_cm=BO_cm.*(1/255);
% interpolate color gradient
interpsteps=[0 .125 .25 .375 .5 .625 .75 .875 1];
BO_cm=interp1(interpsteps,BO_cm,linspace(0,1,255));
custommap=BO_cm;
%custommap=colormap(inferno);

%%% matches circular hist
% for 180 degree max
%roybigbl_cm=inferno(6);
%roybigbl_cm(1,:)=[0, 0, 255];
%roybigbl_cm(2,:)=[0, 255, 255];
%roybigbl_cm(3,:)=[116, 192, 68];
%roybigbl_cm(4,:)=[246, 235, 20];
%roybigbl_cm(5,:)=[255, 165, 0];
%roybigbl_cm(6,:)=[255, 0, 0];
% scale to 1
%roybigbl_cm=roybigbl_cm.*(1/255);
% interpolate color gradient
%interpsteps=[0 .2 .4 .6 .8 1];
%roybigbl_cm=interp1(interpsteps,roybigbl_cm,linspace(0,1,255));
% make circular with flipud
%custommap=vertcat(flipud(roybigbl_cm),roybigbl_cm);

% for 90 degree max
%roybigbl_cm=inferno(3);
% blue
%roybigbl_cm(1,:)=[0, 0, 255];
% cyan
%roybigbl_cm(2,:)=[0, 255, 255];
% green
%roybigbl_cm(3,:)=[116, 192, 68];
% yellow
%roybigbl_cm(4,:)=[246, 235, 20];
% scale to 1
%roybigbl_cm=roybigbl_cm.*(1/255);
% interpolate color gradient
%interpsteps=[0 .33333 .66666 1];
%interpsteps=[0 .5 1];
%roybigbl_cm=interp1(interpsteps,roybigbl_cm,linspace(0,1,255));
% make circular with flipud
%custommap=vertcat(flipud(roybigbl_cm),roybigbl_cm);
%custommap=roybigbl_cm;
%custommap=colormap(parula);
% mw to black
%custommap(1,:)=[0 0 0];

% load in left cortical surface
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/lh.inflated']);
F_L=faces;
V_L=vertices;
% LOAD IN DMN STUFF
network=load(['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/Smooth_Nets_fs4.mat']);
n_LH=network.nets.Lnets(:,1);
n_RH=network.nets.Rnets(:,1);

% use same DMN mask used in analysis
n_LHBool=ones(1,2562);
n_LHBool(n_LH<.3)=0;
data(~logical(n_LHBool))=0;
% reset mincol here
mincol=-.012;
maxcol=.012;
figure
asub = subaxis(2,2,1, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);

aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
view([90 0]);
colormap(custommap)
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
shading flat;
camlight;
	alpha(1)

set(gca,'CLim',[mincol,maxcol]);
%bplot=quiver3D(vertices(VertVecL>0,1),vertices(VertVecL>0,2),vertices(VertVecL>0,3),ret(VertVecL>0,1), ret(VertVecL>0,2), ret(VertVecL>0,3),[.3 .5 .7])
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');

asub = subaxis(2,2,4, 'sh', 0.00, 'sv', 0.00, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
%bplot=quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),[.3 .5 .7])
view([90 0]);
rotate(aplot, [0 0 1], 180)
%rotate(bplot, [0 0 1], 180)
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
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');


%%% right hemisphere
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/rh.inflated']);
V_R=vertices;
F_R=faces;

data=zeros(1,2562);
data(mwIndVec_r)=VertVecR;
% for medial wall plotting
%data=VertVecR;

% use same DMN mask used in analysis
n_RHBool=ones(1,2562);
n_RHBool(n_RH<.3)=0;
data(~logical(n_RHBool))=0;

asub = subaxis(2,2,2, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0,'Holdaxis',1);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
%bplot = quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),[.3 .5 .7])
view([90 0]);
rotate(aplot, [0 0 1], 180)
%rotate(bplot, [0 0 1], 180)
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
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');

asub = subaxis(2,2,3, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
%bplot=quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),[.3 .5 .7])
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
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');
%colorbar
c=colorbar
c.Location='southoutside'

colormap(custommap)
print(Fn,'-dpng','-r800')
