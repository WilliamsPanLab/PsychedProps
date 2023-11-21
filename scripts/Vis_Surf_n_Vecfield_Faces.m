function Vis_Surf_n_Vecfield_Faces(surfl,surfr,vecl,vecr,Fn) 
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))
% designed for fsaverage4 surface
VertVecL=surfl;
VertVecR=surfr;
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
% Get incenters of triangles.
TR = TriRep(faces_l, vx_l);
P = TR.incenters;

% faces_R
F_R=faces_r;
% vertices V
V_R=vx_r;
TR = TriRep(faces_r, vx_r);
Pr = TR.incenters;

%%%%%%%%%%%%%%%%%%%%%%%%
% set data in the plotting script's terms
data=VertVecL;
ret=VecNormalize(vecl);
% replace nans with 0
ret(isnan(ret))=0;
% vector scaling factor
scalingfactor=1;
% colors
mincol=0;
maxcol=max(VertVecL(:));
% scale RGB values to max
minVal = min(VertVecL(:));
maxVal = max(VertVecL(:));
% consider scaling RGB colors by this value as well


%%% for red/blue 0-centered
%mincol=-9;
%maxcol=9;
%custommap=colormap(b2r(mincol,maxcol));
% abscense of color to gray to accom. lighting "none"
%custommap(126,:)=[.5 .5 .5];
%custommap=colormap(jet)

% blue-orange color scheme
%BO_cm=inferno(9);
%BO_cm(1,:)=[49 197 244];
%BO_cm(2,:)=[71 141 203];
%BO_cm(3,:)=[61 90 168];
%BO_cm(4,:)=[64 104 178];
%BO_cm(5,:)=[126 126 126];
%BO_cm(6,:)=[240 74 35];
%BO_cm(7,:)=[243 108 33];
%BO_cm(8,:)=[252 177 11];
%BO_cm(9,:)=[247 236 31];
% scale to 1
%BO_cm=BO_cm.*(1/255);
% interpolate color gradient
%interpsteps=[0 .125 .25 .375 .5 .625 .75 .875 1];
%BO_cm=interp1(interpsteps,BO_cm,linspace(0,1,255));
%custommap=BO_cm;
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
custommap=colormap(parula);
%custommap = [
%    1 0 0;    % Red
%    0 1 0;    % Green
%    0 0 1;    % Blue
%    1 1 0;    % Yellow
%    0 1 1;    % Cyan
%    1 0 1;    % Magenta
%    0 0 0;    % Black
%    1 1 1     % White
%];
% 0 to gray
%custommap(1,:)=[.5 .5 .5];
% medial left hemisphere
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/lh.sphere']);

% begin figure
figure
asub = subaxis(2,2,1, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
hold on;
% P from incenters up top
quiver3D(P(:,1),P(:,2),P(:,3),ret(:,1), ret(:,2), ret(:,3),[0 0 0],scalingfactor)
view([90 0]);
colormap(custommap)
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
shading flat;
camlight;

set(gca,'CLim',[mincol,maxcol]);
aplot.FaceVertexCData=VertVecL;
aplot.FaceAlpha=.7;

% other view of left hemisphere (lateral)
asub = subaxis(2,2,4, 'sh', 0.00, 'sv', 0.00, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
view([90 0]);
hold on;
colormap(custommap);
caxis([mincol; maxcol]);
bplot=quiver3D(P(:,1),P(:,2),P(:,3),ret(:,1), ret(:,2), ret(:,3),[0 0 0],scalingfactor)
daspect([1 1 1]);
axis tight;
rotate(aplot, [0 0 1], 180)
rotate(bplot, [0 0 1], 180)

axis vis3d off;
lighting none;
material metal %shiny %metal;
shading flat;
camlight;

% insert RGB colors onto surface
aplot.FaceVertexCData=VertVecL;
aplot.FaceAlpha=.7;

% right hemisphere
ret=VecNormalize(vecr);
% replace nans with 0
ret(isnan(ret))=0;
% vector scaling factor
scalingfactor=1;
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/rh.sphere']);
asub = subaxis(2,2,2, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
hold on;
% P from incenters up top
quiver3D(Pr(:,1),Pr(:,2),Pr(:,3),ret(:,1), ret(:,2), ret(:,3),[0 0 0],scalingfactor)
view([90 0]);
colormap(custommap)
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
shading flat;
camlight;

set(gca,'CLim',[mincol,maxcol]);
aplot.FaceVertexCData=VertVecR;
aplot.FaceAlpha=.7;

% other view of right hemisphere (lateral)
asub = subaxis(2,2,3, 'sh', 0.00, 'sv', 0.00, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
view([90 0]);
hold on;
colormap(custommap);
caxis([mincol; maxcol]);
bplot=quiver3D(Pr(:,1),Pr(:,2),Pr(:,3),ret(:,1), ret(:,2), ret(:,3),[0 0 0],scalingfactor)
daspect([1 1 1]);
axis tight;
rotate(aplot, [0 0 1], 180)
rotate(bplot, [0 0 1], 180)

axis vis3d off;
lighting none;
material metal %shiny %metal;
shading flat;
camlight;

% insert RGB colors onto surface
aplot.FaceVertexCData=VertVecR;
aplot.FaceAlpha=.7;


% printout
print(Fn,'-dpng','-r1000')





% [vertices, faces] = freesurfer_read_surf([SubjectsFolder '/rh.inflated']);

% asub = subaxis(2,2,2, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0,'Holdaxis',1);
% aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
% view([90 0]);
% rotate(aplot, [0 0 1], 180)
%colormap(custommap)
% caxis([mincol; maxcol]);
% daspect([1 1 1]);
% axis tight;
% axis vis3d off;
% lighting none
% material metal %shiny %metal;%shading flat;
% shading flat;
% camlight;
%  pos = get(asub, 'Position');
%  posnew = pos; posnew(1) = posnew(1) - 0.11; set(asub, 'Position', posnew);
% alpha(1)


% set(gca,'CLim',[mincol,maxcol]);
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');

% asub = subaxis(2,2,3, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);
% aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
% view([90 0]);
%colormap(custommap)
% caxis([mincol; maxcol]);
% daspect([1 1 1]);
% axis tight;
% axis vis3d off;
% lighting none;
% material metal %shiny %metal;
% shading flat;
% camlight;
% alpha(1)
% pos = get(asub, 'Position');
%  posnew = pos; posnew(2) = posnew(2) + 0.13; set(asub, 'Position', posnew);
% set(gcf,'Color','w')


% set(gca,'CLim',[mincol,maxcol]);
%%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');
% colorbar
% c=colorbar
% c.Location='southoutside'

% colormap(custommap)

% print(Fn,'-dpng')
