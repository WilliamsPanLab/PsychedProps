function Vis_VertScorevec(subj,sesh) 
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))
% recreate VertVecScore FP from subj and sesh
VectorSscoresFP=['/scratch/users/apines/SimStreams/' subj '_' sesh '_vectorScores.mat'];
VectorSscores=load(VectorSscoresFP).MWvectorscores;
VertVecL=VectorSscores;
% output filename
Fn=['/scratch/users/apines/SimStreams/' subj '_' sesh '_vectorScores.png'];

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

%%%%%%%%%%%%%%%%%%%%%%%%
data=VertVecL;

% colors
mincol=0;
maxcol=max(VertVecL(:));

% scale RGB values to max
minVal = min(VertVecL(:));
maxVal = max(VertVecL(:));
% need everything to be above 0 (+absminval) but scaled 0-1 (./maxval+absminval)
RGBValues=(VertVecL+abs(minVal))./(maxVal+abs(minVal));



% load in number of TRs to scale vectors by
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh];
numTRs=0;
for task = ["rs1" "rs2" "emotion" "gambling" "wm"];
       task=char(task);
	CSIfp=[childfp '/' subj '_' sesh '_task-' task '_ValidSegments_Trunc.txt'];
	CSI = importdata(CSIfp);
	% assure that TR count is the same between time series and valid segments txt
	SegNum=size(CSI);
	SegNum=SegNum(1);
	% trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
	numTRsVS=CSI(SegNum,1)+CSI(SegNum,2)-1
	numTRs=numTRs+numTRsVS;
end
% scale vectors to be smaller if more TRs included
scalingfactor=4000/numTRs;
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
%custommap=colormap(parula);
% mw to black
%custommap(1,:)=[0 0 0];

custommap = [
    1 0 0;    % Red
    0 1 0;    % Green
    0 0 1;    % Blue
    1 1 0;    % Yellow
    0 1 1;    % Cyan
    1 0 1;    % Magenta
    0 0 0;    % Black
    1 1 1     % White
];


figure
% medial left hemisphere
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/lh.inflated']);
asub = subaxis(2,2,1, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);

aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
% add vector field with OG vector values
ret=VertVecL;
% * 1.02 to make vectors visible above  surface
bplot=quiver3D(vertices(:,1)*1.02,vertices(:,2)*1.02,vertices(:,3)*1.02,ret(:,1), ret(:,2), ret(:,3),[1 1 1],scalingfactor)

view([90 0]);
colormap(custommap)
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
shading flat;
camlight;

set(gca,'CLim',[mincol,maxcol]);
aplot.FaceVertexCData=RGBValues;
aplot.FaceAlpha=.8;

% other view of left hemisphere (lateral)
asub = subaxis(2,2,4, 'sh', 0.00, 'sv', 0.00, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
view([90 0]);
rotate(aplot, [0 0 1], 180)
%colormap(custommap)
caxis([mincol; maxcol]);
 pos = get(asub, 'Position');
 posnew = pos; posnew(2) = posnew(2) + 0.13; posnew(1) = posnew(1) -.11; set(asub, 'Position', posnew);
set(gcf,'Color','w')

set(gca,'CLim',[mincol,maxcol]);
aplot.FaceVertexCData=RGBValues;

% add vector field with OG vector values
ret=VertVecL;
% * 1.1 to make vectors visible above  surface
bplot=quiver3D(vertices(:,1)*1.02,vertices(:,2)*1.02,vertices(:,3)*1.02,ret(:,1), ret(:,2), ret(:,3),[1 1 1],scalingfactor)
rotate(bplot, [0 0 1], 180)

daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
material metal %shiny %metal;
shading flat;
camlight;

% insert RGB colors onto surface
aplot.FaceVertexCData=RGBValues;
aplot.FaceAlpha=.8;
% printout
print(Fn,'-dpng','-r2000')




%%% right hemisphere
% data=VertVecR;

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
