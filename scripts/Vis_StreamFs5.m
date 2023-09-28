function Vis_Vertvec(StreamMatL,StreamMatR,Fn) 

addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))

%%% Load in surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/fs5surf';
surfL = [SubjectsFolder '/lh.inflated'];
surfR = [SubjectsFolder '/rh.inflated'];
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
%%%% medial wall stuff
% use native freesurfer command for mw mask indices
surfML = [SubjectsFolder '/lh.Medial_wall.label'];
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = [SubjectsFolder '/rh.Medial_wall.label'];
mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "isn't medial wall" vector for vertices
mw_L=ones(1,10242);
mw_L(mwIndVec_l)=0;
mw_R=ones(1,10242);
mw_R(mwIndVec_r)=0;

% get number of instances where streamMat is connected: lets start with over 30
overThresh=sum(exceedingStreams)>30;
numExceeding=sum(overThresh);
exceedingInd=find(overThresh);



vertex1_index = 24
% extract this row
rowOfInt=StreamMatL(vertex1_index,:);
cellsOfInt=find(rowOfInt==100);

%% LOOP OVER ALL CELLSOFINT
vertex2_index = cellsOfInt(2);

% Get the coordinates of the two vertices
vertex1_coords = V_L(vertex1_index, :);
vertex2_coords = V_L(vertex2_index, :);

% Create the surface plot
figure;
aplot=trisurf(F_L, V_L(:, 1), V_L(:, 2), V_L(:, 3), 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none');

hold on;  % Enable "hold on" to keep the surface plot visible

% Plot a line connecting the two vertices
line_coords = [vertex1_coords; vertex2_coords];
plot3(line_coords(:, 1), line_coords(:, 2), line_coords(:, 3), 'r', 'LineWidth', 2);

view([90 0])
rotate(aplot, [0 0 1], 180)







%%%%%%%%%%%%%%%%%%%%%%%%
plotdata=zeros(1,10242);
plotdata(logical(mw_L))=VertVecL;
%plotdata=VertVecL;

%%% circular
mincol=0;
maxcol=8;

figure
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/lh.inflated']);
asub = subaxis(2,2,1, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);

aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),plotdata)
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
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');

asub = subaxis(2,2,4, 'sh', 0.00, 'sv', 0.00, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),plotdata)
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
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');


c=colorbar

colormap(custommap)

print(Fn,'-dpng')
