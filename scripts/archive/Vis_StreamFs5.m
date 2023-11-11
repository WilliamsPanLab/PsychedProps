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
V_Lmasked=V_L(logical(mw_L'),:);  
% get number of instances where streamMat is connected: lets start with over 30
overThresh=sum(exceedingStreams)>25;
numExceeding=sum(overThresh);
exceedingInd=find(overThresh);

for v=exceedingInd
% set vertex index
vertex1_index = v

% extract this row
rowOfInt=StreamMatL(vertex1_index,:);
cellsOfInt=find(rowOfInt>1);

% init figure
figure
aplot=trisurf(F_L, V_L(:, 1), V_L(:, 2), V_L(:, 3), 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none','FaceAlpha',0.5);

hold on;  % Enable "hold on" to keep the surface plot visible
%% LOOP OVER ALL CELLSOFINT
for c=cellsOfInt(1:30);
vertex2_index = c;

% Get the coordinates of the two vertices
vertex1_coords = V_Lmasked(vertex1_index, :);
vertex2_coords = V_Lmasked(vertex2_index, :);

% Plot a line connecting the two vertices
line_coords = [vertex1_coords; vertex2_coords]
plot3(line_coords(:, 1), line_coords(:, 2), line_coords(:, 3), 'r', 'LineWidth', 2);

hold on;  % Enable "hold on" to keep the surface plot visible

end
axis tight;
axis vis3d off;
view([90 0])
daspect([1 1 1]);
rotate(aplot, [0 0 1], 180)
print(['~/testStreams' num2str(v) '.png'],'-dpng')
end

% print some cases where it is below thresh





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
