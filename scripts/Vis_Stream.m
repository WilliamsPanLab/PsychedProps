function Vis_Vertvec(StreamMatL,rowOfInt,Fn) 

addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))

%%% Load in surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
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
mw_L=ones(1,2562);
mw_L(mwIndVec_l)=0;
mw_R=ones(1,2562);
mw_R(mwIndVec_r)=0;
V_Lmasked=V_L(logical(mw_L'),:);  
nonMW_L=setdiff(1:2562,mwIndVec_l);

% get row of interest
OGrow=StreamMatL(rowOfInt,:);

% get next iteration of rows of interst
rowsOfInt=find(OGrow>200);

% initialize figure
figure
% Create an axis explicitly for the surface plot
ax = axes;
aplot=trisurf(F_L, V_L(:, 1), V_L(:, 2), V_L(:, 3), 'FaceColor', [0.5, 0.5, 0.5],'Parent', ax, 'EdgeColor', 'none','FaceAlpha',0.5);
hold on;

% get connections
for r=rowsOfInt;
	% set vertex index
	vertex1_index = r

	% extract this row
	Matrow=StreamMatL(vertex1_index,:);
	cellsOfInt=find(Matrow>100);

	%% LOOP OVER ALL CELLSOFINT
	for c=cellsOfInt;
		vertex2_index = c;
	
		% Get the coordinates of the two vertices
		vertex1_coords = V_Lmasked(vertex1_index, :);
		vertex2_coords = V_Lmasked(vertex2_index, :);

		% Plot a line connecting the two vertices
		line_coords = [vertex1_coords; vertex2_coords];

		alpha_value = Matrow(c) / max(StreamMatL(:)); % Scale alpha based on the maximum value in StreamMatL

		plot3(line_coords(:, 1), line_coords(:, 2), line_coords(:, 3), 'LineWidth', 1,'Color', [1 0 0 alpha_value]);
		hold on;

	end
	hold on;  % Enable "hold on" to keep the surface plot visible
end
axis vis3d off;
axis tight;
rotate(ax, [0 0 1], 180)
view(ax,[90 0])
daspect(ax,[1 1 1]);
print(['~/testStreams' num2str(rowOfInt) '.png'],'-dpng')
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
