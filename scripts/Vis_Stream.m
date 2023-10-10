function Vis_Stream(subj,StreamMatL,rowOfInt) 
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
% load subject-specific curvature to map onto surface
subjCurvfp=['/scratch/users/apines/data/mdma/' subj '/' subj '_fs4space_curv.L.func.gii'];
curv=gifti(subjCurvfp);
% get row of interest
OGrow=StreamMatL(rowOfInt,:);

% get next iteration of rows of interst
rowsOfInt=find(OGrow>100);

% initialize figure
figure
% Create an axis explicitly for the surface plot
ax1 = subaxis(2,1,1, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);
aplot1=trisurf(F_L, V_L(:, 1), V_L(:, 2), V_L(:, 3),curv.cdata,'EdgeColor', 'none','FaceAlpha',0.5);
hold on;

highlightRadius = 5;
highlightVertex = rowOfInt;

% Get the coordinates of the highlighted vertex
highlightX = V_Lmasked(highlightVertex, 1);
highlightY = V_Lmasked(highlightVertex, 2);
highlightZ = V_Lmasked(highlightVertex, 3);
% rotate it for this side
rotation_angle = deg2rad(180);
rotated_vertex_coords = [cos(rotation_angle), -sin(rotation_angle), 0; sin(rotation_angle), cos(rotation_angle), 0; 0, 0, 1] * [highlightX highlightY highlightZ]';

% Create a bright point (sphere) at the specified vertex
scatter3(ax1, rotated_vertex_coords(1), rotated_vertex_coords(2), rotated_vertex_coords(3), highlightRadius, 'blue', 'filled');

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
	
		% Rotate the coordinates manually (180 degrees around Z-axis)
        	rotated_vertex1_coords = [cos(rotation_angle), -sin(rotation_angle), 0; sin(rotation_angle), cos(rotation_angle), 0; 0, 0, 1] * vertex1_coords';
	        rotated_vertex2_coords = [cos(rotation_angle), -sin(rotation_angle), 0; sin(rotation_angle), cos(rotation_angle), 0; 0, 0, 1] * vertex2_coords';
	
		% Plot a line connecting the two vertices
		line_coords = [rotated_vertex1_coords'; rotated_vertex2_coords'];
		alpha_value = Matrow(c) / max(StreamMatL(:)); % Scale alpha based on the maximum value in StreamMatL
		plot3(line_coords(:, 1), line_coords(:, 2), line_coords(:, 3), 'LineWidth', 1,'Color', [1 0 0 alpha_value]);
		hold on;

	end
	hold on;  % Enable "hold on" to keep the surface plot visible
end
axis vis3d off;
axis tight;
rotate(aplot1, [0 0 1], 180)
view(ax1,[90 0])
daspect(ax1,[1 1 1]);

% rotated plot
ax2 = subaxis(2,1,2, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);
% aplot2=trisurf(F_L, V_L(:, 1), V_L(:, 2), V_L(:, 3), 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none','FaceAlpha',0.5);
aplot2=trisurf(F_L, V_L(:, 1), V_L(:, 2), V_L(:, 3),curv.cdata,'EdgeColor', 'none','FaceAlpha',0.5);
hold on;

% highlight current assayed vertex
highlightVertex=rowOfInt;

% Get the coordinates of the highlighted vertex
highlightX = V_Lmasked(highlightVertex, 1);
highlightY = V_Lmasked(highlightVertex, 2);
highlightZ = V_Lmasked(highlightVertex, 3);

% Create a bright point (sphere) at the specified vertex
scatter3(ax2, highlightX, highlightY, highlightZ, highlightRadius, 'blue', 'filled');

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
rotate(aplot2,[0 0 1],0);
view(ax2,[90 0]);
daspect(ax2,[1 1 1]);
% save outi - CONSIDER SAVING TO SUBJECT DIRECTORY
rowChar=char(num2str(rowOfInt));
print(['~/SigStreams_' subj '_'  rowChar '.png'],'-dpng','-r600')
