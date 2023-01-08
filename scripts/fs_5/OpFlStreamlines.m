function OpFlStreamlines(subj,sesh,filename)

% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
% load in optical flow output
data=load([childfp '/' subj '_' sesh '_OpFl_rs_fs5.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/fs5surf';
% for surface data
surfL = [SubjectsFolder '/lh.sphere'];
surfR = [SubjectsFolder '/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
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

%%% E as left Hemi
E=data.us.vf_left;
% initialize vector field to be plotted
plotVF=zeros(n,3);
% loop over to extract from unfortunate cell structure
for i=1:length(E);
	plotVF=plotVF+E{i};
end
% start with mean of each x y and z component
plotVF=plotVF./length(E);

% consider circular mean: might require cart2sphvec

% Set seed points.
[X, Y] = meshgrid(-1:0.04:1, -1:0.04:1);
idx = find(X.^2 + Y.^2 <= 1);
S = [X(idx), Y(idx)];

% Set parameters.
nmax = max(sqrt(sum((plotVF).^2, 2)));
h = 0.1/nmax;
maxit = 70;
lw = 1;

% Streamlines for first component.
v = plotVF;

F = createFigure('summer', -1, 1, -1, 1);
% read in curvature as data for background surf
curv=read_curv('/share/software/user/open/freesurfer/6.0.0/subjects/fsaverage5/surf/lh.curv']);
% scale cortical mantle
vx_l=vx_l./100;
% add cortical mantle %
aplot = trisurf(faces_l, vx_l(:,1), vx_l(:,2), vx_l(:,3),curv)
streamlines2(P, v, S, h, maxit, 'summer', lw);
%adjustFigure;
% one rotation for insula
savefigure(F, fullfile(childfp, filename), '-png', '-r600');
%%%%%%%%%%%%%%







F = createFigure3('summer');
% Create white sphere so that manifold is not transparent.
[x,y,z] = sphere;
idx = all(z >= 0, 2);
surf(x(idx, :), y(idx, :), z(idx, :), 'FaceColor', 'white', 'EdgeColor', 'white');
% Plot streamlines.
streamlines3(P, v, S, h, maxit, 'summer', lw);
set(gca, 'ZLim', [0, 1]);
set(gca, 'XLim', [-1, 1]);
set(gca, 'YLim', [-1, 1]);
adjustFigure3;
view(3);
savefigure(F, fullfile(renderPath, 'streamu3', sprintf('%s-%i-%i-600dpi.png', filename, k, l)), '-png', '-r600');
savefigure(F, fullfile(renderPath, 'streamu3', sprintf('%s-%i-%i-1200dpi.png', filename, k, l)), '-png', '-r1200');
savefigure(F, fullfile(renderPath, 'streamu3', sprintf('%s-%i-%i-600dpi.jpg', filename, k, l)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile(renderPath, 'streamu3', sprintf('%s-%i-%i-1200dpi.jpg', filename, k, l)), '-jpg', '-r1200', '-q100');

% Rotate by pi.
[az, el] = view;
view(az + 180, el);
savefigure(F, fullfile(renderPath, 'streamu3', sprintf('%s-%i-%i-rotated-600dpi.png', filename, k, l)), '-png', '-r600');
savefigure(F, fullfile(renderPath, 'streamu3', sprintf('%s-%i-%i-rotated-1200dpi.png', filename, k, l)), '-png', '-r1200');
savefigure(F, fullfile(renderPath, 'streamu3', sprintf('%s-%i-%i-rotated-600dpi.jpg', filename, k, l)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile(renderPath, 'streamu3', sprintf('%s-%i-%i-rotated-1200dpi.jpg', filename, k, l)), '-jpg', '-r1200', '-q100');

