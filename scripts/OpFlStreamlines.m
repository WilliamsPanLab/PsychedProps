function OpFlStreamlines(subj,sesh)

% load in optical flow output

%%% NEED D, T.INCENTERS, E

% calculate normative directionality streamlines for this subject and session
n = size(D.Faces, 1);
T = TriRep(D.Faces, D.Verts);
P = T.incenters;

% Set seed points.
[X, Y] = meshgrid(-1:0.05:1, -1:0.05:1);
idx = find(X.^2 + Y.^2 <= 1);
S = [X(idx), Y(idx)];

% Set parameters.
nmax = max(sqrt(sum((E{k}.U1 + E{k}.U2).^2, 2)));
h = 0.1/nmax;
maxit = 50;
lw = 1;

% Streamlines for first component.
v = E{k}.U1;

F = createFigure('summer', -1, 1, -1, 1);
streamlines2(P, v, S, h, maxit, 'summer', lw);
adjustFigure;
savefigure(F, fullfile(renderPath, 'streamu2', sprintf('%s-%i-%i-600dpi.png', filename, k, l)), '-png', '-r600');
savefigure(F, fullfile(renderPath, 'streamu2', sprintf('%s-%i-%i-1200dpi.png', filename, k, l)), '-png', '-r1200');
savefigure(F, fullfile(renderPath, 'streamu2', sprintf('%s-%i-%i-600dpi.jpg', filename, k, l)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile(renderPath, 'streamu2', sprintf('%s-%i-%i-1200dpi.jpg', filename, k, l)), '-jpg', '-r1200', '-q100');

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

