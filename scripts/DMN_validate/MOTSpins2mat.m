% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))

% convert downsampled DMN spins back to a .mat for easy loading in matlab
in_dir = '/scratch/users/apines/MOTspins_fs4_gii/';
n_perm = 10000;

% load first file as template
g = gifti(fullfile(in_dir, 'perm00001.L.3k.func.gii'));
n_vert = length(g.cdata);

% initialize
bigrotl_3k = zeros(n_perm, n_vert);
bigrotr_3k = zeros(n_perm, n_vert);
for p = 1:n_perm
    p
    % load left
    gL = gifti(fullfile(in_dir, sprintf('perm%05d.L.3k.func.gii', p)));
    bigrotl_3k(p, :) = gL.cdata(:)';
    % load right
    gR = gifti(fullfile(in_dir, sprintf('perm%05d.R.3k.func.gii', p)));
    bigrotr_3k(p, :) = gR.cdata(:)';
end
% save
save('/oak/stanford/groups/leanew1/users/apines/MOTspins_3k.mat', 'bigrotl_3k', 'bigrotr_3k', '-v7.3');

