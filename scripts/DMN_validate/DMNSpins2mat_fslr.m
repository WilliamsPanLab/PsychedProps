% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))

% convert downsampled DMN spins back to a .mat for easy loading in matlab
in_dir = '/scratch/users/apines/DMNspins_fslr_gii/';
n_perm = 10000;

% load first file as template
g = gifti(fullfile(in_dir, 'perm00001.L.32k.func.gii'));
n_vert = length(g.cdata);

% initialize
bigrotl_32k = zeros(n_perm, n_vert);
bigrotr_32k = zeros(n_perm, n_vert);
for p = 1:n_perm
    p
    % load left
    gL = gifti(fullfile(in_dir, sprintf('perm%05d.L.32k.func.gii', p)));
    bigrotl_32k(p, :) = gL.cdata(:)';
    % load right
    gR = gifti(fullfile(in_dir, sprintf('perm%05d.R.32k.func.gii', p)));
    bigrotr_32k(p, :) = gR.cdata(:)';
end
% save
save('/oak/stanford/groups/leanew1/users/apines/DMNspins_32k.mat', 'bigrotl_32k', 'bigrotr_32k', '-v7.3');

