addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% save out spins as giftis for downsampling
load('/oak/stanford/groups/leanew1/users/apines/DMNspins.mat');
% output in scratch
out_dir = '/scratch/users/apines/DMNspins_fs4_tmp';
mkdir(out_dir);
% write each permutation to a gifti
for p = 1:10000
    p
    % write left
    gL = gifti;
    gL.cdata = bigrotl(p,:)';
    save(gL, fullfile(out_dir, sprintf('perm%05d.L.func.gii', p)));

    % write right
    gR = gifti;
    gR.cdata = bigrotr(p,:)';
    save(gR, fullfile(out_dir, sprintf('perm%05d.R.func.gii', p)));
end
