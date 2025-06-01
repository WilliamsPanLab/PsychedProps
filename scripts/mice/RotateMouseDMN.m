% load in mouse DMN
networks = load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
Dnet = networks.components_matrix';
Mask = networks.mask;

% re-boolean and transpose mask for consistency
Mask = Mask == 1;
Mask = Mask';

% treshold DMN mask (threshold > 0.6)
DMN_bool = Dnet > 0.6;
DMN_bool(Mask == 0) = 0;
DMN_bool = logical(DMN_bool);

% Apply masking rules
Dnet(~DMN_bool) = 0;        % Gray for non-DMN brain areas

% same rotation parameters to initiate
nPerms = 10000;                  
% need most of rotated DMN to fall within brainmask
minInsideFrac = 0.5;
permCount = 0;
attempts = 0;

% and initialize output as height width nperms
[H, W] = size(Dnet);
% padding for crop issues
padSize = 20;
paddedSize = [H + 2*padSize, W + 2*padSize];
DMN_perms = zeros(H, W, nPerms);

% loop over nperms
while permCount < nPerms
    attempts = attempts + 1;
    % random rotation angle between 10 and 350 degrees
    angle = randi([10, 350]);
    % build 2D rotation transform (around image center)
    tform_rot = affine2d([cosd(angle) -sind(angle) 0;
                          sind(angle)  cosd(angle) 0;
                          0            0           1]);
    % pad it with 0s so cropping doesnt make incompat. image sizes
    Dnet_padded = zeros(paddedSize);
    Dnet_padded(padSize+1:padSize+H, padSize+1:padSize+W) = Dnet;
    % eotate padded DMN
    DMN_rot_full = imwarp(Dnet_padded, tform_rot, 'FillValues', 0);
    % annoying but we need to center-crop the result back to original image size
    [H, W] = size(Dnet);
    [H2, W2] = size(DMN_rot_full);
    x_start = floor((W2 - W)/2) + 1;
    y_start = floor((H2 - H)/2) + 1;
    x_end = x_start + W - 1;
    y_end = y_start + H - 1;
    % crop to original size
    DMN_rot = DMN_rot_full(y_start:y_end, x_start:x_end);
    % random translation (at least 5 pixels in both directions)
    x_shift = sign(randn()) * (5 + randi(10));
    y_shift = sign(randn()) * (5 + randi(10));
    % translate rotated DMN
    DMN_shifted = imtranslate(DMN_rot, [x_shift y_shift], 'OutputView', 'same');
    % Check overlap: How much of the shifted DMN falls inside the brain mask?
    overlap = DMN_shifted & Mask;
    frac_inside = nnz(overlap) / nnz(DMN_shifted);
    % Only accept permutations with enough overlap
    if frac_inside > minInsideFrac
        permCount = permCount + 1
	DMN_perms(:, :, permCount) = DMN_shifted;
        % visualize first 10
	if permCount <11
		% build background image for visualization
        	baseImage = zeros(size(Mask));
        	baseImage(Mask == 1 & DMN_shifted == 0) = 0;  % Gray for non-DMN brain
		baseImage(DMN_shifted>0) = DMN_shifted(DMN_shifted>0);
        	baseImage(Mask == 0) = -0.1;   % White background
        	% Plot and save the permutation
        	f = figure('Visible', 'off');
        	imagesc(baseImage);
        	axis image off;
        	colormap([1 1 1; 0.5 0.5 0.5; hot]);  % white, gray, then hot colormap
        	caxis([-0.1 1]);
        	% include random rotation and translation info
        	title(sprintf('Perm %d | Rot %.1f°, Shift (%d,%d), %.1f%% inside mask', ...
        	      permCount, angle, x_shift, y_shift, frac_inside * 100));
        	saveas(f, sprintf('~/test_spin%d.png', permCount));
        	close(f);
    	end
    end
end
% save it out
save('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_permutations_10k.mat', 'DMN_perms', '-v7.3');
