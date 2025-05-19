% 📂load DMN components and mask
networks = load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
Dnet = networks.components_matrix';
Mask = networks.mask;

% ✅re-boolean and transpose mask for consistency
Mask = Mask == 1;
Mask = Mask';

% 🧪hreshold DMN mask (threshold > 0.6)
DMN_bool = Dnet > 0.6;
DMN_bool(Mask == 0) = 0;
DMN_bool = logical(DMN_bool);

% Apply masking rules
Dnet(~DMN_bool) = 0;        % Gray for non-DMN brain areas
Dnet(Mask == 0) = -0.1;     % White for non-brain areas

% Upsample Dnet and mask by a factor of 2
upsample_factor = 2;
Dnet_upsampled = imresize(Dnet, upsample_factor, 'bilinear');
Mask_upsampled = imresize(Mask, upsample_factor, 'nearest');  % Preserve binary mask


% set colormap

% Custom colormap as specified
roybigbl_cm = inferno(16);
roybigbl_cm(1,:) = [255, 255, 0];
roybigbl_cm(2,:) = [255, 200, 0];
roybigbl_cm(3,:) = [255, 120, 0];
roybigbl_cm(4,:) = [255, 0, 0];
roybigbl_cm(5,:) = [200, 0, 0];
roybigbl_cm(6,:) = [150, 0, 0];
roybigbl_cm(7,:) = [200, 200, 200];
roybigbl_cm(8,:) = [200, 200, 200]; % Gray tone for masked areas
roybigbl_cm(9,:) = [200, 200, 200]; % Gray tone
roybigbl_cm(10,:) = [200, 200, 200]; % Gray tone
roybigbl_cm(11,:) = [200, 200, 200]; % Gray tone
roybigbl_cm(12,:) = [200, 200, 200]; % Gray tone
roybigbl_cm(13,:) = [200, 200, 200]; % Gray tone
roybigbl_cm(14,:) = [200, 200, 200]; % Gray tone
roybigbl_cm(15,:) = [200, 200, 200]; % Gray tone
roybigbl_cm(16,:) = [255, 255, 255]; % Gray tone

roybigbl_cm = roybigbl_cm .* (1/255);
interpsteps = [0 0.0666 0.1333 0.2 0.2666 0.3333 0.4 0.4666 0.5333 0.6 0.6666 0.7333 0.8 0.86666 0.9333 1];
roybigbl_cm = interp1(interpsteps, roybigbl_cm, linspace(0, 1, 255));
roybigbl_cm = flipud(roybigbl_cm); % Yellow as high
custommap = roybigbl_cm(1:255, :);


% Visualization
figure('Visible', 'off');
imagesc(Dnet_upsampled);
colormap(custommap);     % Use inferno colormap
caxis([-0.1 1]);            % White (-0.1–0), Gray (0–0.6), Inferno (0.6–1)
set(gca, 'Color', [1 1 1]); % Background white
axis off;
set(gca, 'XColor', 'none', 'YColor', 'none');

%  Save high-resolution image
output_img = '/scratch/users/apines/Mouse_DMN_Upsampled_MATLAB.png';
print(output_img, '-dpng', '-r600');
close;

