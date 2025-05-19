% Load DMN components and mask
networks = load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
Dnet = networks.components_matrix';
Mask = networks.mask;

% Re-boolean from interpolation and transpose to match Python-MATLAB orientation
Mask = Mask == 1;
Mask = Mask';

% Compute gradient of the DMN values
[gradX, gradY] = gradient(Dnet);

reate the DMN mask (threshold >0.6)
DMN_bool = Dnet > 0.6;
DMN_bool(Mask == 0) = 0;
DMN_bool = logical(DMN_bool);

% Zero gradient vectors outside of the DMN mask
gradX(DMN_bool == 0) = 0;
gradY(DMN_bool == 0) = 0;

% Set non-brain regions to gray by assigning a specific low value
Dnet(~DMN_bool) = 0;
Dnet(Mask == 0) = -0.1;  % Adjust this value if needed for the colormap mapping

% Prepare grid for quiver plot
[xGrid, yGrid] = meshgrid(1:size(Dnet, 2), 1:size(Dnet, 1));

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

%  Plot DMN underlay with vector field overlay
figure('Visible', 'off'); % For headless cluster environments
imagesc(Dnet);
colormap(custommap);
caxis([-0.1 1]); % Adjust if needed based on Dnet range
set(gca, 'Color', [0.8 0.8 0.8]); % Gray background for non-masked areas
hold on;

% Overlay the gradient as a vector field (black arrows)
quiver(xGrid, yGrid, gradX, gradY, 'k');
axis equal;

% Overlay DMN contour
contour(DMN_bool, [0.5 0.5], 'k', 'LineWidth', 1.5); % Black contour for DMN boundary

% Clean up figure: no axes, ticks, or labels
axis off;
set(gca, 'XColor', 'none', 'YColor', 'none');
set(gca, 'xtick', [], 'ytick', []);

% Save the final figure
output_file = '/scratch/users/apines/DMN_gradient_vector_custommap.png';
print(output_file, '-dpng', '-r600');
close;
