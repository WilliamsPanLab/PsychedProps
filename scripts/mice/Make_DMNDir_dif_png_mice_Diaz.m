% load data
resultant = readtable('~/mouse_resultant_vector_differences_Diaz.csv');
% load in mask
networks=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
Dnet=networks.components_matrix';
Mask=networks.mask;
% need to re-boolean from the interpolation in the downsample
Mask=Mask==1;
% need to transpose mask to match python-matlab difference: to be confirmed via visualization
Mask=Mask';
% get size of mask and indices to include
MaskSize=sum(sum(Mask));
MaskInds=find(Mask);
% load in dmn gradient for reference indices to mimic extract relative angles
nets=Dnet;
% create pixel-wise network mask
DMN_bool=Dnet;
DMN_bool(DMN_bool>.6)=1;
DMN_bool(DMN_bool<.6)=0;
% adding mask into this step 7/1/24
DMN_bool(Mask==0)=0;
DMN_bool=logical(DMN_bool);
DMN_inds = find(DMN_bool);
% convert tos
[rows, cols] = ind2sub(size(DMN_bool), DMN_inds);
% combine coords into one variable
DMN_coords = [rows, cols];

% Initialize vector field for the difference
vecMap_delta = zeros(size(DMN_bool,1), size(DMN_bool,2), 2);

% Populate difference map at DMN locations
for i = 1:length(DMN_inds)
    r = rows(i); c = cols(i);

    dx = resultant.Delta_X(i);  % Cartesian x-component
    dy = resultant.Delta_Y(i);  % Cartesian y-component

    vecMap_delta(r, c, :) = [dx, dy];
end

% Extract x/y delta components
U_delta = vecMap_delta(:,:,1);
V_delta = vecMap_delta(:,:,2);

% Create coordinate grid
[xx, yy] = meshgrid(1:size(DMN_bool,2), 1:size(DMN_bool,1));

% gray base 
gray_base = ones(size(Mask)) * 126/255;  % default gray background
gray_img = repmat(gray_base, [1 1 3]);   % RGB image
% dmn region to white
% Overwrite DMN region with white
for ch = 1:3
    tmp = gray_img(:,:,ch);
    tmp(DMN_bool) = .9;  % set to white
    gray_img(:,:,ch) = tmp;
end

% Plot
fig = figure('Color','w','Position',[100 100 700 700]);

%show custom RGB image
image(gray_img); axis equal tight ij; hold on;

% Add contour at Mask boundary
contour(Mask, [0.5 0.5], 'k', 'LineWidth', 1.5);

quiver(xx(Mask), yy(Mask), ...
       U_delta(Mask), V_delta(Mask), ...
       0.75, 'Color', [0.9373, 0.5843, 0], 'LineWidth', 1.5);

% Save
print(fig, '/scratch/users/apines/mouse_resultant_vector_difference_Diaz.png', '-dpng', '-r400');
close(fig);
