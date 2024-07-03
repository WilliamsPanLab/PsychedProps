% load in t-stats
DrugTs=readmatrix('/scratch/users/apines/taskVerts/valid_Drug_Ts_pixels.csv');

% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))

% load in supplemental info
% load in mask
networks=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
Dnet=networks.components_matrix';
Mask=networks.mask;
% need to re-boolean from the interpolation in the downsample
Mask=Mask==1;
% need to transpose mask to match python-matlab difference: to be confirmed via visualization
Mask=Mask';
% load in DMN
networks=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
Dnet=networks.components_matrix';
net=Dnet;
[nGx, nGy] =imgradientxy(net);

% Define the colormap limits
mincol = -9.2;
maxcol = 9.2;

% Create a custom colormap
custommap = b2r(mincol, maxcol);
whiteColor = [1, 1, 1];    % Define white color
grayColor = [0.7, 0.7, 0.7];  % Define gray color

% Add white and gray colors to the colormap
custommap = [whiteColor; custommap; grayColor];

% Define the outlier value
outlierValue = 999;

% Normalize the values in DrugTs to the range [0, 1], excluding outliers and zeros
normalizedDrugTs = (DrugTs - mincol) / (maxcol - mincol);
normalizedDrugTs(normalizedDrugTs < 0) = 0;  % Clip values below 0
normalizedDrugTs(normalizedDrugTs > 1) = 1;  % Clip values above 1

% Convert normalizedDrugTs to indices in the colormap
colormapIndices = round(normalizedDrugTs * (size(custommap, 1) - 3)) + 2; % +2 to account for white and gray

% Assign the unique colormap index to zero values (white) and outlier values (gray)
whiteIndex = 1;  % Index for the white color
grayIndex = size(custommap, 1);  % Index for the gray color

colormapIndices(DrugTs == 0) = whiteIndex;
colormapIndices(DrugTs == outlierValue) = grayIndex;

% actual figure code
fig=figure;
imagesc(colormapIndices)
colormap(custommap);
hold on;
% Create a grid for the quiver plot
[x, y] = meshgrid(1:size(net, 2), 1:size(net, 1));
quiver(x, y, nGx, nGy);
print('~/MouseLSD_pixelwiseStats.png','-dpng','-r600');
% Save the image using imwrite
%imwrite(colormapIndices, custommap, '~/MouseLSD_pixelwiseStats.png');
