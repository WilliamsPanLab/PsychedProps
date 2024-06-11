function Viz_Vecfields(subj,sesh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% Load in flatmouse opflow calc
% note this is for LSD only! Adapt recording date for ketamine if needed
childfp='/scratch/users/apines/p50_mice/proc/20200228/'
datafp=[childfp subj '_vf_out_' num2str(sesh) '.mat']
data=load(datafp);
vf=data.outputs.velocityFields;
% and actual signal
signalGrid=data.outputs.filteredSignal;
networks=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
Dnet=networks.components_matrix';
Mask=networks.mask;
% need to re-boolean from the interpolation in the downsample
Mask=Mask==1;
% need to transpose mask to match python-matlab difference: to be confirmed via visualization
Mask=Mask';

Vx = real(vf(:,:,600)); % X component
Vy = imag(vf(:,:,600)); % Y component

Vx(Mask==0)=0;
Vy(Mask==0)=0;

% Define grid for quiver plot
[xGrid, yGrid] = meshgrid(1:size(Vx, 2), 1:size(Vx, 1));

% Plot quiver plot with black vectors
figure;
quiver(xGrid, yGrid, Vx, Vy, 'k');
title('Velocity Field Quiver Plot');
xlabel('X Axis');
ylabel('Y Axis');
axis equal;

% Save the figure
filename = ['/scratch/users/apines/' subj '_' num2str(sesh) '_quiver_plot'];
print(filename, '-dpng', '-r600');
