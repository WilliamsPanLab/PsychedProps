function Viz_Vecfields(subj,sesh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% Load in flatmouse opflow calc
% note this is for LSD only!
childfp='/oak/stanford/groups/leanew1/users/apines/p50_mice/proc2/proc/20200228/'
datafp=[childfp subj '_vf_out_' num2str(sesh) '.mat']
data=load(datafp);
vf=data.outputs.velocityFields;
% and actual signal
signalGrid=data.outputs.filteredSignal;
% load in mouse DMN
networks=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
Dnet=networks.components_matrix';
Mask=networks.mask;
% need to re-boolean from the interpolation in the downsample
Mask=Mask==1;
% need to transpose mask to match python-matlab difference: to be confirmed via visualization
Mask=Mask';
% set signal colormap
% pull in connectome workbench colormap
roybigbl_cm=inferno(16);
roybigbl_cm(1,:)=[255, 255, 0 ];
roybigbl_cm(2,:)=[255, 200, 0];
roybigbl_cm(3,:)=[255, 120, 0];
roybigbl_cm(4,:)=[255, 0, 0 ];
roybigbl_cm(5,:)=[200, 0, 0 ];
roybigbl_cm(6,:)=[150, 0, 0 ];
roybigbl_cm(7,:)=[100, 0, 0 ];
roybigbl_cm(8,:)=[60, 0, 0 ];
roybigbl_cm(9,:)=[0, 0, 80 ];
roybigbl_cm(10,:)=[0, 0, 170];
roybigbl_cm(11,:)=[75, 0, 125];
roybigbl_cm(12,:)=[125, 0, 160];
roybigbl_cm(13,:)=[75, 125, 0];
roybigbl_cm(14,:)=[0, 200, 0];
roybigbl_cm(15,:)=[0, 255, 0];
roybigbl_cm(16,:)=[0, 255, 255];
roybigbl_cm=roybigbl_cm.*(1/255);
% interpolate color gradient
interpsteps=[0 0.0666 0.1333 0.2 0.2666 0.3333 0.4 0.4666 0.5333 0.6 0.6666 0.7333 0.8 0.86666 0.9333 1];
roybigbl_cm=interp1(interpsteps,roybigbl_cm,linspace(0,1,255));
% yellow as high
roybigbl_cm=flipud(roybigbl_cm);
% reduce just a little bit on the close-to-white coloring
custommap=roybigbl_cm(15:240,:);
% adding white at the beginning for zero values
%custommap = [1, 1, 1; custommap]; % Adding white at the beginning
for timepoint=20:50
	Vx = real(vf(:,:,timepoint)); % X component
	Vy = imag(vf(:,:,timepoint)); % Y component
	Vx(Mask==0)=0;
	Vy(Mask==0)=0;

	% Define grid for quiver plot
	[xGrid, yGrid] = meshgrid(1:size(Vx, 2), 1:size(Vx, 1));


	% pull ca2+ data
	current_data=signalGrid(:,:,timepoint);
	zscored_data = (current_data - mean(current_data(:), 'omitnan')) / std(current_data(:), 'omitnan');
	zscored_data(Mask == 0) = NaN;

	% make figure
	figure
	imagesc(zscored_data, 'AlphaData', ~isnan(zscored_data));
	colormap(custommap)
	caxis([-3 3]);
	hold on
	% pull OF data
	quiver(xGrid, yGrid, Vx, Vy, 'k');
	axis equal;
	contour(Mask, [0.5 0.5], 'k', 'LineWidth', 1.5);
	set(gca, 'Color', [1 1 1]); % Set background to white
	set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axis lines
	set(gca, 'xtick', [], 'ytick', []);
	% Save the figure
	filename = ['/scratch/users/apines/' subj '_' num2str(sesh) '_quiver_plot' num2str(timepoint) '.png'];
	print(filename, '-dpng', '-r600');
end


