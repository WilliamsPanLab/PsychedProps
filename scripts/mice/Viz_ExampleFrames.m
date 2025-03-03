function Viz_ExampleFrames(subj,sesh)
run=sesh;
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts')
% load in time series
basefp='/oak/stanford/groups/leanew1/users/apines/p50_mice/proc2/proc/20200228/'

% load in specified scan (only set to run with scan 1 for now)
if run==1
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_preLSD0p3mgkg_1/masked_dff_Full_BP_Smoothed_NoDS.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        	foratlas=h5read([basefp 'thy1gc6s_0p3mgkg_' subj '_preLSD0p3mgkg_1/masked_dff.h5'],'/atlas');
	else
                disp('no run found')
        end
end
%% add if/else
% post 1
if run==2
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_0/masked_dff.h5']
        if exist(fn)
               data=h5read(fn, '/vid');
                foratlas=h5read(fn, '/atlas');
        else
                disp('no run found')
        end
end

% post 2
if run==3
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_5/masked_dff.h5']
        if exist(fn)
                data=h5read(fn, '/vid');
                foratlas=h5read(fn, '/atlas');
        else
                disp('no run found')
        end
end
% post 3
if run==4
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_10/masked_dff.h5']
        if exist(fn)
                data=h5read(fn, '/vid');
                foratlas=h5read(fn, '/atlas');
        else
                disp('no run found')
        end
end
% post 4
if run==5
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_15/masked_dff.h5']
        if exist(fn)
                data=h5read(fn, '/vid');
                foratlas=h5read(fn, '/atlas');
        else
                disp('no run found')
        end
end
% post 5
if run==6
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_20/masked_dff.h5']
        if exist(fn)
                data=h5read(fn, '/vid');
                foratlas=h5read(fn, '/atlas');
        else
                disp('no run found')
        end
end

% pull out mask
Mask=foratlas<1;

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
% scale to 1
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

% for some timepoints
for t=40:50
    filename=['/scratch/users/apines/' subj '_' num2str(sesh) '_Signal_t' num2str(t)]; 
    % save out png (imagesc of signal in inferno)
    figure;
    current_data = data(:,:,t);
    % set zero values to an out-of-range value for colormap
 %   current_data(mask) = NaN;
    % create image
    %imagesc(current_data);
    %colormap(custommap);
    %caxis([-.1 .1])
    %caxis([0.0; .015]);
    colorbar
  
	zscored_data = (current_data - mean(current_data(:), 'omitnan')) / std(current_data(:), 'omitnan');
	% mask it
	zscored_data(Mask) = NaN;
	imagesc(zscored_data, 'AlphaData', ~isnan(zscored_data))
	% Plot with color axis set to ±3 SD
	colormap(custommap);
	caxis([-3 3]); % Set color limits to ±3 standard deviations
	% set background to white
	set(gca, 'Color', [1 1 1]); 
	hold on
	contour(Mask, [0.5 0.5], 'k', 'LineWidth', 1.5);
        % remove axis jumbles
        axis off;
        set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axis lines
        set(gca, 'xtick', [], 'ytick', []); 	
	print(filename,'-dpng','-r600')
end

