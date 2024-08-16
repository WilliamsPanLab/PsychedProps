% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% list the Diaz mice and location
parentfolder='/scratch/users/apines/data/mouse/Diaz/';

% list of Diaz mice
mList={'m7520','m7522','m7531','m7576','m7590'};
% initialize DMN mask on mouse brain parallel to matrix construction in extract OpFl features scripts
networks=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
Dnet=networks.components_matrix';
Mask=networks.mask;
% need to re-boolean from the interpolation in the downsample
Mask=Mask==1;
% need to transpose mask to match python-matlab difference: confirmed via visualization in Extract_DMNMag and Extract_relativeangles
Mask=Mask';
% get size of mask and indices to include
MaskSize=sum(sum(Mask));
MaskInds=find(Mask);
% use fullMask to initially mask 67*70 grid, 
fullMask = Mask;
% Now that we have the brain mask, get the DMN-specific mask
nets=Dnet;
% create pixel-wise network mask
DMN_bool=Dnet;
DMN_bool(DMN_bool>.6)=1;
DMN_bool(DMN_bool<.6)=0;
% adding mask into this step 7/1/24
DMN_bool(Mask==0)=0;
DMN_bool=logical(DMN_bool);

% initialize grid for mouse pre Diaz
preDiazgrid=zeros(67,70);
% initialize grid for mouse post Diaz
postDiazgrid=zeros(67,70);
% for each mouse
for i=1:length(mList)
	% retrieve current mouse name
	mouse=mList{i};
	% load in mouse tp1
	mousePreDiazfp=[parentfolder mouse '_1_Mag_TS_dmn.csv'];
	mousePreDiaz=dlmread(mousePreDiazfp);
	% grid to populate and reinitialize iteratively
	iterGrid=zeros(67,70,size(mousePreDiaz,2));
	% for each timepoint
	for tp=1:size(mousePreDiaz,2)
		tempGrid = zeros(67, 70); % create a temporary grid for this time point
	        tempGrid(DMN_bool) = mousePreDiaz(:, tp); % apply mask to the temporary grid
	        iterGrid(:, :, tp) = tempGrid; % assign the temporary grid to the appropriate time slice
	end
	% Calculate the percentage of observations that are < 90 across all timepoints
	iterGrid_pct=mean(iterGrid,3);
	% plug into master grid
	preDiazgrid=preDiazgrid+iterGrid_pct;
end

% average
preDiazgrid=preDiazgrid./(length(mList));

% load in mouse post
% for each mouse
for i=1:length(mList)
        % retrieve current mouse name
        mouse=mList{i};
	seshs=[3 4 5];
	% initialize iterGrid_pct for across sessions
	iterGrid_pct=zeros(67,70);
	% for each session
        for s=seshs;
		% -2 because list starts with 3
		sesh=seshs(s-2);
		% load in mouse tp1
        	mousePostDiazfp=[parentfolder mouse '_' num2str(sesh) '_Mag_TS_dmn.csv'];
        	mousePostDiaz=dlmread(mousePostDiazfp);
        	% grid to populate and reinitialize iteratively
        	iterGrid=zeros(67,70,size(mousePostDiaz,2));
        	% for each timepoint
        	for tp=1:size(mousePostDiaz,2)
        	        tempGrid = zeros(67, 70); % create a temporary grid for this time point
        	        tempGrid(DMN_bool) = mousePostDiaz(:, tp); % apply mask to the temporary grid
        	        iterGrid(:, :, tp) = tempGrid; % assign the temporary grid to the appropriate time slice
        	end
        	iterGrid_pct=iterGrid_pct+mean(iterGrid,3);
        end
	% average by number of sessions
	iterGrid_pct=iterGrid_pct./length(seshs);	
	% plug into master grid
        postDiazgrid=postDiazgrid+iterGrid_pct;
end
% average for each mouse
postDiazgrid=postDiazgrid./(length(mList));

% get subtraction of pre - post 
preMinPostDiaz=preDiazgrid-postDiazgrid;

% add in colormap info
% blue-orange color scheme
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts')
BO_cm=inferno(9);
BO_cm(1,:)=[49 197 244];
BO_cm(2,:)=[71 141 203];
BO_cm(3,:)=[61 90 168];
BO_cm(4,:)=[64 104 178];
BO_cm(5,:)=[126 126 126];
BO_cm(6,:)=[240 74 35];
BO_cm(7,:)=[243 108 33];
BO_cm(8,:)=[252 177 11];
BO_cm(9,:)=[247 236 31];
% scale to 1
BO_cm=BO_cm.*(1/255);
% interpolate color gradient
interpsteps=[0 .125 .25 .375 .5 .625 .75 .875 1];
BO_cm=interp1(interpsteps,BO_cm,linspace(0,1,255));
custommap=BO_cm;

% plot pre - post
figure
imagesc(preMinPostDiaz);
colormap(custommap);
maxcol=1.3;
mincol=-1.3;
caxis([mincol maxcol]); 
colorbar
hold on;
contour(Mask, [0.5 0.5], 'k', 'LineWidth', 1.5); % Add a black contour line
print('~/Diaz_preMinPost_Mag.png','-dpng')
