% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% list the LSD mice and location (note dex and diaz are in sep. folders)
parentfolder='/scratch/users/apines/data/mouse/';

% list of LSD mice
mList={'m2000','m7507','m7520','m7522','m7589','m7594'};

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

% initialize grid for mouse pre LSD
preLSDgrid=zeros(67,70);
% initialize grid for mouse post LSD
postLSDgrid=zeros(67,70);
% for each mouse
for i=1:length(mList)
	% retrieve current mouse name
	mouse=mList{i};
	% load in mouse tp1 (LSD)
	mousePreLSDfp=[parentfolder mouse '_1_Mag_TS_dmn.csv'];
	mousePreLSD=dlmread(mousePreLSDfp);
	% grid to populate and reinitialize iteratively
	iterGrid=zeros(67,70,size(mousePreLSD,2));
	% for each timepoint
	for tp=1:size(mousePreLSD,2)
		tempGrid = zeros(67, 70); % create a temporary grid for this time point
	        tempGrid(DMN_bool) = mousePreLSD(:, tp); % apply mask to the temporary grid
	        iterGrid(:, :, tp) = tempGrid; % assign the temporary grid to the appropriate time slice
	end
	% Calculate the percentage of observations that are < 90 across all timepoints
	iterGrid_pct=mean(iterGrid,3);
	% plug into master grid
	preLSDgrid=preLSDgrid+iterGrid_pct;
end

% average
preLSDgrid=preLSDgrid./(length(mList));

% load in mouse post (LSD)
% for each mouse
for i=1:length(mList)
        % retrieve current mouse name
        mouse=mList{i};
	% get all LSD scans for each mouse
	if mouse=='m7507'
		seshs=[2 3 4 5];
	else
		seshs=[2 3 4 5 6];
	end
	% initialize iterGrid_pct for across sessions
	iterGrid_pct=zeros(67,70);
	% for each session
        for s=seshs;
		% -1 because list starts with 2
		sesh=seshs(s-1);
		% load in mouse tp1 (LSD)
        	mousePostLSDfp=[parentfolder mouse '_' num2str(sesh) '_Mag_TS_dmn.csv'];
        	mousePostLSD=dlmread(mousePostLSDfp);
        	% grid to populate and reinitialize iteratively
        	iterGrid=zeros(67,70,size(mousePostLSD,2));
        	% for each timepoint
        	for tp=1:size(mousePostLSD,2)
        	        tempGrid = zeros(67, 70); % create a temporary grid for this time point
        	        tempGrid(DMN_bool) = mousePostLSD(:, tp); % apply mask to the temporary grid
        	        iterGrid(:, :, tp) = tempGrid; % assign the temporary grid to the appropriate time slice
        	end
		% Calculate the percentage of observations that are < 90 across all timepoints
		iterGrid_sesh=mean(iterGrid,3);
        	% Calculate the percentage of observations that are < 90 across all timepoints
        	iterGrid_pct = iterGrid_pct + iterGrid_sesh;
		
        end
	% average by number of sessions
	iterGrid_pct=iterGrid_pct./length(seshs);	
	% plug into master grid
        postLSDgrid=postLSDgrid+iterGrid_pct;
end
% average for each mouse
postLSDgrid=postLSDgrid./(length(mList));

% get subtraction of pre - post 
preMinPostLSD=preLSDgrid-postLSDgrid;

% new mouselist: Dex Pre
mDexList={'m1','m824','m961','m962','m964'};
% new mouselist: Diaz Pre
mDiazList={'m7520','m7522','m7531','m7576','m7590'};
% initialize grid for mouse pre LSD
preDexgrid=zeros(67,70);
% initialize grid for mouse pre LSD
preDiazgrid=zeros(67,70);
% for each Dex mouse
for i=1:length(mDexList)
        % retrieve current mouse name
        mouse=mDexList{i};
        % load in mouse tp1 (Dex)
        mousePreDexfp=[parentfolder 'Dex/' mouse '_1_Mag_TS_dmn.csv'];
        mousePreDex=dlmread(mousePreDexfp);
        % grid to populate and reinitialize iteratively
        iterGrid=zeros(67,70,size(mousePreDex,2));
        % for each timepoint
        for tp=1:size(mousePreDex,2)
                tempGrid = zeros(67, 70); % create a temporary grid for this time point
                tempGrid(DMN_bool) = mousePreDex(:, tp); % apply mask to the temporary grid
                iterGrid(:, :, tp) = tempGrid; % assign the temporary grid to the appropriate time slice
        end
	iterGrid_pct=mean(iterGrid,3);
        % plug into master grid
        preDexgrid=preDexgrid+iterGrid_pct;
end
% average
preDexgrid=preDexgrid./(length(mDexList));

% Diazepam
% for each Diaz mouse
for i=1:length(mDiazList)
        % retrieve current mouse name
        mouse=mDiazList{i};
        % load in mouse tp1 (Diazepam)
        mousePreDiazfp=[parentfolder 'Diaz/' mouse '_1_Mag_TS_dmn.csv'];
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
preDiazgrid=preDiazgrid./(length(mDiazList));

% average the 3 pre's
preMaster=preLSDgrid+preDexgrid+preDiazgrid;
preMaster=preMaster./3;

% get substraction of nodrug - LSD
NoDrugMinPostLSD=preMaster-postLSDgrid;

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
imagesc(preMinPostLSD);
colormap(custommap);
maxcol=(max(max(preMinPostLSD)));
mincol=-maxcol;
caxis([mincol maxcol]); 
colorbar
hold on;
contour(Mask, [0.5 0.5], 'k', 'LineWidth', 1.5); % Add a black contour line
print('~/LSD_preMinPost_Mag.png','-dpng')

% plot nodrug - LSD
figure
imagesc(NoDrugMinPostLSD);
colormap(custommap);
% maintain same color axis
caxis([mincol maxcol]);
colorbar
hold on;
contour(Mask, [0.5 0.5], 'k', 'LineWidth', 1.5); % Add a black contour line
print('~/LSD_NoDrugMinPost_Mag.png','-dpng')
