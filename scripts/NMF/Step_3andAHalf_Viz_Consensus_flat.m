% load in consensus and visualize it
ProjectFolder = '/oak/stanford/groups/leanew1/users/apines/data'
resultantFolder = [ProjectFolder '/RobustInitialization'];
Consensus=load([resultantFolder '/init.mat']);

% load in mask to repopulate consensus onto masked grid
fn='/scratch/users/apines/p50_mice/proc/20200228/thy1gc6s_0p3mgkg_m2000_preLSD0p3mgkg_1/masked_dff_Gro_Masked_BP_Smoothed.h5'
formask=h5read(fn, '/mask');
mask = cellfun(@(x) strcmpi(x, 'TRUE'), formask);

% blue-orange color scheme
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


% visualize each component
for k=1:13
	figure
	% make outgrid to saveout with imwrite
	outgrid=zeros(200,208);
	outgrid(mask)=Consensus.initV(:,k);
	outgrid(outgrid<.2)=0;
	outputfilename=['~/Component_' num2str(k) '.png'];
	imagesc(outgrid);
	colormap(custommap);
	maxcol=1;
	mincol=-1;
	caxis([mincol maxcol]);
	hold on;
	contour(mask, [0.5 0.5], 'k', 'LineWidth', 1.5);
	colorbar
	print(['~/Component_' num2str(k) '.png'],'-dpng')
end
