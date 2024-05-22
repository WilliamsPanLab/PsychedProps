function Vis_Surf_n_Vecfield(surfl,surfr,vecl,vecr,Fn,Coloration) 
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))
% designed for fsaverage4 surface
VertVecL=surfl;
VertVecR=surfr;
%%% Load in surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfL = [SubjectsFolder '/lh.inflated'];
surfR = [SubjectsFolder '/rh.inflated'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% faces_L
F_L=faces_l;
% vertices V
V_L=vx_l;
% faces_R
F_R=faces_r;
% vertices V
V_R=vx_r;

%%% load in medial wall and SNR mask (to zero out medial wall in trisurf and quvier3)
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
mw_L=zeros(1,2562);
mw_L(mwAndTSNR_L==1)=1;
mw_R=zeros(1,2562);
mw_R(mwAndTSNR_R==1)=1;
mw_L=logical(mw_L);
mw_R=logical(mw_R);
%%%%%%%%%%%%%%%%%%%%%%%%
% set data in the plotting script's terms
data=VertVecL;
ret=vecl;
% mask with medial wall... leading to rendering issue in vector fields? commented out for now
%data(mw_L)=0;
%ret(mw_L,:)=0;
% vector scaling factor
scalingfactor=1;
% colors
mincol=0;
if strcmp(Coloration,'Directional')
	maxcol=max(vecl(:));
	% scale RGB values to max
	minVal = min(vecl(:));
	maxVal = max(vecl(:));
	% consider scaling RGB colors by this value as well
	mincol=0;
	maxcol=max(vecl(:));
	% scale RGB values to max
	minVal = min(vecl(:));
	maxVal = max(vecl(:));
	% need everything to be above 0 (+absminval) but scaled 0-1 (./maxval+absminval)
	RGBValues_L=(vecl+abs(minVal))./(maxVal+abs(minVal));
	% repeat for R
	maxcol=max(vecl(:));
        % scale RGB values to max
        minVal = min(vecr(:));
        maxVal = max(vecr(:));
        % consider scaling RGB colors by this value as well
        mincol=0;
        maxcol=max(vecr(:));
        % scale RGB values to max
        minVal = min(vecr(:));
        maxVal = max(vecr(:));
        % need everything to be above 0 (+absminval) but scaled 0-1 (./maxval+absminval)
        RGBValues_R=(vecr+abs(minVal))./(maxVal+abs(minVal));
elseif strcmp(Coloration,'BOLD')
	%mincol=min(min([VertVecL VertVecR]))
	%maxcol=max(max([VertVecL VertVecR]))
	mincol=-10;
	maxcol=10;
end

% replace nans with 0
ret(isnan(ret))=0;
if strcmp(Coloration,'Directional')
	% if colormap=directional, Rgb color map
	custommap = [
	    1 0 0;    % Red
	    0 1 0;    % Green
	    0 0 1;    % Blue
	    1 1 0;    % Yellow
	    0 1 1;    % Cyan
	    1 0 1;    % Magenta
	    0 0 0;    % Black
	    1 1 1     % White
	];
elseif strcmp(Coloration,'BOLD')
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
end
% medial left hemisphere
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/lh.inflated']);
% begin figure
figure
asub = subaxis(2,2,1, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);
if strcmp(Coloration,'Directional')
	aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
	set(aplot, 'EdgeColor', 'none');
	set(aplot, 'FaceColor', 'w');
	set(aplot, 'LineStyle', 'none');
elseif strcmp(Coloration,'BOLD')
	aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
	aplot.FaceVertexCData=data;
end
hold on;
% create a quiver if this is a directional plot, fill-in if not
if strcmp(Coloration,'Directional')
	quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),RGBValues_L,scalingfactor)
elseif strcmp(Coloration,'BOLD')
	quiver3D(0,0,0,0, 0, 0,0,scalingfactor)
end
view([90 0]);
colormap(custommap)
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
shading flat;
set(gca,'CLim',[mincol,maxcol]);
% other view of left hemisphere (lateral)
asub = subaxis(2,2,4, 'sh', 0.00, 'sv', 0.00, 'padding', 0, 'margin', 0);
if strcmp(Coloration,'Directional')
        aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
	set(aplot, 'EdgeColor', 'none');
	set(aplot, 'FaceColor', 'w');
	set(aplot, 'LineStyle', 'none');
elseif strcmp(Coloration,'BOLD')
        aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
        aplot.FaceVertexCData=data;
end
view([90 0]);
hold on;
colormap(custommap);
caxis([mincol; maxcol]);
% removing upscaling vector field so vectors are locked to vertices
if strcmp(Coloration,'Directional')
	bplot=quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),RGBValues_L,scalingfactor)
elseif strcmp(Coloration,'BOLD')
        bplot=quiver3D(0,0,0,0, 0, 0,0,scalingfactor)
end
daspect([1 1 1]);
axis tight;
rotate(aplot, [0 0 1], 180)
rotate(bplot, [0 0 1], 180)

axis vis3d off;
lighting none;
shading flat;


% RIGHT hemisphere
data=VertVecR;
ret=vecr;
% mask with medial wall
%data(mw_R)=0;
%ret(mw_R,:)=0;
% medial left hemisphere
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/rh.inflated']);
% begin figure
asub = subaxis(2,2,2, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);
if strcmp(Coloration,'Directional')
        aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
        set(aplot, 'EdgeColor', 'none');
        set(aplot, 'FaceColor', 'w');
        set(aplot, 'LineStyle', 'none');
elseif strcmp(Coloration,'BOLD')
        aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
        aplot.FaceVertexCData=data;
end
hold on;
% create a quiver if this is a directional plot, fill-in if not
if strcmp(Coloration,'Directional')
        bplot=quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),RGBValues_R,scalingfactor)
elseif strcmp(Coloration,'BOLD')
        bplot=quiver3D(0,0,0,0, 0, 0,0,scalingfactor)
end
view([90 0]);
rotate(aplot, [0 0 1], 180)
rotate(bplot, [0 0 1], 180)
colormap(custommap)
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
shading flat;
set(gca,'CLim',[mincol,maxcol]);
% other view of left hemisphere (lateral)
asub = subaxis(2,2,3, 'sh', 0.00, 'sv', 0.00, 'padding', 0, 'margin', 0);
if strcmp(Coloration,'Directional')
        aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
        set(aplot, 'EdgeColor', 'none');
        set(aplot, 'FaceColor', 'w');
        set(aplot, 'LineStyle', 'none');
elseif strcmp(Coloration,'BOLD')
        aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
        aplot.FaceVertexCData=data;
end
hold on;
colormap(custommap);
caxis([mincol; maxcol]);
% removing upscaling vector field so vectors are locked to vertices
if strcmp(Coloration,'Directional')
        bplot=quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),RGBValues_R,scalingfactor)
elseif strcmp(Coloration,'BOLD')
        bplot=quiver3D(0,0,0,0, 0, 0,0,scalingfactor)
end
view([90 0]);
daspect([1 1 1]);
axis tight;

axis vis3d off;
lighting none;
shading flat;


% printout
print(Fn,'-dpng','-r1000')




%%% right hemisphere
% data=VertVecR;

% [vertices, faces] = freesurfer_read_surf([SubjectsFolder '/rh.inflated']);

% asub = subaxis(2,2,2, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0,'Holdaxis',1);
% aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
% view([90 0]);
%colormap(custommap)
% caxis([mincol; maxcol]);
% daspect([1 1 1]);
% axis tight;
% axis vis3d off;
% lighting none
% material metal %shiny %metal;%shading flat;
% shading flat;
% camlight;
%  pos = get(asub, 'Position');
%  posnew = pos; posnew(1) = posnew(1) - 0.11; set(asub, 'Position', posnew);
% alpha(1)


% set(gca,'CLim',[mincol,maxcol]);
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');

% asub = subaxis(2,2,3, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);
% aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
% view([90 0]);
%colormap(custommap)
% caxis([mincol; maxcol]);
% daspect([1 1 1]);
% axis tight;
% axis vis3d off;
% lighting none;
% material metal %shiny %metal;
% shading flat;
% camlight;
% alpha(1)
% pos = get(asub, 'Position');
%  posnew = pos; posnew(2) = posnew(2) + 0.13; set(asub, 'Position', posnew);
% set(gcf,'Color','w')


% set(gca,'CLim',[mincol,maxcol]);
%%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');
% colorbar
% c=colorbar
% c.Location='southoutside'

% colormap(custommap)

% print(Fn,'-dpng')
