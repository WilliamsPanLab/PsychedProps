function Vis_Surf_n_Vecfield_Sphere_DMNGrad(Fn) 

addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))
% load in DMNGrad for VertVec values
network=load(['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/gro_Nets_fs4.mat']);
n_LH=network.nets.Lnets(:,1);
n_RH=network.nets.Rnets(:,1);
% designed for fsaverage4 surface
VertVecL=n_LH;
VertVecR=n_RH;
vecl=VertVecL;
vecr=VertVecR;
%%% Load in surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfL = [SubjectsFolder '/lh.sphere'];
surfR = [SubjectsFolder '/rh.sphere'];
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
% incenters
TR = TriRep(faces_l, vx_l);
P = TR.incenters;
TR = TriRep(faces_r, vx_r);
Pr = TR.incenters;


% set directional for vector field
Coloration='Directional';
% set for DMN loadings
Coloration='BOLD';

%%% load in medial wall and SNR mask (to zero out medial wall in trisurf and quvier3)
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
mw_L=zeros(1,2562);
mw_L(mwAndTSNR_L>.5)=1;
mw_R=zeros(1,2562);
mw_R(mwAndTSNR_R>.5)=1;
mw_L=logical(mw_L);
mw_R=logical(mw_R);
% "good" indices
nonMW_L=ones(1,2562);
nonMW_R=ones(1,2562);
nonMW_L(mw_L)=0;
nonMW_R(mw_R)=0;
nonMW_L=logical(nonMW_L);
nonMW_R=logical(nonMW_R);

% set MW for faces
% add TSNR mask, includes medial wall
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
mw_L=zeros(1,2562);
mw_L(mwAndTSNR_L==1)=1;
mw_R=zeros(1,2562);
mw_R(mwAndTSNR_R==1)=1;
% convert to faces
F_MW_L=sum(mw_L(faces_l),2)./3;
F_MW_R=sum(mw_R(faces_r),2)./3;
% convert "partial" medial wall to medial wall
F_MW_L=ceil(F_MW_L);
F_MW_R=ceil(F_MW_R);
% face mask indices
fmwIndVec_l=find(F_MW_L);
fmwIndVec_r=find(F_MW_R);
% make medial wall vector
g_noMW_combined_L=setdiff([1:5120],fmwIndVec_l);
g_noMW_combined_R=setdiff([1:5120],fmwIndVec_r);
% face-wise DMN mask
nets_LH=n_LH;
nets_RH=n_RH;
DMN_bool_L=sum(nets_LH(faces_l),2)./3;
DMN_bool_R=sum(nets_RH(faces_r),2)./3;
DMN_bool_L(DMN_bool_L>.3)=1;
DMN_bool_R(DMN_bool_R>.3)=1;
DMN_bool_L(DMN_bool_L<.3)=0;
DMN_bool_R(DMN_bool_R<.3)=0;
DMN_bool_L=logical(DMN_bool_L);
DMN_bool_R=logical(DMN_bool_R);
% and face-wise DMN values for coloration
DMN_L=sum(nets_LH(faces_l),2)./3;
DMN_R=sum(nets_RH(faces_r),2)./3;
DMN_R(DMN_R<.3)=0;

% get network gradient for vector field to plot
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/lukaslang-ofd-614a2ffc50d6'))
% calculate network gradients on sphere
ng_L = grad(F_L, V_L, n_LH);
ng_R = grad(F_R, V_R, n_RH);
ret=ng_L;
data=n_LH;
data(data<.3)=0;
scalingfactor=1;
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
	mincol=-.6;
	maxcol=.6;
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
	roybigbl_cm(8,:)=[200, 200, 200 ];
	roybigbl_cm(9,:)=[200, 200, 200 ];
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
	custommap=roybigbl_cm(1:255,:);
end

% medial left hemisphere
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/lh.sphere']);
% begin figure
figure
asub = subaxis(2,2,1, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);
if strcmp(Coloration,'Directional')
	aplot = trisurf(faces, vertices(:,1)*.99, vertices(:,2)*.99, vertices(:,3)*.99)
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
	quiver3D(P(g_noMW_combined_L,1),P(g_noMW_combined_L,2),P(g_noMW_combined_L,3),ret(g_noMW_combined_L,1), ret(g_noMW_combined_L,2), ret(g_noMW_combined_L,3),[0 0 0],scalingfactor)
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
        aplot = trisurf(faces, vertices(:,1)*.99, vertices(:,2)*.99, vertices(:,3)*.99)
	set(aplot, 'EdgeColor', 'none');
	set(aplot, 'FaceColor', 'w');
	set(aplot, 'LineStyle', 'none');
elseif strcmp(Coloration,'BOLD')
        aplot = trisurf(faces, vertices(:,1)*.99, vertices(:,2)*.99, vertices(:,3)*.99)
        aplot.FaceVertexCData=data;
end
view([90 0]);
hold on;
colormap(custommap);
caxis([mincol; maxcol]);
% removing upscaling vector field so vectors are locked to vertices
if strcmp(Coloration,'Directional')
	bplot=quiver3D(P(DMN_bool_L,1),P(DMN_bool_L,2),P(DMN_bool_L,3),ret(DMN_bool_L,1), ret(DMN_bool_L,2), ret(DMN_bool_L,3),[0 0 0],scalingfactor)
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
data=n_RH;
% lower thresh to account for face/vert dif
data(data<.3)=0;
ret=ng_R;
% mask with medial wall
%data(mw_R)=0;
%ret(mw_R,:)=0;
% medial left hemisphere
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/rh.sphere']);
% begin figure
asub = subaxis(2,2,2, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);
if strcmp(Coloration,'Directional') % trying *.99 so arrows get less cut off
        aplot = trisurf(faces, vertices(:,1)*.99, vertices(:,2)*.99, vertices(:,3)*.99)
        set(aplot, 'EdgeColor', 'none');
        set(aplot, 'FaceColor', 'w');
        set(aplot, 'LineStyle', 'none');
elseif strcmp(Coloration,'BOLD')
        aplot = trisurf(faces, vertices(:,1)*.99, vertices(:,2)*.99, vertices(:,3)*.99)
        %aplot.FaceVertexCData=data;
	set(aplot,'FaceColor','flat','FaceVertexCData',DMN_R,'CDataMapping','scaled');
end
hold on;
% create a quiver if this is a directional plot, fill-in if not
if strcmp(Coloration,'Directional')
        bplot=quiver3D(Pr(DMN_bool_R,1),Pr(DMN_bool_R,2),Pr(DMN_bool_R,3),ret(DMN_bool_R,1), ret(DMN_bool_R,2), ret(DMN_bool_R,3),[0 0 0],scalingfactor)
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
        aplot = trisurf(faces, vertices(:,1)*.99, vertices(:,2)*.99, vertices(:,3)*.99)
        set(aplot, 'EdgeColor', 'none');
        set(aplot, 'FaceColor', 'w');
        set(aplot, 'LineStyle', 'none');
elseif strcmp(Coloration,'BOLD')
        aplot = trisurf(faces, vertices(:,1)*.99, vertices(:,2)*.99, vertices(:,3)*.99)
        aplot.FaceVertexCData=data;
end
hold on;
colormap(custommap);
caxis([mincol; maxcol]);
% removing upscaling vector field so vectors are locked to vertices
if strcmp(Coloration,'Directional')
        bplot=quiver3D(Pr(:,1),Pr(:,2),Pr(:,3),ret(:,1), ret(:,2), ret(:,3),[0 0 0],scalingfactor)
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
print(Fn,'-dpng','-r800')




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
