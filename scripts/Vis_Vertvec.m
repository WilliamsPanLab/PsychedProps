function Vis_Vertvec(VertVecL,VertVecR,Fn) 

addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))

%%% Load in surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfL = [SubjectsFolder '/lh.sphere'];
surfR = [SubjectsFolder '/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
%faces_l = faces_l + 1;
%faces_r = faces_r + 1;
% faces_L
%F_L=faces_l;
% vertices V
%V_L=vx_l;
% faces_R
%F_R=faces_r;
% vertices V
%V_R=vx_r;

% use native freesurfer command for mw mask indices
%surfML = [SubjectsFolder '/lh.Medial_wall.label'];
%mwIndVec_l = read_medial_wall_label(surfML);
%surfMR = [SubjectsFolder '/rh.Medial_wall.label'];
%mwIndVec_r = read_medial_wall_label(surfMR);
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
mw_L=ones(1,2562);
mw_L(mwAndTSNR_L==1)=0;
mw_R=ones(1,2562);
mw_R(mwAndTSNR_R==1)=0;
mwIndVec_l=find(mw_L);
mwIndVec_r=find(mw_R);


%%%%%%%%%%%%%%%%%%%%%%%%
data=zeros(1,2562);
data(mwIndVec_l)=VertVecL;
%data(mwIndVec_l)=0;

%%%%%%% fixed colorscale varities

%%% circular
mincol=min(VertVecL);
maxcol=max(VertVecL);
% for resultant vector distance mapping (Calc_cirdist)
mincol=.005;
maxcol=.025;
mincol=-7;
maxcol=7;
% for nmf networks
%mincol=0;
%maxcol=1;
% for t-stats
%  mincol=-9;
%  maxcol=9;
%%% for red/blue 0-centered
%mincol=-9;
%maxcol=9;
%custommap=colormap(b2r(mincol,maxcol));
% abscense of color to gray to accom. lighting "none"
%  grayColor = [0.7, 0.7, 0.7];  % Define gray color
% Add gray color to the colormap
%   custommap = [custommap; grayColor];
%custommap(126,:)=[.5 .5 .5];
custommap=colormap(jet);

% blue-orange color scheme
%BO_cm=inferno(9);
%BO_cm(1,:)=[49 197 244];
%BO_cm(2,:)=[71 141 203];
%BO_cm(3,:)=[61 90 168];
%BO_cm(4,:)=[64 104 178];
%BO_cm(5,:)=[126 126 126];
%BO_cm(6,:)=[240 74 35];
%BO_cm(7,:)=[243 108 33];
%BO_cm(8,:)=[252 177 11];
%BO_cm(9,:)=[247 236 31];
% scale to 1
%BO_cm=BO_cm.*(1/255);
% interpolate color gradient
%interpsteps=[0 .125 .25 .375 .5 .625 .75 .875 1];
%BO_cm=interp1(interpsteps,BO_cm,linspace(0,1,255));
%custommap=BO_cm;
%custommap=colormap(inferno);

%%% matches circular hist
% for 180 degree max
%roybigbl_cm=inferno(6);
%roybigbl_cm(1,:)=[0, 0, 255];
%roybigbl_cm(2,:)=[0, 255, 255];
%roybigbl_cm(3,:)=[116, 192, 68];
%roybigbl_cm(4,:)=[246, 235, 20];
%roybigbl_cm(5,:)=[255, 165, 0];
%roybigbl_cm(6,:)=[255, 0, 0];
% scale to 1
%roybigbl_cm=roybigbl_cm.*(1/255);
% interpolate color gradient
%interpsteps=[0 .2 .4 .6 .8 1];
%roybigbl_cm=interp1(interpsteps,roybigbl_cm,linspace(0,1,255));
% make circular with flipud
%custommap=vertcat(flipud(roybigbl_cm),roybigbl_cm);

% for 90 degree max
%roybigbl_cm=inferno(3);
% blue
%roybigbl_cm(1,:)=[0, 0, 255];
% cyan
%roybigbl_cm(2,:)=[0, 255, 255];
% green
%roybigbl_cm(3,:)=[116, 192, 68];
% yellow
%roybigbl_cm(4,:)=[246, 235, 20];
% scale to 1
%roybigbl_cm=roybigbl_cm.*(1/255);
% interpolate color gradient
%interpsteps=[0 .33333 .66666 1];
%interpsteps=[0 .5 1];
%roybigbl_cm=interp1(interpsteps,roybigbl_cm,linspace(0,1,255));
% make circular with flipud
%custommap=vertcat(flipud(roybigbl_cm),roybigbl_cm);
%custommap=roybigbl_cm;
%custommap=colormap(parula);
% mw to black
%custommap(1,:)=[0 0 0];

%%%%%%% THIS TRICK IS TO ADD IN THE GRADIENT OF THE SMOOTH DMN ONTO THE INFLATED SURFACE
% load in left cortical surface
% add lukas lang functions
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/lukaslang-ofd-614a2ffc50d6'))
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/lh.inflated']);
F_L=faces;
V_L=vertices;
% LOAD IN DMN STUFF
network=load(['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/Smooth_Nets_fs4.mat']);
n_LH=network.nets.Lnets(:,1);
n_RH=network.nets.Rnets(:,1);
% calculate network gradients on sphere
ng_L = grad(F_L, V_L, n_LH);
% convert both back to vertices for angular comparisons
vertwise_grad_L=zeros(2562,3);
% for each vertex, grab adjacent face values and merge em
for v=1:2562
	[InvolvedFaces_l,~]=find(F_L==v);
	vertwise_grad_L(v,:)=mean(ng_L(InvolvedFaces_l,:),1);
end
% ORTHOGONALIZE TO SURFACE OF INFLATED
ret=vertwise_grad_L;
% for each vector, subtract weighted surface-orthogonal component from original vector
for v=1:length(vertices)
        % retrieve original vector
        OGvec=ret(v,:);
        % find the three faces involved in this vertex 
        [InvolvedFaces,~]=find(faces==v);
        % get normal vectors of each involved face
        normalVectors = cross(vertices(faces(InvolvedFaces, 2), :) - vertices(faces(InvolvedFaces, 1), :), vertices(faces(InvolvedFaces, 3), :) - vertices(faces(InvolvedFaces, 1), :));
        % find vector as close as possible to orthogonal from these three faces
        meanNormalVector = mean(normalVectors, 1);
        % normalize normal vector
        meanNormalVector=VecNormalize(meanNormalVector);
        % get dot product of orthogonal vector and original vector
        OGvecOrthogonal = dot(OGvec, meanNormalVector) * meanNormalVector;
        % subtract orthogonal component of original vector from original vector
        modVec = OGvec - OGvecOrthogonal;;
        % add modified vector to initialized matrix
        ret(v,:)=modVec;
end
% normalize vectors for equal length
ret=VecNormalize(ret);
ret(mwIndVec_l,:)=0;
%%% END ORTHOGONALIZATION OF LEFT DMN GRAD TO INFLATED SURFACE

figure
asub = subaxis(2,2,1, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);

aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
view([90 0]);
colormap(custommap)
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
shading flat;
camlight;
	alpha(1)

set(gca,'CLim',[mincol,maxcol]);
% COMMENTED OUT VFS FOR NOW
%bplot=quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),[.3 .5 .7])
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');

asub = subaxis(2,2,4, 'sh', 0.00, 'sv', 0.00, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
%bplot=quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),[.3 .5 .7])
view([90 0]);
rotate(aplot, [0 0 1], 180)
%rotate(bplot, [0 0 1], 180)
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
material metal %shiny %metal;
shading flat;
camlight;
alpha(1)
 pos = get(asub, 'Position');
 posnew = pos; posnew(2) = posnew(2) + 0.13; posnew(1) = posnew(1) -.11; set(asub, 'Position', posnew);
set(gcf,'Color','w')

set(gca,'CLim',[mincol,maxcol]);
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');


%%% right hemisphere
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/rh.inflated']);
V_R=vertices;
F_R=faces;

data=zeros(1,2562);
data(mwIndVec_r)=VertVecR;
% calculate network gradients on sphere
ng_R = grad(F_R, V_R, n_RH);
% convert both back to vertices for angular comparisons
vertwise_grad_R=zeros(2562,3);
% for each vertex, grab adjacent face values and merge em
for v=1:2562
        [InvolvedFaces_r,~]=find(F_R==v);
        vertwise_grad_R(v,:)=mean(ng_R(InvolvedFaces_r,:),1);
end
% ORTHOGONALIZE TO SURFACE OF INFLATED
ret=vertwise_grad_R;
% for each vector, subtract weighted surface-orthogonal component from original vector
for v=1:length(vertices)
        % retrieve original vector
        OGvec=ret(v,:);
        % find the three faces involved in this vertex 
        [InvolvedFaces,~]=find(faces==v);
        % get normal vectors of each involved face
        normalVectors = cross(vertices(faces(InvolvedFaces, 2), :) - vertices(faces(InvolvedFaces, 1), :), vertices(faces(InvolvedFaces, 3), :) - vertices(faces(InvolvedFaces, 1), :));
        % find vector as close as possible to orthogonal from these three faces
        meanNormalVector = mean(normalVectors, 1);
        % normalize normal vector
        meanNormalVector=VecNormalize(meanNormalVector);
        % get dot product of orthogonal vector and original vector
        OGvecOrthogonal = dot(OGvec, meanNormalVector) * meanNormalVector;
        % subtract orthogonal component of original vector from original vector
        modVec = OGvec - OGvecOrthogonal;;
        % add modified vector to initialized matrix
        ret(v,:)=modVec;
end
% normalize vectors for equal length
ret=VecNormalize(ret);
ret(mwIndVec_r,:)=0;

%mincol=min(data);
%maxcol=max(data);

asub = subaxis(2,2,2, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0,'Holdaxis',1);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
%bplot = quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),[.3 .5 .7])
view([90 0]);
rotate(aplot, [0 0 1], 180)
%rotate(bplot, [0 0 1], 180)
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none
material metal %shiny %metal;%shading flat;
shading flat;
camlight;
 pos = get(asub, 'Position');
 posnew = pos; posnew(1) = posnew(1) - 0.11; set(asub, 'Position', posnew);
alpha(1)


set(gca,'CLim',[mincol,maxcol]);
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');

asub = subaxis(2,2,3, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
%bplot=quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),[.3 .5 .7])
view([90 0]);
colormap(custommap)
caxis([mincol; maxcol]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
material metal %shiny %metal;
shading flat;
camlight;
alpha(1)
 pos = get(asub, 'Position');
 posnew = pos; posnew(2) = posnew(2) + 0.13; set(asub, 'Position', posnew);
set(gcf,'Color','w')


set(gca,'CLim',[mincol,maxcol]);
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');
%colorbar
c=colorbar
%c.Location='southoutside'

colormap(custommap)

print(Fn,'-dpng','-r300')
