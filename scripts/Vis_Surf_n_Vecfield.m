function Vis_Surf_n_Vecfield(surfl,surfr,vecl,vecr,Fn) 
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))
% designed for fsaverage4 surface
VertVecL=surfl;
%%% Load in surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfL = [SubjectsFolder '/lh.pial'];
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

%%%%%%%%%%%%%%%%%%%%%%%%
% set data in the plotting script's terms
data=VertVecL;
ret=vecl;
% vector scaling factor
scalingfactor=1;
% colors
mincol=0;
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
RGBValues=(vecl+abs(minVal))./(maxVal+abs(minVal));

%%% for red/blue 0-centered
%mincol=-9;
%maxcol=9;
%custommap=colormap(b2r(mincol,maxcol));
% abscense of color to gray to accom. lighting "none"
%custommap(126,:)=[.5 .5 .5];
%custommap=colormap(jet)

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

% flatten vectors to surface
% medial left hemisphere
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/lh.inflated']);
% make VertVecL parallel to pial surface as "ret"
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
        ret(v,:)=VecNormalize(modVec);
end
% replace nans with 0
ret(isnan(ret))=0;
% set Rgb color map
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
% 0 to gray
%custommap=colormap('jet');
% medial left hemisphere
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/lh.inflated']);
% begin figure
figure
asub = subaxis(2,2,1, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
hold on;
quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),RGBValues,scalingfactor)
view([90 0]);
colormap(custommap)
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
shading flat;
camlight;

set(gca,'CLim',[mincol,maxcol]);
aplot.FaceVertexCData=data;
aplot.FaceAlpha=.3;

% other view of left hemisphere (lateral)
asub = subaxis(2,2,4, 'sh', 0.00, 'sv', 0.00, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3))
view([90 0]);
hold on;
colormap(custommap);
caxis([mincol; maxcol]);
% removing upscaling vector field so vectors are locked to vertices
bplot=quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),RGBValues,scalingfactor)
daspect([1 1 1]);
axis tight;
rotate(aplot, [0 0 1], 180)
rotate(bplot, [0 0 1], 180)

axis vis3d off;
lighting none;
material metal %shiny %metal;
shading flat;
camlight;

% insert RGB colors onto surface
aplot.FaceVertexCData=data;
aplot.FaceAlpha=.3;
% printout
print(Fn,'-dpng','-r2000')




%%% right hemisphere
% data=VertVecR;

% [vertices, faces] = freesurfer_read_surf([SubjectsFolder '/rh.inflated']);

% asub = subaxis(2,2,2, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0,'Holdaxis',1);
% aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
% view([90 0]);
% rotate(aplot, [0 0 1], 180)
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
