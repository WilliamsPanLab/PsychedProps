function Vis_Vec_only_fs5(VertVecL,VertVecR,Fn) 


% comment out trisurfs (aplot) if you need matlab to actually render everything. Appears to be pushing it's limits with this 3rd part vector field plotting.

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
mw_L=zeros(1,2562);
mw_L(mwAndTSNR_L>0)=1;
mw_R=zeros(1,2562);
mw_R(mwAndTSNR_R>0)=1;
mwIndVec_l=find(mw_L);
mwIndVec_r=find(mw_R);



%%%%%%% THIS TRICK IS TO ADD IN THE GRADIENT OF THE SMOOTH DMN ONTO THE INFLATED SURFACE
% add lukas lang functions
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/lukaslang-ofd-614a2ffc50d6'))
% load in left cortical surface
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/lh.inflated']);
F_L=faces;
V_L=vertices;
% ORTHOGONALIZE TO SURFACE OF INFLATED
ret=VertVecL;
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
%ret=VecNormalize(ret);
ret(mwIndVec_l,:)=0;
%%% END ORTHOGONALIZATION OF LEFT DMN GRAD TO INFLATED SURFACE

figure
asub = subaxis(2,2,1, 'sh', 0, 'sv', 0, 'padding', .04, 'margin', 0);
% scaled to be a little smaller with .99
aplot = trisurf(faces, vertices(:,1).*.99, vertices(:,2).*.99, vertices(:,3).*.99, ...
                ones(1,2562), 'FaceColor', 'w', 'EdgeColor', 'w', 'LineStyle', 'none');
view([90 0]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
shading flat;
camlight;
	alpha(1)
% map surface
hold on;
bplot=quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),[.5 .5 .5],2)
% next panel
asub = subaxis(2,2,4, 'sh', 0.00, 'sv', 0.00, 'padding', .04, 'margin', 0);
% scaled to be a little smaller with .99
aplot = trisurf(faces, vertices(:,1).*.99, vertices(:,2).*.99, vertices(:,3).*.99, ...
                ones(1,2562), 'FaceColor', 'w', 'EdgeColor', 'w', 'LineStyle', 'none');
hold on;
view([90 0]);
rotate(bplot, [0 0 1], 180)
rotate(aplot,[0 0 1],180);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
material metal %shiny %metal;
shading flat;
camlight;
alpha(1)
 pos = get(asub, 'Position');
 posnew = pos; posnew(2) = posnew(2) - 0.01; posnew(1) = posnew(1) +.03; set(asub, 'Position', posnew);
bplot=quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),[.5 .5 .5],2)
hold on;

%set(gca,'CLim',[mincol,maxcol]);
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');


%%% right hemisphere
[vertices, faces] = freesurfer_read_surf([SubjectsFolder '/rh.inflated']);
V_R=vertices;
F_R=faces;

% ORTHOGONALIZE TO SURFACE OF INFLATED
ret=VertVecR;
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
% mw mask
%ret(mwIndVec_r,:)=0;

asub = subaxis(2,2,2, 'sh', 0.0, 'sv', 0.0, 'padding', 0.04, 'margin', 0,'Holdaxis',1);
% scaled to be a little smaller with .99
aplot = trisurf(faces, vertices(:,1).*.99, vertices(:,2).*.99, vertices(:,3).*.99, ...
                ones(1,2562), 'FaceColor', 'w', 'EdgeColor', 'w', 'LineStyle', 'none');
view([90 0]);
rotate(aplot, [0 0 1], 180)
rotate(bplot, [0 0 1], 180)
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none
material metal %shiny %metal;%shading flat;
shading flat;
camlight;
 pos = get(asub, 'Position');
 posnew = pos; posnew(1) = posnew(1) + 0.03; set(asub, 'Position', posnew);
alpha(1)
bplot=quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),[.5 .5 .5],2)
hold on;
%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');

% next panel
asub = subaxis(2,2,3, 'sh', 0.0, 'sv', 0.0, 'padding', 0.04, 'margin', 0);
% scaled to be a little smaller with .99
aplot = trisurf(faces, vertices(:,1).*.99, vertices(:,2).*.99, vertices(:,3).*.99, ...
                ones(1,2562), 'FaceColor', 'w', 'EdgeColor', 'w', 'LineStyle', 'none');

view([90 0]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
material metal %shiny %metal;
shading flat;
camlight;
alpha(1)
 pos = get(asub, 'Position');
 posnew = pos; posnew(2) = posnew(2) - 0.01; set(asub, 'Position', posnew);
bplot=quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),[.5 .5 .5],2)
hold on;
%set(gcf,'Color','w')


%set(aplot,'FaceColor','flat','FaceVertexCData',data','CDataMapping','scaled');
%colorbar
%c=colorbar
%c.Location='southoutside'

print(Fn,'-dpng','-r1200')
