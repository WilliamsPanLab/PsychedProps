% calculate reference streams

% addpaths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))

% need TT - transverse temporal
% CS - central sulcus
% Calcarine


% load in fs4 and labels
%%% Load in surface data LEFT
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfL = [SubjectsFolder '/lh.pial'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
% faces_L
F_L=faces_l;
% vertices V
V_L=vx_l;
% read table
[v,label,ct]=read_annotation('/share/software/user/open/freesurfer/7.4.1/subjects/fsaverage4/label/lh.aparc.a2009s.annot');
% temporal transverse is last (76)
A1lab=ct.table(76,5);
% central sulcus is 47
CSlab=ct.table(47,5);
V1lab=ct.table(46,5);
% get indices
A1Inds=find(label==A1lab);
CSInds=find(label==CSlab);
V1inds=find(label==V1lab);
% together... they are approximately primary cortices
PrimaryInds=[A1Inds' CSInds' V1inds'];
% get geodesic distances
% convert faces into edge list
edges = [F_L(:,1), F_L(:,2); F_L(:,2), F_L(:,3); F_L(:,3), F_L(:,1)];
edges = unique(sort(edges, 2), 'rows'); % Remove duplicates
% compute euclidean edge weights: direct cortical distances but on cortical surface
edge_weights = vecnorm(V_L(edges(:,1), :) - V_L(edges(:,2), :), 2, 2);
% create a weighted graph representing path length on the surface
G = graph(edges(:,1), edges(:,2), edge_weights);
% intialize distances from primary cortex
geo_distances_L = inf(2562, 1); % Initialize distances to infinity
% compute shortest distance to any PrimaryInd vertex
for i = 1:2562
    geo_distances_L(i) = min(distances(G, i, PrimaryInds));
end

%%% Load in surface data RIGHT
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfR = [SubjectsFolder '/rh.pial'];
% surface topography
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_r = faces_r + 1;
% faces_L
F_R=faces_r;
% vertices V
V_R=vx_r;
% read table
[v,label,ct]=read_annotation('/share/software/user/open/freesurfer/7.4.1/subjects/fsaverage4/label/rh.aparc.a2009s.annot');
% temporal transverse is last (76)
A1lab=ct.table(76,5);
% central sulcus is 47
CSlab=ct.table(47,5);
V1lab=ct.table(46,5);
% get indices
A1Inds=find(label==A1lab);
CSInds=find(label==CSlab);
V1inds=find(label==V1lab);
% together... they are approximately primary cortices
PrimaryInds=[A1Inds' CSInds' V1inds'];
% get geodesic distances
% convert faces into edge list
edges = [F_R(:,1), F_R(:,2); F_R(:,2), F_R(:,3); F_R(:,3), F_R(:,1)];
edges = unique(sort(edges, 2), 'rows'); % Remove duplicates
% compute euclidean edge weights: direct cortical distances but on cortical surface
edge_weights = vecnorm(V_R(edges(:,1), :) - V_R(edges(:,2), :), 2, 2);
% create a weighted graph representing path length on the surface
G = graph(edges(:,1), edges(:,2), edge_weights);
% intialize distances from primary cortex
geo_distances_R = inf(2562, 1); % Initialize distances to infinity
% compute shortest distance to any PrimaryInd vertex
for i = 1:2562
    geo_distances_R(i) = min(distances(G, i, PrimaryInds));
end


% load in TSNR
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
mw_L=zeros(1,2562);
mw_L(mwAndTSNR_L>.1)=1;
mw_R=zeros(1,2562);
mw_R(mwAndTSNR_R>.1)=1;
mwIndVec_l=find(mw_L);
mwIndVec_r=find(mw_R);

% load in DMN



% data is geodesic distance, surface is pial
faces=F_L;
vertices=V_L;
data=geo_distances_L;
data(mwIndVec_l)=0;
% LEFT
figure;
asub = subaxis(2,2,1, 'sh', 0, 'sv', 0, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
view([90 0]);
daspect([1 1 1]);
axis tight;
axis vis3d off;
lighting none;
shading flat;
camlight;
        alpha(1)
asub = subaxis(2,2,4, 'sh', 0.00, 'sv', 0.00, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
view([90 0]);
rotate(aplot, [0 0 1], 180)
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
%%% right hemisphere
% data is geodesic distance, surface is pail
faces=F_R;
vertices=V_R;
data=geo_distances_R;
data(mwIndVec_r)=0;
asub = subaxis(2,2,2, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0,'Holdaxis',1);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
%bplot = quiver3D(vertices(:,1),vertices(:,2),vertices(:,3),ret(:,1), ret(:,2), ret(:,3),[.3 .5 .7])
view([90 0]);
rotate(aplot, [0 0 1], 180)
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
asub = subaxis(2,2,3, 'sh', 0.0, 'sv', 0.0, 'padding', 0, 'margin', 0);
aplot = trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3),data)
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
 posnew = pos; posnew(2) = posnew(2) + 0.13; set(asub, 'Position', posnew);
c=colorbar
c.Location='southoutside'
Fn=['~/DistanceFromPrimary.png'];
print(Fn,'-dpng','-r800')

