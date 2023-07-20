function OpFlStreamlines(subj,sesh,filename)

% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))

% set child filepath
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
% load in optical flow output
data=load([childfp '/' subj '_' sesh '_OpFl_rs.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfL = [SubjectsFolder '/lh.sphere'];
surfR = [SubjectsFolder '/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;


% use native freesurfer command for mw mask indices
surfML = [SubjectsFolder '/lh.Medial_wall.label'];
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = [SubjectsFolder '/rh.Medial_wall.label'];
mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "is medial wall" vector for vertices
mw_L=zeros(1,2562);
mw_L(mwIndVec_l)=1;
mw_R=zeros(1,2562);
mw_R(mwIndVec_r)=1;
% convert to faces
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



% normalize verts to unit sphere
% left
numV=length(vx_l);
vx_l(numV+1:end, :) = VecNormalize(vx_l(numV+1:end, :));
% right
numV=length(vx_r);
vx_r(numV+1:end, :) = VecNormalize(vx_r(numV+1:end, :));

% calculate normative directionality streamlines for this subject and session
n = size(faces_l, 1);
T = TriRep(faces_l, vx_l);
P = T.incenters./100;

%%% E as left Hemi
E=data.us.vf_left;
% initialize vector field to be plotted
plotVF=zeros(n,3);
% loop over to extract from unfortunate cell structure
for i=1:length(E);
	plotVF=plotVF+E{i};
end
% start with mean of each x y and z component
plotVF=plotVF./length(E);

% consider circular mean: might require cart2sphvec

% Set seed points.
[X, Y] = meshgrid(-1:0.04:1, -1:0.04:1);
idx = find(X.^2 + Y.^2 <= 1);
S = [X(idx), Y(idx)];

% Set parameters.
nmax = max(sqrt(sum((plotVF).^2, 2)));
h = 0.1/nmax;
maxit = 50;
lw = .7;

% Streamlines for first component.
v = plotVF;

%F = createFigure('summer', -1, 1, -1, 1);
% read in DMN as background
networks=load(['/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/group_Nets_fs4.mat']);
nets_LH=networks.nets.Lnets(:,2);

% scale cortical mantle
vx_l=vx_l./101;
% move z coord of cortical mantle backwards
% vx_l(:,3)=vx_l(:,3)-1;
figure
% add cortical mantle %
aplot = trisurf(faces_l, vx_l(:,1), vx_l(:,2), vx_l(:,3),nets_LH)
colormap('gray')
freezeColors;
hold on
streamlines3(P, v, S, h, maxit, 'summer', lw);
view(3);
daspect([1 1 1]);
%h=get(gca,'Children');
%set(gca,'Children',[h(76009) h(1:76008)]);
%adjustFigure;
% one rotation for insula
%savefigure(F, fullfile(childfp, filename), '-png', '-r600');
print(['~/streams/' filename],'-dpng')
%%%%%%%%%%%%%%


