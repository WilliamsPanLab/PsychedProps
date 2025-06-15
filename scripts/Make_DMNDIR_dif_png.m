% Convert mean theta CSVs to 2D cartesian vectors and visualize them using Vis_FaceVec
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts')

%%% PREPARE AVERAGE ANGLES
% Load theta angles
theta_L = readmatrix('~/mean_theta_by_drug_LH.csv');
theta_R = readmatrix('~/mean_theta_by_drug_RH.csv');
% Columns: [faceID, DrugTheta, NoDrugTheta]
theta_drug_L = theta_L(:,2);
theta_nodrug_L = theta_L(:,3);
theta_drug_R = theta_R(:,2);
theta_nodrug_R = theta_R(:,3);
% Convert to unit vectors in 2D tangent plane
[xD_L, yD_L] = pol2cart(theta_drug_L, 1);
[xND_L, yND_L] = pol2cart(theta_nodrug_L, 1);
[xD_R, yD_R] = pol2cart(theta_drug_R, 1);
[xND_R, yND_R] = pol2cart(theta_nodrug_R, 1);
azel_Drug_L = [xD_L yD_L];
azel_NoDrug_L = [xND_L yND_L];
azel_Drug_R = [xD_R yD_R];
azel_NoDrug_R = [xND_R yND_R];

% get spherical coordinates of orignal surface
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
surfL = ['/oak/stanford/groups/leanew1/users/apines/surf/lh.sphere'];
surfR = ['/oak/stanford/groups/leanew1/users/apines/surf/rh.sphere'];
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
% get incenters of triangles
TR_L = TriRep(F_L,V_L);
P_L = TR_L.incenters;
TR_R = TriRep(F_R,V_R);
P_R = TR_R.incenters;
[az_L, el_L, ~] = cart2sph(P_L(:,1), P_L(:,2), P_L(:,3));
[az_R, el_R, ~] = cart2sph(P_R(:,1), P_R(:,2), P_R(:,3));


%%%%% LOAD mask from derivation
master_L = readmatrix('~/MasterMask_L_1.csv');
master_R = readmatrix('~/MasterMask_R_1.csv');
incl_L = find(master_L);
incl_R = find(master_R);

%% Reconstruct vectors in original 3D surface orientation
v3D_Drug_L = zeros(length(incl_L), 3);
v3D_NoDrug_L = zeros(length(incl_L), 3);
v3D_Drug_R = zeros(length(incl_R), 3);
v3D_NoDrug_R = zeros(length(incl_R), 3);

for i = 1:length(incl_L)
    f = incl_L(i);
    v3D_Drug_L(i,:)   = sph2cartvec([azel_Drug_L(i,:) 0]', rad2deg(az_L(f)), rad2deg(el_L(f)));
    v3D_NoDrug_L(i,:) = sph2cartvec([azel_NoDrug_L(i,:) 0]', rad2deg(az_L(f)), rad2deg(el_L(f)));
end

for i = 1:length(incl_R)
    f = incl_R(i);
    v3D_Drug_R(i,:)   = sph2cartvec([azel_Drug_R(i,:) 0]', rad2deg(az_R(f)), rad2deg(el_R(f)));
    v3D_NoDrug_R(i,:) = sph2cartvec([azel_NoDrug_R(i,:) 0]', rad2deg(az_R(f)), rad2deg(el_R(f)));
end
% load in DMN
network = load('/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/gro_Nets_fs4.mat');
dmn_L = network.nets.Lnets(:,1);
dmn_R = network.nets.Rnets(:,1);




print(f, '~/DMN_Direction_Drug.png', '-dpng', '-r400');


