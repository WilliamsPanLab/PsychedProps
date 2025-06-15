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
% load in DMN
network = load('/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/gro_Nets_fs4.mat');
dmn_L = network.nets.Lnets(:,1);
dmn_R = network.nets.Rnets(:,1);
% load TSNR mask
% Load medial wall mask (TSNR-based)
mwAndTSNR_L = gifti('/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii');
mw_L = zeros(2562, 1);
mw_L(mwAndTSNR_L.cdata(:,1) == 1) = 1;
% desaturate jet for vector visibility
jetmap = jet(256);
desaturated_jet = jetmap * 0.5 + 0.35;
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
ret     = zeros(size(V_L));     % [n_vertices x 3]
retD    = zeros(size(V_L));
weights = zeros(size(V_L,1), 1);  % track how many times each vertex gets hit

% orthogonalize to surface of mesh
for i = 1:length(incl_L)
    face_idx = incl_L(i);                     % face index
    verts = F_L(face_idx, :);                 % 3 vertex indices of this face
    % --- NoDrug ---
    OGvec = v3D_NoDrug_L(i,:);
    v1 = V_L(verts(1), :);
    v2 = V_L(verts(2), :);
    v3 = V_L(verts(3), :);
    normal = cross(v2 - v1, v3 - v1);
    normal = normal / norm(normal);
    proj = dot(OGvec, normal) * normal;
    modVec = OGvec - proj;
    for v = verts
        ret(v, :) = ret(v, :) + modVec;
        weights(v) = weights(v) + 1;
    end
    % --- Drug ---
    OGvecD = v3D_Drug_L(i,:);
    projD = dot(OGvecD, normal) * normal;
    modVecD = OGvecD - projD;

    for v = verts
        retD(v, :) = retD(v, :) + modVecD;
    end
end

% DIFFERENCE IN ANGLES MASK
% Compute angular difference between normalized vectors
unit_ret  = ret  ./ vecnorm(ret, 2, 2);
unit_retD = retD ./ vecnorm(retD, 2, 2);

% Clamp dot product to valid acos range
dotprod = sum(unit_ret .* unit_retD, 2);
dotprod = max(min(dotprod, 1), -1);

% Angular difference in degrees
angle_diff = acosd(dotprod);

% Threshold: keep only vertices where angle ≥ 35 degrees
angle_mask = angle_diff >= 35;

% Apply mask: zero out small-difference vectors
%ret(~angle_mask, :)  = 0;
retD(~angle_mask, :) = 0;
nonzero=vecnorm(ret,2,2)>0;
nonzero_D = vecnorm(retD, 2, 2) > 0;
ret(nonzero,:) = ret(nonzero,:) ./ vecnorm(ret(nonzero,:),2,2);
retD(nonzero_D, :) = retD(nonzero_D, :) ./ vecnorm(retD(nonzero_D, :), 2, 2);

% make figure
figure('Color','w','Position',[100 100 1600 800]);
subplot(1,2,1)
trisurf(F_L, V_L(:,1), V_L(:,2), V_L(:,3), ...
    dmn_L, 'EdgeColor', 'none', 'FaceAlpha', 1);  % Use DMN as vertex color
colormap(desaturated_jet);
hold on;

% Overlay NoDrug vector field
quiver3D(V_L(:,1), V_L(:,2), V_L(:,3), ...
         ret(:,1), ret(:,2), ret(:,3), ...
         repmat([0.2 0.2 0.9], size(V_L,1), 1), 2);
hold on;
% Overlay Drug vector field
quiver3D(V_L(:,1), V_L(:,2), V_L(:,3), ...
         retD(:,1), retD(:,2), retD(:,3), ...
         repmat([0.9 0.2 0.2], size(V_L,1), 1), 2);
hold on;
view([90 0]);  % Medial LH
axis equal off tight;
lighting none; camlight headlight; material dull;

% ----------- Lateral View (Right Panel) -----------
subplot(1,2,2)
trisurf(F_L, V_L(:,1), V_L(:,2), V_L(:,3), ...
    dmn_L, 'EdgeColor', 'none', 'FaceAlpha', 1);  % Reuse DMN colors
colormap(desaturated_jet);
hold on;

% Overlay NoDrug vector field again
quiver3D(V_L(:,1), V_L(:,2), V_L(:,3), ...
         ret(:,1), ret(:,2), ret(:,3), ...
         repmat([0.2 0.2 0.9], size(V_L,1), 1), 2);
hold on;
% Overlay NoDrug vector field again
quiver3D(V_L(:,1), V_L(:,2), V_L(:,3), ...
         retD(:,1), retD(:,2), retD(:,3), ...
         repmat([0.9 0.2 0.2], size(V_L,1), 1), 2);

view([-90 0]);  % Lateral LH
axis equal off tight;
lighting none; camlight headlight; material dull;

print('/scratch/users/apines/DrugVNoDrug_L.png','-dpng','-r400');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now do right hemisphere
ret     = zeros(size(V_R));     % [n_vertices x 3]
retD    = zeros(size(V_R));
weights = zeros(size(V_R,1), 1);  % track how many times each vertex gets hit

% orthogonalize to surface of mesh
for i = 1:length(incl_R)
    face_idx = incl_R(i);                     % face index
    verts = F_R(face_idx, :);                 % 3 vertex indices of this face
    % --- NoDrug ---
    OGvec = v3D_NoDrug_R(i,:);
    v1 = V_R(verts(1), :);
    v2 = V_R(verts(2), :);
    v3 = V_R(verts(3), :);
    normal = cross(v2 - v1, v3 - v1);
    normal = normal / norm(normal);
    proj = dot(OGvec, normal) * normal;
    modVec = OGvec - proj;
    for v = verts
        ret(v, :) = ret(v, :) + modVec;
        weights(v) = weights(v) + 1;
    end
    % --- Drug ---
    OGvecD = v3D_Drug_R(i,:);
    projD = dot(OGvecD, normal) * normal;
    modVecD = OGvecD - projD;

    for v = verts
        retD(v, :) = retD(v, :) + modVecD;
    end
end

% DIFFERENCE IN ANGLES MASK
% Compute angular difference between normalized vectors
unit_ret  = ret  ./ vecnorm(ret, 2, 2);
unit_retD = retD ./ vecnorm(retD, 2, 2);
% Clamp dot product to valid acos range
dotprod = sum(unit_ret .* unit_retD, 2);
dotprod = max(min(dotprod, 1), -1);
% Angular difference in degrees
angle_diff = acosd(dotprod);
% Threshold: keep only vertices where angle ≥ 35 degrees
angle_mask = angle_diff >= 35;
% Apply mask: zero out small-difference vectors
%ret(~angle_mask, :)  = 0;
retD(~angle_mask, :) = 0;
nonzero=vecnorm(ret,2,2)>0;
nonzero_D = vecnorm(retD, 2, 2) > 0;
ret(nonzero,:) = ret(nonzero,:) ./ vecnorm(ret(nonzero,:),2,2);
retD(nonzero_D, :) = retD(nonzero_D, :) ./ vecnorm(retD(nonzero_D, :), 2, 2);

% make figure
figure('Color','w','Position',[100 100 1600 800]);
subplot(1,2,1)
trisurf(F_R, V_R(:,1), V_R(:,2), V_R(:,3), ...
    dmn_R, 'EdgeColor', 'none', 'FaceAlpha', 1);  % Use DMN as vertex color
colormap(desaturated_jet);
hold on;

% Overlay NoDrug vector field
quiver3D(V_R(:,1), V_R(:,2), V_R(:,3), ...
         ret(:,1), ret(:,2), ret(:,3), ...
         repmat([0.2 0.2 0.9], size(V_R,1), 1), 2);
hold on;
% Overlay Drug vector field
quiver3D(V_R(:,1), V_R(:,2), V_R(:,3), ...
         retD(:,1), retD(:,2), retD(:,3), ...
         repmat([0.9 0.2 0.2], size(V_R,1), 1), 2);
hold on;
view([90 0]);  % Medial LH
axis equal off tight;
lighting none; camlight headlight; material dull;

% ----------- Lateral View (Right Panel) -----------
subplot(1,2,2)
trisurf(F_R, V_R(:,1), V_R(:,2), V_R(:,3), ...
    dmn_R, 'EdgeColor', 'none', 'FaceAlpha', 1);  % Reuse DMN colors
colormap(desaturated_jet);
hold on;

% Overlay NoDrug vector field again
quiver3D(V_R(:,1), V_R(:,2), V_R(:,3), ...
         ret(:,1), ret(:,2), ret(:,3), ...
         repmat([0.2 0.2 0.9], size(V_R,1), 1), 2);
hold on;
% Overlay NoDrug vector field again
quiver3D(V_R(:,1), V_R(:,2), V_R(:,3), ...
         retD(:,1), retD(:,2), retD(:,3), ...
         repmat([0.9 0.2 0.2], size(V_R,1), 1), 2);

view([-90 0]);  % Lateral LH
axis equal off tight;
lighting none; camlight headlight; material dull;
print('/scratch/users/apines/DrugVNoDrug_R.png','-dpng','-r400');


%%%%%%%%%%%%%%%%%%%% print out medial wall
% Load surface
[V_L, F_L] = read_surf('/oak/stanford/groups/leanew1/users/apines/surf/lh.sphere');
F_L = F_L + 1;

% Load medial wall mask (vertex-wise)
mwAndTSNR_L = gifti('/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii');
mw_vtx_L = zeros(2562,1);
mw_vtx_L(mwAndTSNR_L.cdata == 1) = 1;

% Convert vertex-wise to face-wise mask
F_MW_L = sum(mw_vtx_L(F_L), 2) ./ 3;
F_MW_L = ceil(F_MW_L);  % any face with ≥1 MW vertex = MW
face_vals_L = ones(size(F_L,1),1);        % default: white
face_vals_L(F_MW_L==1) = 0;               % MW: black

% Plot
figure('Color','w','Position',[100 100 1600 800]);

% Medial view
subplot(1,2,1)
trisurf(F_L, V_L(:,1), V_L(:,2), V_L(:,3), ...
        face_vals_L, 'EdgeColor','none', 'FaceAlpha', 1);
colormap([0 0 0; 1 1 1]); caxis([0 1]);  % 0=black, 1=white
view([90 0]); axis equal off tight;
title('LH Medial - MW in Black');

% Lateral view
subplot(1,2,2)
trisurf(F_L, V_L(:,1), V_L(:,2), V_L(:,3), ...
        face_vals_L, 'EdgeColor','none', 'FaceAlpha', 1);
colormap([0 0 0; 1 1 1]); caxis([0 1]);
view([-90 0]); axis equal off tight;
title('LH Lateral - MW in Black');

print('/scratch/users/apines/medial_wall_faces_LH.png','-dpng','-r400');

% Load surface
[V_R, F_R] = read_surf('/oak/stanford/groups/leanew1/users/apines/surf/rh.sphere');
F_R = F_R + 1;

% Load medial wall mask (vertex-wise)
mwAndTSNR_R = gifti('/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii');
mw_vtx_R = zeros(2562,1);
mw_vtx_R(mwAndTSNR_R.cdata == 1) = 1;

% Convert vertex-wise to face-wise mask
F_MW_R = sum(mw_vtx_R(F_R), 2) ./ 3;
F_MW_R = ceil(F_MW_R);
face_vals_R = ones(size(F_R,1),1);
face_vals_R(F_MW_R==1) = 0;

% Plot
figure('Color','w','Position',[100 100 1600 800]);

% Medial view
subplot(1,2,1)
trisurf(F_R, V_R(:,1), V_R(:,2), V_R(:,3), ...
        face_vals_R, 'EdgeColor','none', 'FaceAlpha', 1);
colormap([0 0 0; 1 1 1]); caxis([0 1]);
view([-90 0]); axis equal off tight;
title('RH Medial - MW in Black');

% Lateral view
subplot(1,2,2)
trisurf(F_R, V_R(:,1), V_R(:,2), V_R(:,3), ...
        face_vals_R, 'EdgeColor','none', 'FaceAlpha', 1);
colormap([0 0 0; 1 1 1]); caxis([0 1]);
view([90 0]); axis equal off tight;
title('RH Lateral - MW in Black');

print('/scratch/users/apines/medial_wall_faces_RH.png','-dpng','-r400');

