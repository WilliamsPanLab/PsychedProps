% Convert mean theta CSVs to 2D cartesian vectors and visualize them using Vis_FaceVec
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts')

%%% PREPARE AVERAGE ANGLES
% Load theta angles
theta_L = readmatrix('~/mean_theta_by_drug_LH.csv');
theta_R = readmatrix('~/mean_theta_by_drug_RH.csv');
rho_L = readmatrix('~/mean_rho_by_drug_LH.csv');  % [Face, Rho_Drug, Rho_NoDrug]
rho_R = readmatrix('~/mean_rho_by_drug_RH.csv');

% Columns: [faceID, DrugTheta, NoDrugTheta]
theta_drug_L = theta_L(:,2);
theta_nodrug_L = theta_L(:,3);
theta_drug_R = theta_R(:,2);
theta_nodrug_R = theta_R(:,3);
rho_drug_L = rho_L(:,2);
rho_nodrug_L = rho_L(:,3);
rho_drug_R = rho_R(:,2);
rho_nodrug_R = rho_R(:,3);

% Convert to unit vectors in 2D tangent plane
[xD_L, yD_L] = pol2cart(theta_drug_L, rho_drug_L);
[xND_L, yND_L] = pol2cart(theta_nodrug_L, rho_nodrug_L);
[xD_R, yD_R] = pol2cart(theta_drug_R, rho_drug_R);
[xND_R, yND_R] = pol2cart(theta_nodrug_R, rho_nodrug_R);
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
desaturated_jet = jetmap * 0.4 + 0.45;
%%%%% LOAD mask from derivation
master_L = readmatrix('~/MasterMask_L_1.csv');
master_R = readmatrix('~/MasterMask_R_1.csv');
incl_L = find(master_L);
incl_R = find(master_R);
%% Reconstruct vectors in original 3D surface orientation
vDelta_L = zeros(length(incl_L), 3);
vDelta_R = zeros(length(incl_R), 3);
for i = 1:length(incl_L)
    f = incl_L(i);
    vDrug = sph2cartvec([azel_Drug_L(i,:) 0]', rad2deg(az_L(f)), rad2deg(el_L(f)));
    vNoDrug = sph2cartvec([azel_NoDrug_L(i,:) 0]', rad2deg(az_L(f)), rad2deg(el_L(f)));
    vDelta_L(i,:) = vDrug - vNoDrug;
end

for i = 1:length(incl_R)
    f = incl_R(i);
    vDrug = sph2cartvec([azel_Drug_R(i,:) 0]', rad2deg(az_R(f)), rad2deg(el_R(f)));
    vNoDrug = sph2cartvec([azel_NoDrug_R(i,:) 0]', rad2deg(az_R(f)), rad2deg(el_R(f)));
    vDelta_R(i,:) = vDrug - vNoDrug;
end


% build vertex - face mapping
nVertices_L = size(V_L, 1);
vDelta_L_vertex = zeros(nVertices_L, 3);
vDelta_L_count = zeros(nVertices_L, 1);
for i = 1:length(incl_L)
    face_idx = incl_L(i);
    vec = vDelta_L(i,:);  % the vector computed for that face
    verts = F_L(face_idx, :);  % the 3 vertices of that face

    for v = verts
        vDelta_L_vertex(v, :) = vDelta_L_vertex(v, :) + vec;
        vDelta_L_count(v) = vDelta_L_count(v) + 1;
    end
end
% Avoid divide-by-zero
vDelta_L_count(vDelta_L_count == 0) = 1;
% Final per-vertex vectors (averaged)
vDelta_L_vertex = vDelta_L_vertex ./ vDelta_L_count;
% right
nVertices_R = size(V_R, 1);
vDelta_R_vertex = zeros(nVertices_R, 3);
vDelta_R_count = zeros(nVertices_R, 1);
for i = 1:length(incl_R)
    face_idx = incl_R(i);
    vec = vDelta_R(i,:);  % the vector computed for that face
    verts = F_R(face_idx, :);  % the 3 vertices of that face
    for v = verts
        vDelta_R_vertex(v, :) = vDelta_R_vertex(v, :) + vec;
        vDelta_R_count(v) = vDelta_R_count(v) + 1;
    end
end
% Avoid divide-by-zero
vDelta_R_count(vDelta_R_count == 0) = 1;
% Final per-vertex vectors (averaged)
vDelta_R_vertex = vDelta_R_vertex ./ vDelta_R_count;

%%%%%%%%%%%%%%%%%% need to construe vectors in terms of vertices they point towards for transformation to inflated
%%% --- BUILD NEIGHBOR LIST: LEFT HEMISPHERE ---
nVertices_L = size(V_L, 1);
neighbors_L = cell(nVertices_L, 1);
for i = 1:size(F_L, 1)
    tri = F_L(i,:);
    for j = 1:3
        vj = tri(j);
        others = tri(setdiff(1:3, j));
        neighbors_L{vj} = unique([neighbors_L{vj}, others]);
    end
end

%%% --- BUILD NEIGHBOR LIST: RIGHT HEMISPHERE ---
nVertices_R = size(V_R, 1);
neighbors_R = cell(nVertices_R, 1);
for i = 1:size(F_R, 1)
    tri = F_R(i,:);
    for j = 1:3
        vj = tri(j);
        others = tri(setdiff(1:3, j));
        neighbors_R{vj} = unique([neighbors_R{vj}, others]);
    end
end
%%%%% now we need neighor wieghtings: vectors in temrs of neighbors
% get unique vertex indices from included faces
verts_from_faces_L = unique(F_L(incl_L, :));
verts_from_faces_R = unique(F_R(incl_R, :));
% initialize cell array to hold soft neighbor weights
neighbor_weights_L = cell(length(verts_from_faces_L), 1);
for ii = 1:length(verts_from_faces_L)
    v_idx = verts_from_faces_L(ii);      % actual vertex ID
    vec = vDelta_L_vertex(v_idx,:);                % row-matched vector (1 x 3)
    mag = norm(vec);
    neighs = neighbors_L{v_idx};         % neighbors of the vertex
    % Vectors from center to neighbors
    neigh_dirs = V_L(neighs,:) - V_L(v_idx,:);  % [nNeigh x 3]
    neigh_dirs_normed = neigh_dirs ./ vecnorm(neigh_dirs, 2, 2);  % normalize each row
    % Dot products: alignment of vec with each neighbor direction
    dots = neigh_dirs_normed * (vec' / mag);  % cosine similarity
    % Scale weights to preserve original magnitude
    weights = dots / norm(dots) * mag;
    % Store as [neighbor_id, weight]
    neighbor_weights_L{ii} = [neighs(:), weights(:)];
end
% right hemi
% initialize cell array to hold soft neighbor weights for right hemisphere
neighbor_weights_R = cell(length(verts_from_faces_R), 1);
for ii = 1:length(verts_from_faces_R)
    v_idx = verts_from_faces_R(ii);      % actual vertex ID
    vec = vDelta_R_vertex(v_idx,:);                % row-matched vector (1 x 3)
    mag = norm(vec);
    neighs = neighbors_R{v_idx};         % neighbors of the vertex
    % Vectors from center to neighbors
    neigh_dirs = V_R(neighs,:) - V_R(v_idx,:);  % [nNeigh x 3]
    neigh_dirs_normed = neigh_dirs ./ vecnorm(neigh_dirs, 2, 2);  % normalize each row
    % Dot products: alignment of vec with each neighbor direction
    dots = neigh_dirs_normed * (vec' / mag);  % cosine similarity
    % Scale weights to preserve original magnitude
    weights = dots / norm(dots) * mag;
    % Store as [neighbor_id, weight]
    neighbor_weights_R{ii} = [neighs(:), weights(:)];
end
% OK, now we need to reconstruct vectors on the inflated surface using these neighbor weights
surfL_inflated = '/oak/stanford/groups/leanew1/users/apines/surf/lh.inflated';
[V_L_inflated, F_L_inf] = read_surf(surfL_inflated);
F_L_inf = F_L_inf + 1;
% left: reconstruct vectors on the inflated surface using neighbor weights
recon_L = zeros(size(V_L_inflated));  % [n_vertices x 3]
for ii = 1:length(verts_from_faces_L)
    v_idx = verts_from_faces_L(ii);               % vertex ID on surface
    weights_and_ids = neighbor_weights_L{ii};     % [neighbor_id, weight]
    neigh_ids = weights_and_ids(:,1);
    weights   = weights_and_ids(:,2);
    % Directions to neighbors on *inflated* surface
    dir_vectors = V_L_inflated(neigh_ids,:) - V_L_inflated(v_idx,:);
    dir_vectors = dir_vectors ./ vecnorm(dir_vectors, 2, 2);  % unit vectors
    % Weighted sum of unit directions
    new_vec = sum(dir_vectors .* weights, 1);
    recon_L(v_idx, :) = new_vec;
end
% right: reconstruct
surfR_inflated = '/oak/stanford/groups/leanew1/users/apines/surf/rh.inflated';
[V_R_inflated, F_R_inf] = read_surf(surfR_inflated);
F_R_inf = F_R_inf + 1;
recon_R = zeros(size(V_R_inflated));  % [n_vertices x 3]
for ii = 1:length(verts_from_faces_R)
    v_idx = verts_from_faces_R(ii);               % vertex ID on surface
    weights_and_ids = neighbor_weights_R{ii};     % [neighbor_id, weight]
    neigh_ids = weights_and_ids(:,1);
    weights   = weights_and_ids(:,2);
    % Directions to neighbors on *inflated* surface
    dir_vectors = V_R_inflated(neigh_ids,:) - V_R_inflated(v_idx,:);
    dir_vectors = dir_vectors ./ vecnorm(dir_vectors, 2, 2);  % unit vectors
    % Weighted sum of unit directions
    new_vec = sum(dir_vectors .* weights, 1);
    recon_R(v_idx, :) = new_vec;
end
% orthogonalize to surface of mesh
% left 
ret_inflated_L = zeros(size(V_L_inflated));  % [n_vertices x 3]
for i = 1:length(incl_L)
    face_idx = incl_L(i);
    verts = F_L(face_idx, :);
    % Get reconstructed vector at center vertex (could choose mean of face too)
    OGvec = recon_L(verts(1), :);  % you can also average the 3 if desired
    % Face geometry on inflated surface
    v1 = V_L_inflated(verts(1), :);
    v2 = V_L_inflated(verts(2), :);
    v3 = V_L_inflated(verts(3), :);
    normal = cross(v2 - v1, v3 - v1);
    normal = normal / norm(normal);
    % Project out normal component
    proj = dot(OGvec, normal) * normal;
    modVec = OGvec - proj;
    for v = verts
        ret_inflated_L(v, :) = ret_inflated_L(v, :) + modVec;
    end
end
% Initialize output array for projected vectors
ret_inflated_R = zeros(size(V_R_inflated));  % [n_vertices x 3]
for i = 1:length(incl_R)
    face_idx = incl_R(i);
    verts = F_R(face_idx, :);
    % Get reconstructed vector at center vertex (could choose mean of face too)
    OGvec = recon_R(verts(1), :);  % you can also average the 3 if desired
    % Face geometry on inflated surface
    v1 = V_R_inflated(verts(1), :);
    v2 = V_R_inflated(verts(2), :);
    v3 = V_R_inflated(verts(3), :);
    normal = cross(v2 - v1, v3 - v1);
    normal = normal / norm(normal);
    % Project out normal component
    proj = dot(OGvec, normal) * normal;
    modVec = OGvec - proj;
    for v = verts
        ret_inflated_R(v, :) = ret_inflated_R(v, :) + modVec;
    end
end

% make figure
figure('Color','w','Position',[100 100 1600 800]);
subplot(1,2,1)
trisurf(F_L_inf, V_L_inflated(:,1), V_L_inflated(:,2), V_L_inflated(:,3), ...
    dmn_L, 'EdgeColor', 'none', 'FaceAlpha', 1);  % Use DMN as vertex color
colormap(desaturated_jet);
hold on;

% Overlay NoDrug vector field
quiver3D(V_L_inflated(:,1), V_L_inflated(:,2), V_L_inflated(:,3), ...
         ret_inflated_L(:,1), ret_inflated_L(:,2), ret_inflated_L(:,3), ...
         repmat([0.7569, 0.1294, 0.2235], size(V_L_inflated,1), 1), 2);
hold on;
view([90 0]);  % Medial LH
axis equal off tight;
lighting none; camlight headlight; material dull;

% ----------- Lateral View (Right Panel) -----------
subplot(1,2,2)
trisurf(F_L_inf, V_L_inflated(:,1), V_L_inflated(:,2), V_L_inflated(:,3), ...
    dmn_L, 'EdgeColor', 'none', 'FaceAlpha', 1);
colormap(desaturated_jet);
hold on;

% Overlay NoDrug vector field again
quiver3D(V_L_inflated(:,1), V_L_inflated(:,2), V_L_inflated(:,3), ...
         ret_inflated_L(:,1), ret_inflated_L(:,2), ret_inflated_L(:,3), ...
         repmat([0.7569, 0.1294, 0.2235], size(V_L_inflated,1), 1), 2);
hold on;
view([-90 0]);  % Lateral LH
axis equal off tight;
lighting none; camlight headlight; material dull;
print('/scratch/users/apines/DrugGrtr_L.png','-dpng','-r400');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now do right hemisphere
% make figure
figure('Color','w','Position',[100 100 1600 800]);
subplot(1,2,1)
trisurf(F_R_inf, V_R_inflated(:,1), V_R_inflated(:,2), V_R_inflated(:,3), ...
    dmn_R, 'EdgeColor', 'none', 'FaceAlpha', 1);  % Use DMN as vertex color
colormap(desaturated_jet);
hold on;

% Overlay NoDrug vector field
quiver3D(V_R_inflated(:,1), V_R_inflated(:,2), V_R_inflated(:,3), ...
         ret_inflated_R(:,1), ret_inflated_R(:,2), ret_inflated_R(:,3), ...
         repmat([0.7569, 0.1294, 0.2235], size(V_R_inflated,1), 1), 2);
hold on;
view([90 0]);  % Medial LH
axis equal off tight;
lighting none; camlight headlight; material dull;

% ----------- Lateral View (Right Panel) -----------
subplot(1,2,2)
trisurf(F_R_inf, V_R_inflated(:,1), V_R_inflated(:,2), V_R_inflated(:,3), ...
    dmn_R, 'EdgeColor', 'none', 'FaceAlpha', 1);  % Reuse DMN colors
colormap(desaturated_jet);
hold on;

% Overlay NoDrug vector field again
quiver3D(V_R_inflated(:,1), V_R_inflated(:,2), V_R_inflated(:,3), ...
         ret_inflated_R(:,1), ret_inflated_R(:,2), ret_inflated_R(:,3), ...
         repmat([0.7569, 0.1294, 0.2235], size(V_R_inflated,1), 1), 2);
hold on;
view([-90 0]);  % Lateral LH
axis equal off tight;
lighting none; camlight headlight; material dull;
print('/scratch/users/apines/DrugGrtr_R.png','-dpng','-r400');


%%%%%%%%%%%%%%%%%%%% print out medial wall
% Load surface
[V_L, F_L] = read_surf('/oak/stanford/groups/leanew1/users/apines/surf/lh.inflated');
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

% Lateral view
subplot(1,2,2)
trisurf(F_L, V_L(:,1), V_L(:,2), V_L(:,3), ...
        face_vals_L, 'EdgeColor','none', 'FaceAlpha', 1);
colormap([0 0 0; 1 1 1]); caxis([0 1]);
view([-90 0]); axis equal off tight;

print('/scratch/users/apines/medial_wall_faces_LH.png','-dpng','-r400');

% Load surface
[V_R, F_R] = read_surf('/oak/stanford/groups/leanew1/users/apines/surf/rh.inflated');
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

% Lateral view
subplot(1,2,2)
trisurf(F_R, V_R(:,1), V_R(:,2), V_R(:,3), ...
        face_vals_R, 'EdgeColor','none', 'FaceAlpha', 1);
colormap([0 0 0; 1 1 1]); caxis([0 1]);
view([90 0]); axis equal off tight;

print('/scratch/users/apines/medial_wall_faces_RH.png','-dpng','-r400');

