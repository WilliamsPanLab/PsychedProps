% SOFTWARE REQUIRED FOR THIS DEMO

% 1. Matlab (tested on R2022b)
% 2. Freesurfer "read_surf" command: see https://surfer.nmr.mgh.harvard.edu/fswiki/
% 3. Spherical optical flow: see https://github.com/lukaslang/ofd
% 4. Cart2sphvec (Mathworks): see https://www.mathworks.com/products/phased-array.html

% to record runtime
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 1: simulate data
% add path to load in required software 2: freesurfer read_surf for matlab 
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% load in surface data: fsaverage4 left hemisphere spherical form
surfL = ['/oak/stanford/groups/leanew1/users/apines/surf/lh.sphere'];
[vx_l, faces_l] = read_surf(surfL);
% add one to faces because indexing starts at 0
faces_l=faces_l+1;
% initialize simulated data structure: a two timepoint time series
TS=zeros(size(vx_l,1),2);
% standard deviation for gaussian distribution
SD = 20;

% pull out y-coordinates of points
y_coords = vx_l(:, 2);

% for timepoint 1, generate gaussian distribution with appropriate center
TS(:,1)=exp(-((y_coords - centers(1)).^2) / (2 * SD^2));

% same thing for timepoint 2: gaussian shifted forward by 20
TS(:,2)=exp(-((y_coords - centers(2)).^2) / (2 * SD^2));

% print out the data as a trisurf 
figure;
subplot(1, 2, 1);
% ts 1 for timepoint 1
trisurf(faces_l, vx_l(:, 1), vx_l(:, 2), vx_l(:, 3), TS(:, 1), 'EdgeColor', 'none');
title('Timepoint 1');
xlabel('X'); ylabel('Y'); zlabel('Z');
colorbar;
% standard low-Y-left view
view(90,0);
% aspect ratio
daspect([1 1 1]);

% timepoint 2
subplot(1, 2, 2);
trisurf(faces_l, vx_l(:, 1), vx_l(:, 2), vx_l(:, 3), TS(:, 2), 'EdgeColor', 'none');
title('Timepoint 2');
xlabel('X'); ylabel('Y'); zlabel('Z');
colorbar;
view(90,0);
% aspect ratio
daspect([1 1 1]);
% print out to home directory
print('~/AP_simulated.png','-dpng','-r600');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 2: run spherical optical flow
% add path to require software 3: spherical optical flow
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% and apply standard parameters to this optical flow run
N = 10; % Vector spherical harmonics degree
h = 1; % finite-difference height vector of triangles
alpha = 1; % scales Tikhonov regularization, > alpha = > spatiotemporal smoothness
s = 1; % R(u), regularizing functional, scales Tikhonov regularization more rapidly via penalization of first spatial and first temporal derivative

% run optical flow- "of"
u = of(N, faces_l, vx_l, TS(:,1), TS(:,2), h, alpha, s);

% save output
save('~/DemoOFoutput.mat','u')

% get incenters of triangles for plotting
TR_l = TriRep(faces_l,vx_l);
IC_l = TR_l.incenters;

% print vector fields on sphere
figure;
% Note arrowheads are not extremely visible in matlab's base quiver call
quiver3(IC_l(:, 1), IC_l(:, 2), IC_l(:, 3), u(:, 1), u(:, 2), u(:, 3))
view(90,0);
% aspect ratio
daspect([1 1 1]);
% print out to home directory
print('~/AP_simulatedVecs.png','-dpng','-r600');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 3: extract angular distance from a vector field pointing A->P
% note Part 3 uses cart2sphvec from https://www.mathworks.com/products/phased-array.html
% use y coordinates as AP values
APValues=y_coords;
% take gradient of these values on a sphere
gAP = grad(faces_l, vx_l, APValues);
% extract x y and z components
gAP_x=gAP(:,1);
gAP_y=gAP(:,2);
gAP_z=gAP(:,3);
% initialize angular distance vectors per face
AngDists=zeros(5120,1);
% we'll need azimuth/elevation/rho coordinates of the spherical points to simplify angular distancs to 2 dimensions
[az_L,el_L,r_L]=cart2sph(IC_l(:,1),IC_l(:,2),IC_l(:,3));
% and convert these coordinates from radians to degrees (rho is negligible in a sphere: all points are equidistant from center)
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
% just as we can describe points on sphere in 2 dimensions, we can describe vectors on the surface of a sphere in 2 dimensions
% so we'll apply the same dimensionality reduction to the vector fields themselves. First to the gradient of AP values
% first initialize the output
gAP_az=zeros(5120,1);
gAP_el=zeros(5120,1);
for f=1:5120;
	gVecs=cart2sphvec(double([gAP_x(f);gAP_y(f);gAP_z(f)]),azd_L(f),eld_L(f));
	gAP_az(f)=gVecs(1);
	gAP_el(f)=gVecs(2);
end
% now we'll apply to same technique to the flow vectors
flow_az=zeros(5120,1);
flow_el=zeros(5120,1);
for f=1:5120;
	% recall that flow vectors are stored in u
        fVecs=cart2sphvec(double([u(f,1);u(f,2);u(f,3)]),azd_L(f),eld_L(f));
        flow_az(f)=fVecs(1);
        flow_el(f)=fVecs(2);
end
% initialize angular distance vector
AngDist=zeros(5120,1);
% extract angular distance of estimated flow vectors from this gradient for each face
for f=1:5120;
	% extract flow vector and gradient-of-A->P-vector
	flowVec=[flow_az(f), flow_el(f)];
	gradientVec=[gAP_az(f), gAP_el(f)];
	% calculate angular distance
	AngDist(f) = acosd(min(1,max(-1, gradientVec(:).' * flowVec(:) / norm(gradientVec) / norm(flowVec) )));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 4: extract vector magnitudes
% calculate vector magnitudes from crest and periphery of wave
CrestIndices = find(IC_l(:, 2) > -20 & IC_l(:, 2) < 20);
PeriphIndices = find(IC_l(:, 2) < -80 | IC_l(:, 2) > 80);
% calculate all magnitudes
magnitudes = sqrt(sum(u.^2, 2));
% pull out core and periphery magnitudes
CrestMags=magnitudes(CrestIndices);
PeriphMags=magnitudes(PeriphIndices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 5: evaluate outputs
% if magnitude reflects volume of coherent signal flow, we'd expect flow vectors to be higher near the crest of the wave than the periphery
mean(CrestMags)
mean(PeriphMags)
mean(CrestMags)>mean(PeriphMags)
% if optical flow directions are random w/r/t to the A->P direction, we'd expect angular distances to be below 90 degrees on average
% (90 degrees out of a range of possible angular distances from 0-180 reflects entropic dispersion of flow vectors w/r/t each other)
mean(AngDist)
mean(AngDist)<90
%%%%%%%%%%%%%%%%%%%%%%%

% just to prove the point, re-run calculation with x_coords (L->R axis) as reference vector field rather than y-coords (a->p axis) 
x_coords=vx_l(:, 1);
% use y coordinates as AP values
APValues=x_coords;
% take gradient of these values on a sphere
gAP = grad(faces_l, vx_l, APValues);
% extract x y and z components
gAP_x=gAP(:,1);
gAP_y=gAP(:,2);
gAP_z=gAP(:,3);
% initialize angular distance vectors per face
AngDists=zeros(5120,1);
% we'll need azimuth/elevation/rho coordinates of the spherical points to simplify angular distancs to 2 dimensions
[az_L,el_L,r_L]=cart2sph(IC_l(:,1),IC_l(:,2),IC_l(:,3));
% and convert these coordinates from radians to degrees (rho is negligible in a sphere: all points are equidistant from center)
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
% just as we can describe points on sphere in 2 dimensions, we can describe vectors on the surface of a sphere in 2 dimensions
% so we'll apply the same dimensionality reduction to the vector fields themselves. First to the gradient of AP values
% first initialize the output
gAP_az=zeros(5120,1);
gAP_el=zeros(5120,1);
for f=1:5120;
        gVecs=cart2sphvec(double([gAP_x(f);gAP_y(f);gAP_z(f)]),azd_L(f),eld_L(f));
        gAP_az(f)=gVecs(1);
        gAP_el(f)=gVecs(2);
end
% initialize angular distance vector
AngDist=zeros(5120,1);
% extract angular distance of estimated flow vectors from this gradient for each face
for f=1:5120;
        % extract flow vector and gradient-of-A->P-vector
        flowVec=[flow_az(f), flow_el(f)];
        gradientVec=[gAP_az(f), gAP_el(f)];
        % calculate angular distance
        AngDist(f) = acosd(min(1,max(-1, gradientVec(:).' * flowVec(:) / norm(gradientVec) / norm(flowVec) )));
end
mean(AngDist)
% display runtime for whole script
toc
% end of demo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

