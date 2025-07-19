% Convert mean theta CSVs to 2D cartesian vectors and visualize them using Vis_FaceVec
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/OpFl_CDys/scripts')

%%% PREPARE AVERAGE mags
% Load theta angles
mag_L = readmatrix('~/mean_Mags_by_drug_LH.csv');
mag_R = readmatrix('~/mean_Mags_by_drug_RH.csv');

% Columns: [faceID, Drugmag, NoDrugmag]
mag_drug_L = mag_L(:,2);
mag_nodrug_L = mag_L(:,3);
mag_drug_R = mag_R(:,2);
mag_nodrug_R = mag_R(:,3);
difL=mag_nodrug_L-mag_drug_L;
difR=mag_nodrug_R-mag_drug_R;
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
% load in DMN
network = load('/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/gro_Nets_fs4.mat');
dmn_L = network.nets.Lnets(:,1);
dmn_R = network.nets.Rnets(:,1);
% load TSNR mask
% Load medial wall mask (TSNR-based)
mwAndTSNR_L = gifti('/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii');
mw_L = zeros(2562, 1);
mw_L(mwAndTSNR_L.cdata(:,1) == 1) = 1;
% blue-orange color scheme
BO_cm=inferno(9);
BO_cm(1,:)=[49 197 244];
BO_cm(2,:)=[71 141 203];
BO_cm(3,:)=[61 90 168];
BO_cm(4,:)=[64 104 178];
BO_cm(5,:)=[126 126 126];
BO_cm(6,:)=[240 74 35];
BO_cm(7,:)=[243 108 33];
BO_cm(8,:)=[252 177 11];
BO_cm(9,:)=[247 236 31];
% scale to 1
BO_cm=BO_cm.*(1/255);
% interpolate color gradient
interpsteps=[0 .125 .25 .375 .5 .625 .75 .875 1];
BO_cm=interp1(interpsteps,BO_cm,linspace(0,1,255));
custommap=BO_cm;
mincol=-.75;
maxcol=.75;
set(gca,'CLim',[mincol,maxcol]);
%%%%% LOAD mask from derivation
master_L = readmatrix('~/MasterMask_L_1.csv');
master_R = readmatrix('~/MasterMask_R_1.csv');
incl_L = find(master_L);
incl_R = find(master_R);
% OK, now we need to reconstruct vectors on the inflated surface using these neighbor weights
surfL_inflated = '/oak/stanford/groups/leanew1/users/apines/surf/lh.inflated';
[V_L_inflated, F_L_inf] = read_surf(surfL_inflated);
F_L_inf = F_L_inf + 1;
dataL=zeros(1,5120);
dataL(incl_L)=difL;
% right: reconstruct
surfR_inflated = '/oak/stanford/groups/leanew1/users/apines/surf/rh.inflated';
[V_R_inflated, F_R_inf] = read_surf(surfR_inflated);
F_R_inf = F_R_inf + 1;
dataR=zeros(1,5120);
dataR(incl_R)=difR;
% make figure
figure('Color','w','Position',[100 100 1600 800]);
subplot(1,2,1)
a1=trisurf(F_L_inf, V_L_inflated(:,1), V_L_inflated(:,2), V_L_inflated(:,3),'EdgeColor', 'none', 'FaceAlpha', 1);
set(a1,'FaceColor','flat','FaceVertexCData',dataL','CDataMapping','scaled');
colormap(custommap)
hold on;
caxis([mincol; maxcol]);
view([90 0]);  % Medial LH
axis equal off tight;
lighting none; camlight headlight; material dull;
% ----------- Lateral View (Right Panel) -----------
subplot(1,2,2)
a2=trisurf(F_L_inf, V_L_inflated(:,1), V_L_inflated(:,2), V_L_inflated(:,3),'EdgeColor', 'none', 'FaceAlpha', 1);
colormap(custommap)
set(a2,'FaceColor','flat','FaceVertexCData',dataL','CDataMapping','scaled');
hold on;
caxis([mincol; maxcol]);
view([-90 0]);  % Lateral LH
axis equal off tight;
lighting none; camlight headlight; material dull;
print('/scratch/users/apines/DrugGrtr_L_mag.png','-dpng','-r400');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now do right hemisphere
% make figure
figure('Color','w','Position',[100 100 1600 800]);
subplot(1,2,1)
a1=trisurf(F_R_inf, V_R_inflated(:,1), V_R_inflated(:,2), V_R_inflated(:,3),'EdgeColor', 'none', 'FaceAlpha', 1);
set(a1,'FaceColor','flat','FaceVertexCData',dataR','CDataMapping','scaled');
colormap(custommap)
caxis([mincol; maxcol]);
hold on;
view([90 0]);  % Medial LH
axis equal off tight;
lighting none; camlight headlight; material dull;

% ----------- Lateral View (Right Panel) -----------
subplot(1,2,2)
a2=trisurf(F_R_inf, V_R_inflated(:,1), V_R_inflated(:,2), V_R_inflated(:,3),'EdgeColor', 'none', 'FaceAlpha', 1);
set(a2,'FaceColor','flat','FaceVertexCData',dataR','CDataMapping','scaled');
colormap(custommap)
caxis([mincol; maxcol]);
hold on;
view([-90 0]);  % Lateral LH
axis equal off tight;
lighting none; camlight headlight; material dull;
print('/scratch/users/apines/DrugGrtr_R_mag.png','-dpng','-r400');


