function viz_vec_fields4(subj,sesh,task)
% visualize calculated vector fields on fsaverage4 sphere

%%%%%% Set paths %%%
% grab tool dir
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% add lukas lang functions
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/lukaslang-ofd-614a2ffc50d6'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load op flo from subject %%%
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
datafp=[childfp '/' subj '_' sesh '_' task '_OpFl.mat'];
OpFl=load(datafp);
OpFl=OpFl.us;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Load TS %%%%%%%%
% load in TRs_l and TRs_r
TRs_lfp=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_L_3k.mgh'];
TRs_rfp=['/cbica/projects/pinesParcels/results/PWs/PreProc/' subj '/' subj '_AggTS_R_3k.mgh'];
% load in time series
fpL=[childfp '/' subj '_' sesh '_task-' task '_p2mm_masked_L.mgh'];
fpR=[childfp '/' subj '_' sesh '_task-' task '_p2mm_masked_R.mgh'];
dataL=MRIread(fpL).vol;
dataR=MRIread(fpR).vol;
% files to data
TRs_l=squeeze(dataL);
TRs_r=squeeze(dataR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% uptake functional data (on surface) %%%%%%
% handle input data
disp('Size of input data')
sizeInDl=size(TRs_l)
sizeInDr=size(TRs_r)
disp('Number of TRs detected L and R')
TR_n=sizeInDl(2)
sizeInDr(2)
assert(sizeInDl(2) == sizeInDr(2), 'Unequal time series length between hemispheres')
% left hemi
disp('converting left hemi to struct')
fl=struct;
% populate struct
for TRP=1:TR_n;
        fl.TRs{TRP}=TRs_l(:,TRP);
end
% r h
disp('converting right hemi to struct')
fr=struct;
for TRP=1:TR_n;
        fr.TRs{TRP}=TRs_r(:,TRP);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in fsaverage4 faces and vertices %%%
% Load in surface data
surfL = ['/oak/stanford/groups/leanew1/users/apines/surf/lh.sphere'];
surfR = ['/oak/stanford/groups/leanew1/users/apines/surf/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;

% Get incenters of triangles.
TR = TriRep(faces_l, vx_l);
P = TR.incenters;
TRr = TriRep(faces_r, vx_r);
Pr = TRr.incenters;

%%%% faces-to-vertices converter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% load in group DMN %%%%%
% load in Networks
networks=load(['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/gro_Nets_fs4.mat']);
% note 1 appears to be insulomotor for smooth solution
%% k = 1 to select DMN. k = 2 if smooth used
Dnet_LH=networks.nets.Lnets(:,1);
Dnet_RH=networks.nets.Rnets(:,1);
% calculate network gradients on sphere
PGg_L = grad(faces_l, vx_l, Dnet_LH);
PGg_R = grad(faces_r, vx_r, Dnet_RH);

% extract face-wise vector cartesian vector components
PGx_L=PGg_L(:,1);
PGy_L=PGg_L(:,2);
PGz_L=PGg_L(:,3);
PGx_R=PGg_R(:,1);
PGy_R=PGg_R(:,2);
PGz_R=PGg_R(:,3);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add TSNR mask, includes medial wall
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
mw_L=zeros(1,2562);
mw_L(mwAndTSNR_L==1)=1;
mw_R=zeros(1,2562);
mw_R(mwAndTSNR_R==1)=1;
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

% mask angular distances
goodIndr=g_noMW_combined_R;



% dividing all coordinates by scaling factor because vector size is fixed too small relative to coordinate space otherwise
scalingfactor=3;

% roy-big-bl palette imitation, inferno is just template
roybigbl_cm=inferno(16);
roybigbl_cm(1,:)=[255, 255, 0 ];
roybigbl_cm(2,:)=[255, 200, 0];
roybigbl_cm(3,:)=[255, 120, 0];
roybigbl_cm(4,:)=[255, 0, 0 ];
roybigbl_cm(5,:)=[200, 0, 0 ];
roybigbl_cm(6,:)=[150, 0, 0 ];
roybigbl_cm(7,:)=[100, 0, 0 ];
roybigbl_cm(8,:)=[60, 0, 0 ];
roybigbl_cm(9,:)=[0, 0, 80 ];
roybigbl_cm(10,:)=[0, 0, 170];
roybigbl_cm(11,:)=[75, 0, 125];
roybigbl_cm(12,:)=[125, 0, 160];
roybigbl_cm(13,:)=[75, 125, 0];
roybigbl_cm(14,:)=[0, 200, 0];
roybigbl_cm(15,:)=[0, 255, 0];
roybigbl_cm(16,:)=[0, 255, 255]; 
% scale to 1
roybigbl_cm=roybigbl_cm.*(1/255);
% interpolate color gradient
interpsteps=[0 0.0666 0.1333 0.2 0.2666 0.3333 0.4 0.4666 0.5333 0.6 0.6666 0.7333 0.8 0.86666 0.9333 1];
roybigbl_cm=interp1(interpsteps,roybigbl_cm,linspace(0,1,255));
% yellow as high
roybigbl_cm=flipud(roybigbl_cm);
% reduce just a little bit on the close-to-white coloring
roybigbl_cm=roybigbl_cm(15:240,:);



%%%%% create conversion vector: the aim is to go from included TRs to included TR pairs (valid adjacent frames).
% the main discrep. comes from discontinuous segments: if frames 1:5 are really 1,2,4,5,6, (3 had motion)
% then  pairs (1,2) (4,5) and (5,6) are analyzed. 
% so pull out non-valid TR pairs (i.e., (2,3)) and setdiff to get index of OpFl estimations w/r/t retained TRs
% the last TR of a continuous segments does not have an opfl vec field ascribed to it
for i=100:125
	u=OpFl.vf_right{i};
	vATTR=fr.TRs{i};
	% z-score
	vATTR=zscore(vATTR);
	% normalize vectors, pull em in
	figure('units','pixels','position',[0 0 8500 8500])
	axis([-30, 30, -30, 30, -30, 30]);
	% dividing all coordinates by scaling factor because vector size is fixed too small relative to coordinate space otherwise
	scalingfactor=5;
	% convert vectors to unit vectors for plotting independent of magnitude
	% thank you https://stackoverflow.com/questions/40778629/matlab-convert-vector-to-unit-vector
	ret = bsxfun(@rdivide, u, sqrt(sum(u'.^2))');
	quiver3D(Pr(g_noMW_combined_R,1)./scalingfactor,Pr(g_noMW_combined_R,2)./scalingfactor,Pr(g_noMW_combined_R,3)./scalingfactor,ret(g_noMW_combined_R,1), ret(g_noMW_combined_R,2), ret(g_noMW_combined_R,3),[.5 .5 .5],.7)
	hold on
	PGG_ret=bsxfun(@rdivide, PGg_R, sqrt(sum(PGg_R'.^2))');
	% insert 0's where nans exist
	PGG_ret(isnan(PGG_ret))=0;
	quiver3D(Pr(g_noMW_combined_R,1)./scalingfactor,Pr(g_noMW_combined_R,2)./scalingfactor,Pr(g_noMW_combined_R,3)./scalingfactor,PGG_ret(g_noMW_combined_R,1), PGG_ret(g_noMW_combined_R,2), PGG_ret(g_noMW_combined_R,3),[0.2,0.2,0.8],.7)
	hold on
	% for OpFl Vecs on PG: shrink surface ever so slightly so all vector show up clean
	trisurf(faces_r, vx_r(:, 1)/scalingfactor, vx_r(:, 2)/scalingfactor, vx_r(:, 3)/scalingfactor, vATTR, 'EdgeColor','none');
	caxis([-2.2,2.2])
	axis equal
	daspect([1, 1, 1]);
	%colormap(roybigbl_cm);
	colormap(roybigbl_cm);
	c=colorbar;
	c.FontSize=55;
	c.LineWidth=3;
	%c.Ticks=[-3 -2 -1 0 1 2 3];
	c.Location='southoutside';
	c.FontName='Arial';
	view(280,185);
	%view(60,190)
	fn=['~/boldvec_v2_' num2str(i) '_PAGoverlay.png'];
	print(fn,'-dpng')
end
	%view(60,190)
fn=['~/boldvec_v2_' num2str(i) '_PAGoverlay.png'];
print(fn,'-dpng')
end
