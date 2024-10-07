function Extract_RelativeAngles_task(subj,sesh,task)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take optical flow results, get a bottom-up and top-down resultant vector in x,y coords for each face. Measured relative to networks.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% not task is in vertex space
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% Load in fsav4 opflow calc
childfp=['/scratch/users/apines/data/psil/' subj '/' sesh ];
datafpL=[childfp '/' subj '_' sesh '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
datafpR=[childfp '/' subj '_' sesh '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
dataL=load(datafpL);
dataR=load(datafpR);
% Load in surface data
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
% vector fields
vfl=dataL.vertWise_Vecs_l;
% faces_R
F_R=faces_r;
% vertices V
V_R=vx_r;
% vector fields
vfr=dataR.vertWise_Vecs_r;
% get position of vertices
TR_L = TriRep(F_L,V_L);
P_L = TR_L.X;
TR_R = TriRep(F_R,V_R);
P_R = TR_R.X;

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
vmwIndVec_l=find(mw_L);
vmwIndVec_r=find(mw_R);
% convert to faces
F_MW_L=sum(mw_L(faces_l),2)./3;
F_MW_R=sum(mw_R(faces_r),2)./3;
% convert "partial" medial wall to medial wall
F_MW_L=ceil(F_MW_L);
F_MW_R=ceil(F_MW_R);
% face mask indices
% vertex version!
%F_MW_L=mw_L;
%F_MW_R=mw_R;
fmwIndVec_l=find(F_MW_L);
fmwIndVec_r=find(F_MW_R);
% make medial wall vector: vertices
g_noMW_combined_L=setdiff([1:2562],vmwIndVec_l);
g_noMW_combined_R=setdiff([1:2562],vmwIndVec_r);
% adaptation to task: medial wall masking to occur later
g_noMW_combined_L=1:2562;
g_noMW_combined_R=1:2562;
% save out mask for reference in python visualization script of ATS (not fs5)
save('/oak/stanford/groups/leanew1/users/apines/fs4surf/medial_wall_vectors.mat', 'g_noMW_combined_L', 'g_noMW_combined_R');

% extract size of time series
NumTRs=size(vfl);
NumTRs=NumTRs(2);
lenOpFl=NumTRs;

% initialize out dataframes
Propvec=[];
stringVec={};

% and azez and els for opflow vectors
OpF_azes_L=zeros(length(g_noMW_combined_L),lenOpFl);
OpF_els_L=zeros(length(g_noMW_combined_L),lenOpFl);
OpF_azes_R=zeros(length(g_noMW_combined_R),lenOpFl);
OpF_els_R=zeros(length(g_noMW_combined_R),lenOpFl);

% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));

% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);

% load in time series, convert angles, then loop over each network for angular distances (network gradients iteratively)
% for each face
for F=g_noMW_combined_L
	% for each timepoint
	for fr=1:lenOpFl
		% current vector field
		relVf_L=vfl(:,fr,:);
		% xyz components
        	xComp_L=relVf_L(F,1);
        	yComp_L=relVf_L(F,2);
        	zComp_L=relVf_L(F,3);
		% convert to spherical coord system
        	vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(F),eld_L(F));
		% convert to spherical coordinates
       		%OpFlVec_L= [vs_L(1) vs_L(2)];
		OpF_azes_L(F,fr)=vs_L(1);
		OpF_els_L(F,fr)=vs_L(2);
		% store in output vector (r is redundant across all vecs, only using az and el)
		%[Thetas_L(F,fr),~]=cart2pol(vs_L(1),vs_L(2));
	end
end
% mask mw out of angles
OpF_azes_L=OpF_azes_L(g_noMW_combined_L,:);
OpF_els_L=OpF_els_L(g_noMW_combined_L,:);

% for each face
for F=g_noMW_combined_R
	% for each timepoint
	for fr=1:lenOpFl
		% current vector field
		relVf_R=vfr(:,fr,:);
		% xyz components
		xComp_R=relVf_R(F,1);
		yComp_R=relVf_R(F,2);
		zComp_R=relVf_R(F,3);
		% convert to spherical coord system
		vs_R=cart2sphvec(double([xComp_R;yComp_R;zComp_R]),azd_R(F),eld_R(F));
		% convert to spherical coordinates
       		%OpFlVec_R= [vs_R(1) vs_R(2)];
		OpF_azes_R(F,fr)=vs_R(1);
		OpF_els_R(F,fr)=vs_R(2);
		% store in output vector (r is redundant across all vecs, only using az and el)
		%[Thetas_R(F,fr),~]=cart2pol(vs_R(1),vs_R(2));
	end
end
% mask mw out of angles
OpF_azes_R=OpF_azes_R(g_noMW_combined_R,:);
OpF_els_R=OpF_els_R(g_noMW_combined_R,:);

% add lukas lang functions
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/lukaslang-ofd-614a2ffc50d6'))

% convert azs and els to only those within MW mask
az_L=az_L(g_noMW_combined_L);
el_L=el_L(g_noMW_combined_L);
az_R=az_R(g_noMW_combined_R);
el_R=el_R(g_noMW_combined_R);

% load in Network
network=load(['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/Smooth_Nets_fs4.mat']);
Dnet_LH=network.nets.Lnets;
Dnet_RH=network.nets.Rnets;
for k=1
	nets_LH=network.nets.Lnets(:,k);
	nets_RH=network.nets.Rnets(:,k);
	% initialize matrix for each face 1 network to saveout to scratch
	faceMatrix=zeros((length(g_noMW_combined_L)+length(g_noMW_combined_R)),1);
        % network of interest
        n_LH=nets_LH;
        n_RH=nets_RH;
        % calculate network gradients on sphere
        ng_L = grad(F_L, V_L, n_LH);
        ng_R = grad(F_R, V_R, n_RH);
	% convert both back to vertices for angular comparisons
	vertwise_grad_L=zeros(2562,3);
	vertwise_grad_R=zeros(2562,3);
	% for each vertex, grab adjacent face values and merge em
	for v=1:2562
		[InvolvedFaces_l,~]=find(faces_l==v);
		[InvolvedFaces_r,~]=find(faces_r==v);
		vertwise_grad_L(v,:)=mean(ng_L(InvolvedFaces_l,:),1);
		vertwise_grad_R(v,:)=mean(ng_R(InvolvedFaces_r,:),1);
	end
        % adaptation for vertex-wise: medial wall masking to occur later
	InclLeft=1:2562;
	InclRight=1:2562;
        
	% extract vertex-wise vector cartesian vector components
        nx_L=vertwise_grad_L(InclLeft,1);
        ny_L=vertwise_grad_L(InclLeft,2);
        nz_L=vertwise_grad_L(InclLeft,3);
        nx_R=vertwise_grad_R(InclRight,1);
        ny_R=vertwise_grad_R(InclRight,2);
        nz_R=vertwise_grad_R(InclRight,3);

        % translate xyz spherical coordinates to az/el/r
        % get spherical coordinates (az/el/r, r equal in sphere) relevant to this network
        az_L_n=az_L(InclLeft);
        el_L_n=el_L(InclLeft);
        az_R_n=az_R(InclRight);
        el_R_n=el_R(InclRight);

        % now same mask for the opfl vectors
        OpF_azes_L_n=OpF_azes_L(InclLeft,:);
        OpF_els_L_n=OpF_els_L(InclLeft,:);
        OpF_azes_R_n=OpF_azes_R(InclRight,:);
        OpF_els_R_n=OpF_els_R(InclRight,:);

        % translate xyz vector components at coordinates to az/el/r
        nazes_L=zeros(length(az_L_n),1);
        nels_L=zeros(length(el_L_n),1);
        for i=1:length(az_L_n)
            nvs_L=cart2sphvec(double([nx_L(i);ny_L(i);nz_L(i)]),az_L_n(i),el_L_n(i));
            nazes_L(i)=nvs_L(1);
            nels_L(i)=nvs_L(2);
        end
        % right hemi
        nazes_R=zeros(length(az_R_n),1);
        nels_R=zeros(length(el_R_n),1);
        for i=1:length(az_R_n)
            nvs_R=cart2sphvec(double([nx_R(i);ny_R(i);nz_R(i)]),az_R_n(i),el_R_n(i));
            nazes_R(i)=nvs_R(1);
            nels_R(i)=nvs_R(2);
        end

        % initialize angular distance vector for each network (l and r) above
        NangDs_L=zeros(length(InclLeft),lenOpFl);
        NangDs_R=zeros(length(InclRight),lenOpFl);
	% initialize circ SD vectors
	SD_L=zeros(1,length(InclLeft));
	SD_R=zeros(1,length(InclRight));
	Thetas_L=zeros(1,lenOpFl);
	Thetas_R=zeros(1,lenOpFl);
	% get angular distance for each face for each timepoint
        for F=1:length(InclLeft);
                % get vector for each face (network vector)
                nVec=[nazes_L(F) nels_L(F)];
                % loop over each tp
                for fr=1:lenOpFl
                        % get optical flow vector
                        OpFlVec=[OpF_azes_L_n(F,fr) OpF_els_L_n(F,fr)];
                        % extract native vector for circ stats (sd)
			OpFlVec_L= [OpFlVec(1) OpFlVec(2)];
			% store in output vector (r is redundant across all vecs, only using az and el)
			[Thetas_L(fr),Mags_L(fr)]=cart2pol(OpFlVec_L(1),OpFlVec_L(2));
			% get angular distance at that timepoint (degrees)
                        a = acosd(min(1,max(-1, nVec(:).' *OpFlVec(:) / norm(nVec) / norm(OpFlVec) )));
                        % populate vector
                        NangDs_L(F,fr)=a;
                % end tp loop
                end
		% get circ SD
		L_CSD=circ_std(Thetas_L');
        	% plop into outut vector for left hemi
		SD_L(F)=L_CSD;
	% end each face loop
        end
        for F=1:length(InclRight);
                % get vector for each face (network vector)
                nVec=[nazes_R(F) nels_R(F)];
                % loop over each tp
                for fr=1:lenOpFl
                        % get optical flow vector
                        OpFlVec=[OpF_azes_R_n(F,fr) OpF_els_R_n(F,fr)];
                        % extract native vector for circ stats (sd)
                        OpFlVec_R= [OpFlVec(1) OpFlVec(2)];
                        % store in output vector (r is redundant across all vecs, only using az and el)
                        [Thetas_R(fr),Mags_R(fr)]=cart2pol(OpFlVec_R(1),OpFlVec_R(2));
			% get angular distance at that timepoint (degrees)
                        a = acosd(min(1,max(-1, nVec(:).' *OpFlVec(:) / norm(nVec) / norm(OpFlVec) )));
                        % populate vector
                        NangDs_R(F,fr)=a;
                % end tp loop
                end
		% get circ SD
                R_CSD=circ_std(Thetas_R');
                % plop into outut vector for left hemi
                SD_R(F)=R_CSD;
        % end each face loop
        end
        % average for this network before proceeding to next network loop
        AllAngs=[NangDs_R(:)' NangDs_L(:)'];
        % save out an AngDistMat for hierarchical distributions
	AngDist=struct;
	AngDist.Left=NangDs_L;
	AngDist.Right=NangDs_R;
	% calc outFP
	outFP=['/scratch/users/apines/data/psil/' subj '/' sesh];
	AngDistFP=[outFP '/' subj '_' sesh '_' task '_k' num2str(k) '_AngDistMat_task.mat'];
	save(AngDistFP,'AngDist')
end
