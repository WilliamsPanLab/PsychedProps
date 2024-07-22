function Extract_RelativeAngles_verts_psil(subj,sesh,task)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take optical flow results, get a bottom-up and top-down resultant vector in x,y coords for each face. Measured relative to Network.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% Load in fsav4 opflow calc
childFP=['/scratch/users/apines/data/psil/' subj '/' sesh];
commonFP=['/scratch/users/apines/data/psil/'];
fpl=[commonFP '/' subj '/' sesh '/' subj '_' sesh '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
fpr=[commonFP '/' subj '/' sesh '/' subj '_' sesh '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
OpFl_L=load(fpl);
OpFl_R=load(fpr);

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
% faces_R
F_R=faces_r;
% vertices V
V_R=vx_r;
% get incenters of triangles
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
mw_L=ones(1,2562);
mw_L(mwAndTSNR_L==0)=0;
mw_R=ones(1,2562);
mw_R(mwAndTSNR_R==0)=0;

% make medial wall vector
g_noMW_combined_L=setdiff([1:2562],find(mw_L));
g_noMW_combined_R=setdiff([1:2562],find(mw_R));

% extract size of time series
vfl=OpFl_L.vertWise_Vecs_l;
vfr=OpFl_R.vertWise_Vecs_r;

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
% for each vertex
for V=1:length(g_noMW_combined_L)
	vertexInd=g_noMW_combined_L(V);
	% for each timepoint
	for fr=1:lenOpFl
		% current vector field
		relVf_L=vfl(vertexInd,fr,:);
		% xyz components
        	xComp_L=relVf_L(1);
        	yComp_L=relVf_L(2);
        	zComp_L=relVf_L(3);
		% convert to spherical coord system
        	vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(vertexInd),eld_L(vertexInd));
		% convert to spherical coordinates
       		%OpFlVec_L= [vs_L(1) vs_L(2)];
		OpF_azes_L(V,fr)=vs_L(1);
		OpF_els_L(V,fr)=vs_L(2);
	end
end

%%% RIGHT HEMISPHERE
% for each vertex
for V=1:length(g_noMW_combined_R)
        vertexInd=g_noMW_combined_R(V);
        % for each timepoint
        for fr=1:lenOpFl
                % current vector field
                relVf_R=vfr(vertexInd,fr,:);
                % xyz components
                xComp_R=relVf_R(1);
                yComp_R=relVf_R(2);
                zComp_R=relVf_R(3);
                % convert to spherical coord system
                vs_R=cart2sphvec(double([xComp_R;yComp_R;zComp_R]),azd_R(vertexInd),eld_R(vertexInd));
                % convert to spherical coordinates
                OpF_azes_R(V,fr)=vs_R(1);
                OpF_els_R(V,fr)=vs_R(2);
        end
end

% add lukas lang functions
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/lukaslang-ofd-614a2ffc50d6'))

% convert azs and els to only those within MW mask
az_L=az_L(g_noMW_combined_L);
el_L=el_L(g_noMW_combined_L);
az_R=az_R(g_noMW_combined_R);
el_R=el_R(g_noMW_combined_R);

% load in Networks
networks=load(['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/gro_Nets_fs4.mat']);
%% k = 1 to select DMN.
Dnet_LH=networks.nets.Lnets(:,1);
Dnet_RH=networks.nets.Rnets(:,1);
for k=1:4
	% get DMN indices (>.3)
	DMN_bool_L=zeros(1,2562);
	DMN_bool_R=zeros(1,2562);
	DMN_bool_L(Dnet_LH>.3)=1;
	DMN_bool_R(Dnet_RH>.3)=1;
	% get DMN inds within valid vertices
	DMN_bool_L=DMN_bool_L(g_noMW_combined_L);
	DMN_bool_R=DMN_bool_R(g_noMW_combined_R);
	DMNInds_L=find(DMN_bool_L);
	DMNInds_R=find(DMN_bool_R);
	% initialize matrix for each face over each of k=4 networks to saveout to scratch
	faceMatrix=zeros((length(g_noMW_combined_L)+length(g_noMW_combined_R)),4);
        % network of interest is DMN
        n_LH=Dnet_LH;
        n_RH=Dnet_RH;
        % calculate network gradients on sphere
        ng_L = grad(F_L, V_L, n_LH);
        ng_R = grad(F_R, V_R, n_RH);
	% convert both back to vertices for angular comparisons
	vertwise_grad_L=zeros(2562,3);
	vertwise_grad_R=zeros(2562,3);
	% for each vertex, grab adjacent face values and merge em
	for v=1:2562
        	[InvolvedFaces_l,~]=find(F_L==v);
        	vertwise_grad_L(v,:)=mean(ng_L(InvolvedFaces_l,:),1);
	        [InvolvedFaces_r,~]=find(F_R==v);
	        vertwise_grad_R(v,:)=mean(ng_R(InvolvedFaces_r,:),1);
	end
        % use medial wall mask as common starting point (from which to mask both opfl vecs and net grads further)
        ng_L=vertwise_grad_L(g_noMW_combined_L,:);
        ng_R=vertwise_grad_R(g_noMW_combined_R,:);
        % get within-DMN verts with DMN_bool
	% extract vertex-wise vector cartesian vector components
        nx_L=ng_L(DMNInds_L,1);
        ny_L=ng_L(DMNInds_L,2);
        nz_L=ng_L(DMNInds_L,3);
        nx_R=ng_R(DMNInds_R,1);
        ny_R=ng_R(DMNInds_R,2);
        nz_R=ng_R(DMNInds_R,3);

        % get spherical coordinates (az/el/r, r equal in sphere) relevant to this network
        az_L_n=az_L(DMNInds_L);
        el_L_n=el_L(DMNInds_L);
        az_R_n=az_R(DMNInds_R);
        el_R_n=el_R(DMNInds_R);

        % now same mask for the opfl vectors
        OpF_azes_L_n=OpF_azes_L(DMNInds_L,:);
        OpF_els_L_n=OpF_els_L(DMNInds_L,:);
        OpF_azes_R_n=OpF_azes_R(DMNInds_R,:);
        OpF_els_R_n=OpF_els_R(DMNInds_R,:);

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
        NangDs_L=zeros(length(DMNInds_L),lenOpFl);
        NangDs_R=zeros(length(DMNInds_R),lenOpFl);

	% get angular distance for each vertex for each timepoint
        for V=1:length(DMNInds_L);
                % get vector for each face (network vector)
                nVec=[nazes_L(V) nels_L(V)];
                % loop over each tp
                for fr=1:lenOpFl
                        % get optical flow vector
                        OpFlVec=[OpF_azes_L_n(V,fr) OpF_els_L_n(V,fr)];
                        % extract native vector for circ stats (sd)
			OpFlVec_L= [OpFlVec(1) OpFlVec(2)];
			% get angular distance at that timepoint (degrees)
                        a = acosd(min(1,max(-1, nVec(:).' *OpFlVec(:) / norm(nVec) / norm(OpFlVec) )));
                        % populate vector
                        NangDs_L(V,fr)=a;
                % end tp loop
                end
        % end each vertex loop
	end
        for V=1:length(DMNInds_R);
                % get vector for each face (network vector)
                nVec=[nazes_R(V) nels_R(V)];
                % loop over each tp
                for fr=1:lenOpFl
                        % get optical flow vector
                        OpFlVec=[OpF_azes_R_n(V,fr) OpF_els_R_n(V,fr)];
                        % extract native vector for circ stats (sd)
                        OpFlVec_R= [OpFlVec(1) OpFlVec(2)];
			% get angular distance at that timepoint (degrees)
                        a = acosd(min(1,max(-1, nVec(:).' *OpFlVec(:) / norm(nVec) / norm(OpFlVec) )));
                        % populate vector
                        NangDs_R(V,fr)=a;
                % end tp loop
                end
        % end each vertex loop
        end
        % average for this network before proceeding to next network loop
        AllAngs=[NangDs_R(:)' NangDs_L(:)'];
        % save out an AngDistMat for hierarchical distributions
	AngDist=struct;
	AngDist.Left=NangDs_L;
	AngDist.Right=NangDs_R;
	% calc outFP
	outFP=['/scratch/users/apines/data/psil/' subj '/' sesh];
	AngDistFP=[outFP '/' subj '_' sesh '_' task '_k' num2str(k) '_AngDistVertsMat.mat'];
	save(AngDistFP,'AngDist')
	% average angular distances across hemispheres
        avgD=mean(AllAngs);
        % 6/8/24: replacing with percentage for attempt at clearer presentation of results
	percBUP=length(AllAngs(AllAngs<90))/(length(AllAngs));
	Propvec=[Propvec percBUP];
        % add label
        stringVec=[stringVec ['AngD' num2str(k)]];
	% save out as csv
	T=table(Propvec','RowNames',stringVec);
	% write out
	writetable(T,[outFP '/' subj '_' sesh '_' task '_k' num2str(k) '_Prop_Feats_verts.csv'],'WriteRowNames',true)
end
