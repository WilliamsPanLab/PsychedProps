function Extract_RelativeAngles_mice(subj,sesh,FB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take optical flow results, get a bottom-up and top-down resultant vector in x,y coords for each pixel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% Load in fsav4 opflow calc
% note this is for LSD only! Adapt recording date for ketamine if needed
childfp='/scratch/users/apines/p50_mice/proc/20200228/'
datafp=[basefp subj '_vf_out_' FB '_' num2str(sesh) '.mat']
data=load(datafp)
% load in mask
fn='/scratch/users/apines/p50_mice/proc/20200228/thy1gc6s_0p3mgkg_m2000_preLSD0p3mgkg_1/masked_dff_DS_BP_Smoothed.h5'
formask=h5read(fn, '/mask');
mask = cellfun(@(x) strcmpi(x, 'TRUE'), formask);
surf=mask;
% vector fields
vf=data.outputs.velocityFields;
% and actual signal
signalGrid=data.outputs.wvcfs;
% get incenters of pixels
TR_L = TriRep(F_L,V_L);

% USED TO get INBETWEEN SPOTS IN VIS SCRIPT
x=linspace(1, size(signalGrid,2), size(vf,2))
y=linspace(1, size(signalGrid,1), size(vf,1)

P_L = TR_L.incenters;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add TSNR mask, includes medial wall
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR_psil.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mw_L=zeros(1,2562);
mw_L(mwAndTSNR_L==1)=1;
% convert to faces
F_MW_L=sum(mw_L(faces_l),2)./3;
% convert "partial" medial wall to medial wall
F_MW_L=ceil(F_MW_L);
% face mask indices
fmwIndVec_l=find(F_MW_L);
% make medial wall vector
g_noMW_combined_L=setdiff([1:5120],fmwIndVec_l);
% USE 2D MASK INSTEAD OF MW MASK






% save out mask for reference in python visualization script of ATS (not fs5)
save('/oak/stanford/groups/leanew1/users/apines/fs5surf/medial_wall_vectors.mat', 'g_noMW_combined_L', 'g_noMW_combined_R');
% CHANGE TO FP FOR MOUSE W/ NO SUBFIELDS)
% extract size of time series
vfl=data.us.vf_left;
NumTRs=size(vfl);
NumTRs=NumTRs(2);
lenOpFl=NumTRs;
% CHANGE TO OF MOUSE VARIABLES ^
% initialize out dataframes
Propvec=[];
stringVec={};

% and azez and els for opflow vectors
OpF_azes_L=zeros(length(g_noMW_combined_L),lenOpFl);
OpF_els_L=zeros(length(g_noMW_combined_L),lenOpFl);
% NOTE THIS WILL BE ONE MAP: NOT SEP R AND L
% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
% ENSURE DIMENSIONS ARE LINED UP ^

% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);

% load in time series, convert angles, then loop over each network for angular distances (network gradients iteratively)
% for each face
for F=g_noMW_combined_L
	% for each timepoint
	for fr=1:lenOpFl
		% current vector field
		relVf_L=vfl{fr};
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

% add lukas lang functions
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/lukaslang-ofd-614a2ffc50d6'))

% convert azs and els to only those within MW mask
az_L=az_L(g_noMW_combined_L);
el_L=el_L(g_noMW_combined_L);

% load in Networks
networks=load(['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/gro_Nets_fs4.mat']);

%% k = 1 to select DMN
Dnet_LH=networks.nets.Lnets(:,1);
Dnet_RH=networks.nets.Rnets(:,1);

% smoothed dmn version
networks=load(['/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/group_Nets_fs4_Smooth.mat']);
% NEED TO CREATE SMOOTHED DMN
%% k = 1 to select DMN. k = 2 if old smooth used
Dnet_LH=networks.nets.Lnets(:,1);

for k=1:4
	nets_LH=networks.nets.Lnets(:,k);
	% create face-wise network mask
	DMN_bool_L=zeros(5120,1);
	DMN_bool_L=sum(nets_LH(faces_l),2)./3;
	DMN_bool_L(DMN_bool_L>.3)=1;
	DMN_bool_L(DMN_bool_L<.3)=0;
	DMN_bool_L=logical(DMN_bool_L);
	% combine with medial wall mask
	MasterMask_L=DMN_bool_L;
	MasterMask_L(fmwIndVec_l)=0;
	% save out for de-masking later
	writematrix(MasterMask_L,['~/MasterMask_L_' num2str(k) '.csv'])
	% initialize matrix for each face over each of k=4 networks to saveout to scratch
	faceMatrix=zeros((length(g_noMW_combined_L)+length(g_noMW_combined_R)),4);
        % network of interest
        n_LH=Dnet_LH;
        % calculate network gradients on sphere
        ng_L = grad(F_L, V_L, n_LH);
        % use medial wall mask as common starting point (from which to mask both opfl vecs and net grads further)
        ng_L=ng_L(MasterMask_L,:);
        % get NA vertices
        sumLeft=sum(ng_L,2);
        % finds 0s in left and right network gradients - this is redundant/a backup but functionally inert at the moment
        emptyLeft=find(~sumLeft);
        InclLeft=find(sumLeft);
        % note InclLeft and Right presume mw mask already applied!	
	% save InclLeft and Right to a reference .mat 
	save(['/oak/stanford/groups/leanew1/users/apines/surf/medial_wall_nullGrad' num2str(k) '_vectors.mat'], 'InclLeft', 'InclRight');

        % mask them out of medial wall mask (medial wall mask indicates what to include, emptyLeft indicates what to exclude. setdiff excludes what should be excluded (from eL) from what should be incl. (noMW)
        %n_and_g_noMW_combined_L=setdiff(g_noMW_combined_L,emptyLeft);
        %n_and_g_noMW_combined_R=setdiff(g_noMW_combined_R,emptyRight);
        % extract face-wise vector cartesian vector components

        nx_L=ng_L(InclLeft,1);
        ny_L=ng_L(InclLeft,2);
        nz_L=ng_L(InclLeft,3);

        % translate xyz spherical coordinates to az/el/r
        %[az_L,el_L,r_L]=cart2sph(P_L(n_and_g_noMW_combined_L,1),P_L(n_and_g_noMW_combined_L,2),P_L(n_and_g_noMW_combined_L,3));
        %[az_R,el_R,r_R]=cart2sph(P_R(n_and_g_noMW_combined_R,1),P_R(n_and_g_noMW_combined_R,2),P_R(n_and_g_noMW_combined_R,3));
        % get spherical coordinates (az/el/r, r equal in sphere) relevant to this network
        az_L_n=az_L(InclLeft);
        el_L_n=el_L(InclLeft);

        % now same mask for the opfl vectors
        OpF_azes_L_n=OpF_azes_L(InclLeft,:);
        OpF_els_L_n=OpF_els_L(InclLeft,:);

        % translate xyz vector components at coordinates to az/el/r
        nazes_L=zeros(length(az_L_n),1);
        nels_L=zeros(length(el_L_n),1);
        for i=1:length(az_L_n)
            nvs_L=cart2sphvec(double([nx_L(i);ny_L(i);nz_L(i)]),az_L_n(i),el_L_n(i));
            nazes_L(i)=nvs_L(1);
            nels_L(i)=nvs_L(2);
        end

        % initialize angular distance vector for each network (l and r) above
        NangDs_L=zeros(length(InclLeft),lenOpFl);
	% initialize circ SD vectors
	SD_L=zeros(1,length(InclLeft));
	Thetas_L=zeros(1,lenOpFl);
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
        % average for this network before proceeding to next network loop
        AllAngs=[NangDs_R(:)' NangDs_L(:)'];
        % average left-hemisphere values over time and plop into facematrix for this participant
        faceMatrix(InclLeft,k)=mean(NangDs_L,2);
        % and time series population
	OutTs_L=NangDs_L;
	% average angular distances across hemispheres
        avgD=mean(AllAngs);
        Propvec=[Propvec avgD];
        % add label
        stringVec=[stringVec ['AngD' num2str(k)]];
	% save out as csv
	T=table(Propvec','RowNames',stringVec);
	% calc outFP
	outFP=['/scratch/users/apines/data/psil/' subj '/' sesh];
	% write out
	writetable(T,[outFP '/' subj '_' sesh '_' task '_k' num2str(k) '_Prop_Feats_gro.csv'],'WriteRowNames',true)
	% save out faceMatrix with subject ID as csv to /scratch/users/apines/gp/PropFeatsTemp
	writematrix(faceMatrix,['/scratch/users/apines/gp/PropFeats/' subj '_' sesh '_' task '_k' num2str(k) '_faceMatrix_gro.csv'])
	% save out time series
	writematrix(OutTs_L,[outFP '/' subj '_' sesh '_' task '_k' num2str(k) '_Prop_TS_dmn_L.csv'])
end
