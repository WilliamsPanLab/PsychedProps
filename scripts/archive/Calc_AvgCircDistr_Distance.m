function Calc_AvgCircDistr_Distance(subj)
% read in subj-session-dose correspondence
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');
% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/')
% get spherical coordinates of each vertex
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
% note X returns vertices rather than incenters of faces
P_L = TR_L.X;
TR_R = TriRep(F_R,V_R);
P_R = TR_R.X;
% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));
% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);
% now we have the definition of each vertex in spherical coordinates, in degrees

% load in medial wall + SNR so we don't have to loop over every single vertex and then mask out later
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
% note this is indexing for VALID vertices as opposed to some other scripts with 1 at INVALID vertices
mw_L=ones(1,2562);
mw_L(mwAndTSNR_L==1)=0;
mw_R=ones(1,2562);
mw_R(mwAndTSNR_R==1)=0;
mw_L=logical(mw_L);
mw_R=logical(mw_R);

% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')

% set common filepath
commonFP=['/scratch/users/apines/data/mdma/'];

% initialize cross-condition distance vectors
pl_80_L=zeros(sum(mw_L),1);
pl_80_R=zeros(sum(mw_R),1);
pl_120_L=zeros(sum(mw_L),1);
pl_120_R=zeros(sum(mw_R),1);
pl_drug_L=zeros(sum(mw_L),1);
pl_drug_R=zeros(sum(mw_R),1);
bv_80_L=zeros(sum(mw_L),1);
bv_80_R=zeros(sum(mw_R),1);
bv_120_L=zeros(sum(mw_L),1);
bv_120_R=zeros(sum(mw_R),1);
bv_drug_L=zeros(sum(mw_L),1);
bv_drug_R=zeros(sum(mw_R),1);
bv_plac_L=zeros(sum(mw_L),1);
bv_plac_R=zeros(sum(mw_R),1);
m80_120_L=zeros(sum(mw_L),1);
m80_120_R=zeros(sum(mw_R),1);

% to inherit same code as group-level scripting
s=find(strcmp(subj,subjList));
% now set subj to acutal subject name for filepaths
% should be redundant with input
subj=subjList(s);
% get session-condition correspondence
seshInfo=subSeshDose{s,2:5};
% set conditions
conditions={'bv','pl','m1','m2'};
% initialize opfl structure to hold each condition/hemisphere
OpFl_rs1 = struct();
OpFl_rs2 = struct();

% for rs1
task='rs1'
% now loop over each condition to load in each and concatenate resting-state angular time series
for c=1:4
	condition=conditions{c};
	fpl=[commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	fpr=[commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	% need to add # or surviving TRs check to EVERY load
        survivingTrsFP=strjoin([commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        % first condition needed because some desolate scans dont even have a second column for # of valid TRs, because no sequences of TRs meeting motion criteria exists
	if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;
		OpFl_rs1.(condition).L=load(strjoin(fpl,'')).vertWise_Vecs_l;
		OpFl_rs1.(condition).R=load(strjoin(fpr,'')).vertWise_Vecs_r;
	else
		OpFl_rs1.(condition).L=0;
		OpFl_rs1.(condition).R=0;
	end
end

% for rs2
task='rs2'
% now loop over each condition to load in each and concatenate resting-state angular time series
for c=1:4
        condition=conditions{c};
        fpl=[commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
        fpr=[commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
        % need to add # or surviving TRs check to EVERY load
        survivingTrsFP=strjoin([commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        % first condition needed because some desolate scans dont even have a second column for # of valid TRs, because no sequences of TRs meeting motion criteria exists
        if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;
                OpFl_rs2.(condition).L=load(strjoin(fpl,'')).vertWise_Vecs_l;
                OpFl_rs2.(condition).R=load(strjoin(fpr,'')).vertWise_Vecs_r;
        else    
                OpFl_rs2.(condition).L=0;
                OpFl_rs2.(condition).R=0;
        end     
end 

% now concatenate the scans into a master OpFl struct
OpFl=struct();
for c=1:4
	% for this condition
	condition=conditions{c};	
	% if rs1 and rs2 are ~=0, concatenate
	if OpFl_rs1.(condition).L(1) ~=0 && OpFl_rs2.(condition).L(1) ~=0
		OpFl.(condition).L=cat(2, OpFl_rs1.(condition).L(mw_L,:,:), OpFl_rs2.(condition).L(mw_L,:,:));
	% if rs2 is 0, but not rs1, use rs1
	elseif OpFl_rs1.(condition).L(1) ~=0 && OpFl_rs2.(condition).L(1) ==0
		OpFl.(condition).L=OpFl_rs1.(condition).L(mw_L,:,:);
	% if rs1 is 0 but not rs2, use rs2
	elseif OpFl_rs1.(condition).L(1)==0 && OpFl_rs2.(condition).L(1) ~=0
		OpFl.(condition).L=OpFl_rs2.(condition).L(mw_L,:,:);
	end
	% right hemisphere
	% if rs1 and rs2 are ~=0, concatenate
        if OpFl_rs1.(condition).R(1) ~=0 && OpFl_rs2.(condition).R(1) ~=0
                OpFl.(condition).R=cat(2, OpFl_rs1.(condition).R(mw_R,:,:), OpFl_rs2.(condition).R(mw_R,:,:));
        % if rs2 is 0, but not rs1, use rs1
        elseif OpFl_rs1.(condition).R(1) ~=0 && OpFl_rs2.(condition).R(1) ==0
                OpFl.(condition).R=OpFl_rs1.(condition).R(mw_R,:,:);
        % if rs1 is 0 but not rs2, use rs2
        elseif OpFl_rs1.(condition).R(1) ==0 && OpFl_rs2.(condition).R(1) ~=0
		OpFl.(condition).R=OpFl_rs2.(condition).R(mw_R,:,:);
	end
            
end

% number of timepoints 
numTPsBV=size(OpFl.bv.L,2);
numTPsPL=size(OpFl.pl.L,2);
numTPsm1=size(OpFl.m1.L,2);
numTPsm2=size(OpFl.m2.L,2);

% initialize 2d vector arrays
bv_vs_L=zeros(sum(mw_L),numTPsBV,2);
bv_vs_R=zeros(sum(mw_R),numTPsBV,2);
pl_vs_L=zeros(sum(mw_L),numTPsPL,2);
pl_vs_R=zeros(sum(mw_R),numTPsPL,2);
m1_vs_L=zeros(sum(mw_L),numTPsm1,2);
m1_vs_R=zeros(sum(mw_R),numTPsm1,2);
m2_vs_L=zeros(sum(mw_L),numTPsm2,2);
m2_vs_R=zeros(sum(mw_R),numTPsm2,2);

% making an indices of valid vertices vector for indexing into within the proceeding for loop
valid_verts_L=find(mw_L);
valid_verts_R=find(mw_R);

% for every vertex, get angles in 2d instead of 3d (3 dimensions is redudnant for spherical surface)
% LEFT HEMISPHERE
for v=1:sum(mw_L);
	% convert to spherical coord system for all 4 NOTE COMP"S" is plural, COMP is indexed in nested for loop
	% BV
	bv_x_comps_L=OpFl.bv.L(v,:,1);
	bv_y_comps_L=OpFl.bv.L(v,:,2);
	bv_z_comps_L=OpFl.bv.L(v,:,3);
	% PL
        pl_x_comps_L=OpFl.pl.L(v,:,1);
        pl_y_comps_L=OpFl.pl.L(v,:,2);
        pl_z_comps_L=OpFl.pl.L(v,:,3);
	% m1
        m1_x_comps_L=OpFl.m1.L(v,:,1);
        m1_y_comps_L=OpFl.m1.L(v,:,2);
        m1_z_comps_L=OpFl.m1.L(v,:,3);
	% m2
        m2_x_comps_L=OpFl.m2.L(v,:,1);
        m2_y_comps_L=OpFl.m2.L(v,:,2);
        m2_z_comps_L=OpFl.m2.L(v,:,3);
	% for each timepoint
	for tp=1:numTPsBV
		bv_x_Comp_L=bv_x_comps_L(tp);
		bv_y_Comp_L=bv_y_comps_L(tp);
		bv_z_Comp_L=bv_z_comps_L(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([bv_x_Comp_L;bv_y_Comp_L;bv_z_Comp_L]),azd_L(valid_verts_L(v)),eld_L(valid_verts_L(v)));
		bv_vs_L(v,tp,:)=azelrho(1:2);
	end
	% placebo
	for tp=1:numTPsPL
                pl_x_Comp_L=pl_x_comps_L(tp);
                pl_y_Comp_L=pl_y_comps_L(tp);
                pl_z_Comp_L=pl_z_comps_L(tp);
                % convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([pl_x_Comp_L;pl_y_Comp_L;pl_z_Comp_L]),azd_L(valid_verts_L(v)),eld_L(valid_verts_L(v)));
		pl_vs_L(v,tp,:)=azelrho(1:2);
        end
	% 80 mg
	for tp=1:numTPsm1
                m1_x_Comp_L=m1_x_comps_L(tp);
                m1_y_Comp_L=m1_y_comps_L(tp);
                m1_z_Comp_L=m1_z_comps_L(tp);
                % convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([m1_x_Comp_L;m1_y_Comp_L;m1_z_Comp_L]),azd_L(valid_verts_L(v)),eld_L(valid_verts_L(v)));
		m1_vs_L(v,tp,:)=azelrho(1:2);
        end
	% 120 mg
	for tp=1:numTPsm2
                m2_x_Comp_L=m2_x_comps_L(tp);
                m2_y_Comp_L=m2_y_comps_L(tp);
                m2_z_Comp_L=m2_z_comps_L(tp);
                % convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([m2_x_Comp_L;m2_y_Comp_L;m2_z_Comp_L]),azd_L(valid_verts_L(v)),eld_L(valid_verts_L(v)));
        	m2_vs_L(v,tp,:)=azelrho(1:2);
	end
% end for every vertex, left
end
% RIGHT HEMISPHERE
for v=1:sum(mw_R);
        % convert to spherical coord system for all 4 NOTE COMP"S" is plural, COMP is indexed in nested for loop
        % BV
        bv_x_comps_R=OpFl.bv.R(v,:,1);
        bv_y_comps_R=OpFl.bv.R(v,:,2);
        bv_z_comps_R=OpFl.bv.R(v,:,3);
        % PL
        pl_x_comps_R=OpFl.pl.R(v,:,1);
        pl_y_comps_R=OpFl.pl.R(v,:,2);
        pl_z_comps_R=OpFl.pl.R(v,:,3);
        % m1
        m1_x_comps_R=OpFl.m1.R(v,:,1);
        m1_y_comps_R=OpFl.m1.R(v,:,2);
        m1_z_comps_R=OpFl.m1.R(v,:,3);
        % m2
        m2_x_comps_R=OpFl.m2.R(v,:,1);
        m2_y_comps_R=OpFl.m2.R(v,:,2);
        m2_z_comps_R=OpFl.m2.R(v,:,3);
        % for each timepoint
        for tp=1:numTPsBV
                bv_x_Comp_R=bv_x_comps_R(tp);
                bv_y_Comp_R=bv_y_comps_R(tp);
                bv_z_Comp_R=bv_z_comps_R(tp);
                % convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
                azelrho=cart2sphvec(double([bv_x_Comp_R;bv_y_Comp_R;bv_z_Comp_R]),azd_R(valid_verts_R(v)),eld_R(valid_verts_R(v)));
		bv_vs_R(v,tp,:)=azelrho(1:2);
        end
        % placebo
        for tp=1:numTPsPL
                pl_x_Comp_R=pl_x_comps_R(tp);
                pl_y_Comp_R=pl_y_comps_R(tp);
                pl_z_Comp_R=pl_z_comps_R(tp);
                % convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
                azelrho=cart2sphvec(double([pl_x_Comp_R;pl_y_Comp_R;pl_z_Comp_R]),azd_R(valid_verts_R(v)),eld_R(valid_verts_R(v)));
		pl_vs_R(v,tp,:)=azelrho(1:2);
        end
        % 80 mg
        for tp=1:numTPsm1
                m1_x_Comp_R=m1_x_comps_R(tp);
                m1_y_Comp_R=m1_y_comps_R(tp);
                m1_z_Comp_R=m1_z_comps_R(tp);
                % convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
                azelrho=cart2sphvec(double([m1_x_Comp_R;m1_y_Comp_R;m1_z_Comp_R]),azd_R(valid_verts_R(v)),eld_R(valid_verts_R(v)));
		m1_vs_R(v,tp,:)=azelrho(1:2);
        end
        % 120 mg
        for tp=1:numTPsm2
                m2_x_Comp_R=m2_x_comps_R(tp);
                m2_y_Comp_R=m2_y_comps_R(tp);
                m2_z_Comp_R=m2_z_comps_R(tp);
                % convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
                azelrho=cart2sphvec(double([m2_x_Comp_R;m2_y_Comp_R;m2_z_Comp_R]),azd_R(valid_verts_R(v)),eld_R(valid_verts_R(v)));
		m2_vs_R(v,tp,:)=azelrho(1:2);
        end
% end for every vertex
end
disp('done with pre-distance mapping')
% now for each vertex AGAIN: this time we'll be calculating tp x tp distance matrices for each vertex in each condition contrast of interest
% LEFT
for v=1:sum(mw_L);
	v
	% pl vs 80 for this vertex
	distanceMatrix=zeros(numTPsPL,numTPsm1);
	for tp1 = 1:numTPsPL
		for tp2 = 1:numTPsm1
			% get xy components for current timepoint
			xy_pl=pl_vs_L(v,tp1,:);
			xy_m1=m1_vs_L(v,tp2,:);
			distance = sqrt(sum((xy_pl - xy_m1) .^ 2));
			distance_matrix(tp1, tp2)=distance;
		end
	end
	% put average of this distance matrix into output vector
	pl_80_L(v)=mean(mean(distance_matrix));

	% pl vs 120 for this vertex
        distanceMatrix=zeros(numTPsPL,numTPsm2);
        for tp1 = 1:numTPsPL
                for tp2 = 1:numTPsm2
                        % get xy components for current timepoint
                        xy_pl=pl_vs_L(v,tp1,:);
                        xy_m2=m2_vs_L(v,tp2,:);
                        distance = sqrt(sum((xy_pl - xy_m2) .^ 2));
                        distance_matrix(tp1, tp2)=distance;
                end
        end

	% consider plotting one distance matrix for the records
	%figure;
	%imagesc(distance_matrix)
	%colorbar;
	%print('~/distancematrix_example.png','-dpng','-r400')	

	% put average of this distance matrix into output vector
        pl_120_L(v)=mean(mean(distance_matrix));

	% bv vs 80 for this vertex
        distanceMatrix=zeros(numTPsBV,numTPsm1);
        for tp1 = 1:numTPsBV
                for tp2 = 1:numTPsm1
                        % get xy components for current timepoint
                        xy_bv=bv_vs_L(v,tp1,:);
                        xy_m1=m1_vs_L(v,tp2,:);
                        distance = sqrt(sum((xy_bv - xy_m1) .^ 2));
                        distance_matrix(tp1, tp2)=distance;
                end
        end

	% put average of this distance matrix into output vector
        bv_80_L(v)=mean(mean(distance_matrix));

	% bv vs 120 for this vertex
        distanceMatrix=zeros(numTPsBV,numTPsm2);
        for tp1 = 1:numTPsBV
                for tp2 = 1:numTPsm2
                        % get xy components for current timepoint
                        xy_bv=bv_vs_L(v,tp1,:);
                        xy_m2=m2_vs_L(v,tp2,:);
                        distance = sqrt(sum((xy_bv - xy_m2) .^ 2));
                        distance_matrix(tp1, tp2)=distance;
                end
        end

	% put average of this distance matrix into output vector
        bv_120_L(v)=mean(mean(distance_matrix));

	% pl vs bv for this vertex
        distanceMatrix=zeros(numTPsPL,numTPsBV);
        for tp1 = 1:numTPsPL
                for tp2 = 1:numTPsBV
                        % get xy components for current timepoint
                        xy_pl=pl_vs_L(v,tp1,:);
                        xy_bv=bv_vs_L(v,tp2,:);
                        distance = sqrt(sum((xy_pl - xy_bv) .^ 2));
                        distance_matrix(tp1, tp2)=distance;
                end
        end

	% put average of this distance matrix into output vector
        pl_bv_L(v)=mean(mean(distance_matrix));

	% 80 vs 120 for this vertex
        distanceMatrix=zeros(numTPsm1,numTPsm2);
        for tp1 = 1:numTPsm1
                for tp2 = 1:numTPsm2
                        % get xy components for current timepoint
                        xy_m1=m1_vs_L(v,tp1,:);
                        xy_m2=m2_vs_L(v,tp2,:);
                        distance = sqrt(sum((xy_m1 - xy_m2) .^ 2));
                        distance_matrix(tp1, tp2)=distance;
                end
        end

	% put average of this distance matrix into output vector
        m80_120_L(v)=mean(mean(distance_matrix));
end
disp('done with left hemi distance mapping')
% RIGHT HEMISPHERE
for v=1:sum(mw_R);
	v
        % pl vs 80 for this vertex
        distanceMatrix=zeros(numTPsPL,numTPsm1);
        for tp1 = 1:numTPsPL
                for tp2 = 1:numTPsm1
                        % get xy components for current timepoint
                        xy_pl=pl_vs_R(v,tp1,:);
                        xy_m1=m1_vs_R(v,tp2,:);
                        distance = sqrt(sum((xy_pl - xy_m1) .^ 2));
                        distance_matrix(tp1, tp2)=distance;
                end
        end
        % put average of this distance matrix into output vector
        pl_80_R(v)=mean(mean(distance_matrix));

        % pl vs 120 for this vertex
        distanceMatrix=zeros(numTPsPL,numTPsm2);
        for tp1 = 1:numTPsPL
                for tp2 = 1:numTPsm2
                        % get xy components for current timepoint
                        xy_pl=pl_vs_R(v,tp1,:);
                        xy_m2=m2_vs_R(v,tp2,:);
                        distance = sqrt(sum((xy_pl - xy_m2) .^ 2));
                        distance_matrix(tp1, tp2)=distance;
                end
        end

        % put average of this distance matrix into output vector
        pl_120_R(v)=mean(mean(distance_matrix));

        % bv vs 80 for this vertex
        distanceMatrix=zeros(numTPsBV,numTPsm1);
        for tp1 = 1:numTPsBV
                for tp2 = 1:numTPsm1
                        % get xy components for current timepoint
                        xy_bv=bv_vs_R(v,tp1,:);
                        xy_m1=m1_vs_R(v,tp2,:);
                        distance = sqrt(sum((xy_bv - xy_m1) .^ 2));
                        distance_matrix(tp1, tp2)=distance;
                end
        end

        % put average of this distance matrix into output vector
        bv_80_R(v)=mean(mean(distance_matrix));

        % bv vs 120 for this vertex
        distanceMatrix=zeros(numTPsBV,numTPsm2);
        for tp1 = 1:numTPsBV
                for tp2 = 1:numTPsm2
                        % get xy components for current timepoint
                        xy_bv=bv_vs_R(v,tp1,:);
                        xy_m2=m2_vs_R(v,tp2,:);
                        distance = sqrt(sum((xy_bv - xy_m2) .^ 2));
                        distance_matrix(tp1, tp2)=distance;
                end
        end

        % put average of this distance matrix into output vector
        bv_120_R(v)=mean(mean(distance_matrix));

        % pl vs bv for this vertex
        distanceMatrix=zeros(numTPsPL,numTPsBV);
        for tp1 = 1:numTPsPL
                for tp2 = 1:numTPsBV
                        % get xy components for current timepoint
                        xy_pl=pl_vs_R(v,tp1,:);
                        xy_bv=bv_vs_R(v,tp2,:);
                        distance = sqrt(sum((xy_pl - xy_bv) .^ 2));
                        distance_matrix(tp1, tp2)=distance;
                end
        end

        % put average of this distance matrix into output vector
        pl_bv_R(v)=mean(mean(distance_matrix));

        % 80 vs 120 for this vertex
        distanceMatrix=zeros(numTPsm1,numTPsm2);
        for tp1 = 1:numTPsm1
                for tp2 = 1:numTPsm2
                        % get xy components for current timepoint
                        xy_m1=m1_vs_R(v,tp1,:);
                        xy_m2=m2_vs_R(v,tp2,:);
                        distance = sqrt(sum((xy_m1 - xy_m2) .^ 2));
                        distance_matrix(tp1, tp2)=distance;
                end
        end

        % put average of this distance matrix into output vector
        m80_120_R(v)=mean(mean(distance_matrix));
end

% average across drug conditions for pl v drug and bv v drug
pl_drug_L=(pl_80_L+pl_120_L)/2;
pl_drug_R=(pl_80_R+pl_120_R)/2;
bv_drug_L=(bv_80_L+bv_120_L)/2;
bv_drug_R=(bv_80_R+bv_120_R)/2;
sob_drug_L=(pl_80_L+bv_80_L+pl_120_L+bv_120_L)/4;
sob_drug_R=(pl_80_R+bv_80_R+pl_120_R+bv_120_R)/4;

% save out angular distance contrasts to scratch
outFP=['/scratch/users/apines/data/mdma/' subj];
save(strjoin([outFP '/pl_80_L.mat'],""),'pl_80_L');
save(strjoin([outFP '/pl_80_R.mat'],""),'pl_80_R');
save(strjoin([outFP '/pl_120_L.mat'],""),'pl_120_L');
save(strjoin([outFP '/pl_120_R.mat'],""),'pl_120_R');
save(strjoin([outFP '/bv_80_L.mat'],""),'bv_80_L');
save(strjoin([outFP '/bv_80_R.mat'],""),'bv_80_R');
save(strjoin([outFP '/bv_120_L.mat'],""),'bv_120_L');
save(strjoin([outFP '/bv_120_R.mat'],""),'bv_120_R');
save(strjoin([outFP '/pl_drug_L.mat'],""),'pl_drug_L');
save(strjoin([outFP '/pl_drug_R.mat'],""),'pl_drug_R');
save(strjoin([outFP '/bv_drug_L.mat'],""),'bv_drug_L');
save(strjoin([outFP '/bv_drug_R.mat'],""),'bv_drug_R');
save(strjoin([outFP '/pl_bv_L.mat'],""),'pl_bv_L');
save(strjoin([outFP '/pl_bv_R.mat'],""),'pl_bv_R');
save(strjoin([outFP '/sob_drug_L.mat'],""),'sob_drug_L');
save(strjoin([outFP '/sob_drug_R.mat'],""),'sob_drug_R');




