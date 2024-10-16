function Calc_AvgMagnitude(subj)
% read in subj-session-dose correspondence
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');
% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
addpath('/oak/stanford/groups/leanew1/users/apines/scripts/')

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
% get mesh triangles
TR_L = TriRep(F_L,V_L);
P_L = TR_L;
TR_R = TriRep(F_R,V_R);
P_R = TR_R;
% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));
% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);

% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')

% set common filepath
commonFP=['/scratch/users/apines/data/mdma/'];

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
	fp=strjoin([commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_' task '_OpFl.mat'],'');
	OpFl=load(fp);
	% need to add # or surviving TRs check to EVERY load
        survivingTrsFP=strjoin([commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        % first condition needed because some desolate scans dont even have a second column for # of valid TRs, because no sequences of TRs meeting motion criteria exists
	if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;
		OpFl_rs1.(condition).L=cat(3,OpFl.us.vf_left{:});
		OpFl_rs1.(condition).R=cat(3,OpFl.us.vf_right{:});
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
        fp=strjoin([commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_' task '_OpFl.mat'],'');
	OpFl=load(fp);
	% need to add # or surviving TRs check to EVERY load
        survivingTrsFP=strjoin([commonFP '/' subj '/' seshInfo{c} '/' subj '_' seshInfo{c} '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        % first condition needed because some desolate scans dont even have a second column for # of valid TRs, because no sequences of TRs meeting motion criteria exists
        if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;
                OpFl_rs2.(condition).L=cat(3,OpFl.us.vf_left{:});
                OpFl_rs2.(condition).R=cat(3,OpFl.us.vf_right{:});
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
		OpFl.(condition).L=cat(3, OpFl_rs1.(condition).L, OpFl_rs2.(condition).L);
	% if rs2 is 0, but not rs1, use rs1
	elseif OpFl_rs1.(condition).L(1) ~=0 && OpFl_rs2.(condition).L(1) ==0
		OpFl.(condition).L=OpFl_rs1.(condition).L;
	% if rs1 is 0 but not rs2, use rs2
	elseif OpFl_rs1.(condition).L(1)==0 && OpFl_rs2.(condition).L(1) ~=0
		OpFl.(condition).L=OpFl_rs2.(condition).L;
	end
	% right hemisphere
	% if rs1 and rs2 are ~=0, concatenate
        if OpFl_rs1.(condition).R(1) ~=0 && OpFl_rs2.(condition).R(1) ~=0
                OpFl.(condition).R=cat(3, OpFl_rs1.(condition).R, OpFl_rs2.(condition).R);
        % if rs2 is 0, but not rs1, use rs1
        elseif OpFl_rs1.(condition).R(1) ~=0 && OpFl_rs2.(condition).R(1) ==0
                OpFl.(condition).R=OpFl_rs1.(condition).R;
        % if rs1 is 0 but not rs2, use rs2
        elseif OpFl_rs1.(condition).R(1) ==0 && OpFl_rs2.(condition).R(1) ~=0
		OpFl.(condition).R=OpFl_rs2.(condition).R;
	end
            
end

% for every face, get angles in 2d instead of 3d (3 dimensions is redudnant for spherical surface)
% LEFT HEMISPHERE
for f=1:5120;
	% convert to spherical coord system for all 4 NOTE COMP"S" is plural, COMP is indexed in nested for loop
	% BV
	bv_x_comps_L=OpFl.bv.L(f,1,:);
	bv_y_comps_L=OpFl.bv.L(f,2,:);
	bv_z_comps_L=OpFl.bv.L(f,3,:);
	% PL
        pl_x_comps_L=OpFl.pl.L(f,1,:);
        pl_y_comps_L=OpFl.pl.L(f,2,:);
        pl_z_comps_L=OpFl.pl.L(f,3,:);
	% m1
        m1_x_comps_L=OpFl.m1.L(f,1,:);
        m1_y_comps_L=OpFl.m1.L(f,2,:);
        m1_z_comps_L=OpFl.m1.L(f,3,:);
	% m2
        m2_x_comps_L=OpFl.m2.L(f,1,:);
        m2_y_comps_L=OpFl.m2.L(f,2,:);
        m2_z_comps_L=OpFl.m2.L(f,3,:);
	% for each timepoint
	for tp=1:size(OpFl.bv.L,3);
		bv_x_Comp_L=bv_x_comps_L(tp);
		bv_y_Comp_L=bv_y_comps_L(tp);
		bv_z_Comp_L=bv_z_comps_L(tp);
		% convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([bv_x_Comp_L;bv_y_Comp_L;bv_z_Comp_L]),azd_L(f),eld_L(f));
		bv_vs_L(f,tp,:)=azelrho(1:2);
	end
	% placebo
	for tp=1:size(OpFl.pl.L,3);
                pl_x_Comp_L=pl_x_comps_L(tp);
                pl_y_Comp_L=pl_y_comps_L(tp);
                pl_z_Comp_L=pl_z_comps_L(tp);
                % convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([pl_x_Comp_L;pl_y_Comp_L;pl_z_Comp_L]),azd_L(f),eld_L(f));
		pl_vs_L(f,tp,:)=azelrho(1:2);
        end
	% 80 mg
	for tp=1:size(OpFl.m1.L,3);
                m1_x_Comp_L=m1_x_comps_L(tp);
                m1_y_Comp_L=m1_y_comps_L(tp);
                m1_z_Comp_L=m1_z_comps_L(tp);
                % convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([m1_x_Comp_L;m1_y_Comp_L;m1_z_Comp_L]),azd_L(f),eld_L(f));
		m1_vs_L(f,tp,:)=azelrho(1:2);
        end
	% 120 mg
	for tp=1:size(OpFl.m2.L,3);
                m2_x_Comp_L=m2_x_comps_L(tp);
                m2_y_Comp_L=m2_y_comps_L(tp);
                m2_z_Comp_L=m2_z_comps_L(tp);
                % convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
		azelrho=cart2sphvec(double([m2_x_Comp_L;m2_y_Comp_L;m2_z_Comp_L]),azd_L(f),eld_L(f));
        	m2_vs_L(f,tp,:)=azelrho(1:2);
	end
% end for every face, left
end
% RIGHT HEMISPHERE
for f=1:5120;
        % convert to spherical coord system for all 4 NOTE COMP"S" is plural, COMP is indexed in nested for loop
        % BV
        bv_x_comps_R=OpFl.bv.R(f,1,:);
        bv_y_comps_R=OpFl.bv.R(f,2,:);
        bv_z_comps_R=OpFl.bv.R(f,3,:);
        % PL
        pl_x_comps_R=OpFl.pl.R(f,1,:);
        pl_y_comps_R=OpFl.pl.R(f,2,:);
        pl_z_comps_R=OpFl.pl.R(f,3,:);
        % m1
        m1_x_comps_R=OpFl.m1.R(f,1,:);
        m1_y_comps_R=OpFl.m1.R(f,2,:);
        m1_z_comps_R=OpFl.m1.R(f,3,:);
        % m2
        m2_x_comps_R=OpFl.m2.R(f,1,:);
        m2_y_comps_R=OpFl.m2.R(f,2,:);
        m2_z_comps_R=OpFl.m2.R(f,3,:);
        % for each timepoint
        for tp=1:size(OpFl.bv.L,3)
                bv_x_Comp_R=bv_x_comps_R(tp);
                bv_y_Comp_R=bv_y_comps_R(tp);
                bv_z_Comp_R=bv_z_comps_R(tp);
                % convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
                azelrho=cart2sphvec(double([bv_x_Comp_R;bv_y_Comp_R;bv_z_Comp_R]),azd_R(f),eld_R(f));
		bv_vs_R(f,tp,:)=azelrho(1:2);
        end
        % placebo
        for tp=1:size(OpFl.pl.L,3)
                pl_x_Comp_R=pl_x_comps_R(tp);
                pl_y_Comp_R=pl_y_comps_R(tp);
                pl_z_Comp_R=pl_z_comps_R(tp);
                % convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
                azelrho=cart2sphvec(double([pl_x_Comp_R;pl_y_Comp_R;pl_z_Comp_R]),azd_R(f),eld_R(f));
		pl_vs_R(f,tp,:)=azelrho(1:2);
        end
        % 80 mg
        for tp=1:size(OpFl.m1.L,3)
                m1_x_Comp_R=m1_x_comps_R(tp);
                m1_y_Comp_R=m1_y_comps_R(tp);
                m1_z_Comp_R=m1_z_comps_R(tp);
                % convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
                azelrho=cart2sphvec(double([m1_x_Comp_R;m1_y_Comp_R;m1_z_Comp_R]),azd_R(f),eld_R(f));
		m1_vs_R(f,tp,:)=azelrho(1:2);
        end
        % 120 mg
        for tp=1:size(OpFl.m2.L,3)
                m2_x_Comp_R=m2_x_comps_R(tp);
                m2_y_Comp_R=m2_y_comps_R(tp);
                m2_z_Comp_R=m2_z_comps_R(tp);
                % convert them to x y tangent coordinates (measured in a 2D tangent plane to each point on the 3D sphere)
                azelrho=cart2sphvec(double([m2_x_Comp_R;m2_y_Comp_R;m2_z_Comp_R]),azd_R(f),eld_R(f));
		m2_vs_R(f,tp,:)=azelrho(1:2);
        end
% end for every face
end
disp('done with magnitude merging')

% one big loop that will just print out avg magnitude for each face per scan 
% initialize output vectors
BV_L=zeros(5120,1);
BV_R=zeros(5120,1);
PL_L=zeros(5120,1);
PL_R=zeros(5120,1);
M1_L=zeros(5120,1);
M1_R=zeros(5120,1);
M2_L=zeros(5120,1);
M2_R=zeros(5120,1);

% LEFT - BASELINE
% loop over each left face
for f=1:5120;
	% for each timepoint
	for tp=1:size(OpFl.bv.L,3)
		% get xy components for current timepoint
		xy=bv_vs_L(f,tp,:);
		% get magnitude
		magnitude = sqrt(sum(xy .^ 2));
		% put magnitude into output vector
		BV_L(f)=BV_L(f)+magnitude;
	end
	% average magnitude for this vertex
	BV_L(f)=BV_L(f)/size(OpFl.bv.L,3);
end

% RIGHT
% loop over each right face
for f=1:5120
	% for each timepoint
	for tp=1:size(OpFl.bv.L,3)
		% get xy components for current timepoint
		xy=bv_vs_R(f,tp,:);
		% get magnitude
		magnitude = sqrt(sum(xy .^ 2));
		% put magnitude into output vector
		BV_R(f)=BV_R(f)+magnitude;
	end
	% average magnitude for this vertex
	BV_R(f)=BV_R(f)/size(OpFl.bv.L,3);
end

% LEFT - PLACEBO
% loop over each left face
for f=1:5120
	% for each timepoint
	for tp=1:size(OpFl.pl.L,3)
		% get xy components for current timepoint
		xy=pl_vs_L(f,tp,:);
		% get magnitude
		magnitude = sqrt(sum(xy .^ 2));
		% put magnitude into output vector
		PL_L(f)=PL_L(f)+magnitude;
	end
	% average magnitude for this vertex
	PL_L(f)=PL_L(f)/size(OpFl.pl.L,3);
end

% RIGHT
% loop over each right face
for f=1:5120
	% for each timepoint
	for tp=1:size(OpFl.pl.L,3)
		% get xy components for current timepoint
		xy=pl_vs_R(f,tp,:);
		% get magnitude
		magnitude = sqrt(sum(xy .^ 2));
		% put magnitude into output vector
		PL_R(f)=PL_R(f)+magnitude;
	end
	% average magnitude for this face
	PL_R(f)=PL_R(f)/size(OpFl.pl.L,3);
end

% LEFT - 80mg
% loop over each left face
for f=1:5120;
	% for each timepoint
	for tp=1:size(OpFl.m1.L,3)
		% get xy components for current timepoint
		xy=m1_vs_L(f,tp,:);
		% get magnitude
		magnitude = sqrt(sum(xy .^ 2));
		% put magnitude into output vector
		M1_L(f)=M1_L(f)+magnitude;
	end
	% average magnitude for this face
	M1_L(f)=M1_L(f)/size(OpFl.m1.L,3);
end

% RIGHT
% loop over each right face
for f=1:5120;
	% for each timepoint
	for tp=1:size(OpFl.m1.L,3)
		% get xy components for current timepoint
		xy=m1_vs_R(f,tp,:);
		% get magnitude
		magnitude = sqrt(sum(xy .^ 2));
		% put magnitude into output vector
		M1_R(f)=M1_R(f)+magnitude;
	end
	% average magnitude for this face
	M1_R(f)=M1_R(f)/size(OpFl.m1.L,3);
end

% LEFT - 120mg
% loop over each left face
for f=1:5120;
	% for each timepoint
	for tp=1:size(OpFl.m2.L,3)
		% get xy components for current timepoint
		xy=m2_vs_L(f,tp,:);
		% get magnitude
		magnitude = sqrt(sum(xy .^ 2));
		% put magnitude into output vector
		M2_L(f)=M2_L(f)+magnitude;
	end
	% average magnitude for this face
	M2_L(f)=M2_L(f)/size(OpFl.m2.L,3);
end

% RIGHT
% loop over each right vertex
for f=1:5120
	% for each timepoint
	for tp=1:size(OpFl.m2.L,3)
		% get xy components for current timepoint
		xy=m2_vs_R(f,tp,:);
		% get magnitude
		magnitude = sqrt(sum(xy .^ 2));
		% put magnitude into output vector
		M2_R(f)=M2_R(f)+magnitude;
	end
	% average magnitude for this face
	M2_R(f)=M2_R(f)/size(OpFl.m2.L,3);
end

% save out average magnitude vectors to scratch
outFP=['/scratch/users/apines/data/mdma/' subj];
save(strjoin([outFP '/AvgMag_BV_L.mat'],""),'BV_L');
save(strjoin([outFP '/AvgMag_BV_R.mat'],""),'BV_R');
save(strjoin([outFP '/AvgMag_PL_L.mat'],""),'PL_L');
save(strjoin([outFP '/AvgMag_PL_R.mat'],""),'PL_R');
save(strjoin([outFP '/AvgMag_M1_L.mat'],""),'M1_L');
save(strjoin([outFP '/AvgMag_M1_R.mat'],""),'M1_R');
save(strjoin([outFP '/AvgMag_M2_L.mat'],""),'M2_L');
save(strjoin([outFP '/AvgMag_M2_R.mat'],""),'M2_R');

