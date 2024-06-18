% read in subj-session-dose correspondence
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');
% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))

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

% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')

% set common filepath
commonFP=['/scratch/users/apines/data/mdma/'];

% for each subject
for s=[1 2 3 5 7 8 9 11 12 13 14 15 16 17];
	% display s
	disp(s)
	% pull subject
	subj=subjList(s);
	% get session info
	seshInfo=subSeshDose{s,2:5};
	% write each session out instead of for loop for sesh's
	task='rs1';
	% load in disributions
	bvfpl=[commonFP '/' subj '/' seshInfo{1} '/' subj '_' seshInfo{1} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	bvfpr=[commonFP '/' subj '/' seshInfo{1} '/' subj '_' seshInfo{1} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	bvopFl_L1=load(strjoin(bvfpl),'');
	bvopFl_R1=load(strjoin(bvfpr),'');
	% placebo
	pfpl=[commonFP '/' subj '/' seshInfo{2} '/' subj '_' seshInfo{2} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	pfpr=[commonFP '/' subj '/' seshInfo{2} '/' subj '_' seshInfo{2} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	popFl_L1=load(strjoin(pfpl),'');
	popFl_R1=load(strjoin(pfpr),'');
	% 80 mg
	m1fpl=[commonFP '/' subj '/' seshInfo{3} '/' subj '_' seshInfo{3} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	m1fpr=[commonFP '/' subj '/' seshInfo{3} '/' subj '_' seshInfo{3} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	m1pFl_L1=load(strjoin(m1fpl),'');
	m1pFl_R1=load(strjoin(m1fpr),'');
        % 120 mg
	m2fpl=[commonFP '/' subj '/' seshInfo{4} '/' subj '_' seshInfo{4} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	m2fpr=[commonFP '/' subj '/' seshInfo{4} '/' subj '_' seshInfo{4} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	m2pFl_L1=load(strjoin(m2fpl),'');
	m2pFl_R1=load(strjoin(m2fpr),'');
	% repeat for rs2
	task='rs2';
	% load in disributions
	bvfpl=[commonFP '/' subj '/' seshInfo{1} '/' subj '_' seshInfo{1} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	bvfpr=[commonFP '/' subj '/' seshInfo{1} '/' subj '_' seshInfo{1} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	bvopFl_L2=load(strjoin(bvfpl),'');
	bvopFl_R2=load(strjoin(bvfpr),'');
	% placebo
	pfpl=[commonFP '/' subj '/' seshInfo{2} '/' subj '_' seshInfo{2} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	pfpr=[commonFP '/' subj '/' seshInfo{2} '/' subj '_' seshInfo{2} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	popFl_L2=load(strjoin(pfpl),'');
	popFl_R2=load(strjoin(pfpr),'');
	% 80 mg
	m1fpl=[commonFP '/' subj '/' seshInfo{3} '/' subj '_' seshInfo{3} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	m1fpr=[commonFP '/' subj '/' seshInfo{3} '/' subj '_' seshInfo{3} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	m1pFl_L2=load(strjoin(m1fpl),'');
	m1pFl_R2=load(strjoin(m1fpr),'');
	% 120 mg
	m2fpl=[commonFP '/' subj '/' seshInfo{4} '/' subj '_' seshInfo{4} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	m2fpr=[commonFP '/' subj '/' seshInfo{4} '/' subj '_' seshInfo{4} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	m2pFl_L2=load(strjoin(m2fpl),'');
	m2pFl_R2=load(strjoin(m2fpr),'');
	% combine distributions across resting states
	bvopFl_L=cat(2,bvopFl_L1.vertWise_Vecs_l,bvopFl_L2.vertWise_Vecs_l);
	bvopFl_R=cat(2,bvopFl_R1.vertWise_Vecs_r,bvopFl_R2.vertWise_Vecs_r);
	popFl_L=cat(2,popFl_L1.vertWise_Vecs_l,popFl_L2.vertWise_Vecs_l);
	popFl_R=cat(2,popFl_R1.vertWise_Vecs_r,popFl_R2.vertWise_Vecs_r);
	m1pFl_L=cat(2,m1pFl_L1.vertWise_Vecs_l,m1pFl_L2.vertWise_Vecs_l);
	m1pFl_R=cat(2,m1pFl_R1.vertWise_Vecs_r,m1pFl_R2.vertWise_Vecs_r);
	m2pFl_L=cat(2,m2pFl_L1.vertWise_Vecs_l,m2pFl_L2.vertWise_Vecs_l);
	m2pFl_R=cat(2,m2pFl_R1.vertWise_Vecs_r,m2pFl_R2.vertWise_Vecs_r);
	% initialize 2562 by timepoint arrays for thetas
	timepointsBV=size(bvopFl_L,2);
	timepointsPL=size(popFl_L,2);
	timepointsM1=size(m1pFl_L,2);
	timepointsM2=size(m2pFl_L,2);
	bvOpFl_Lthetas=zeros(2562,timepointsBV);
	bvOpFl_Rthetas=zeros(2562,timepointsBV);
	pOpFl_Lthetas=zeros(2562,timepointsPL);
	pOpFl_Rthetas=zeros(2562,timepointsPL);
	m1pFl_Lthetas=zeros(2562,timepointsM1);
	m1pFl_Rthetas=zeros(2562,timepointsM1);
	m2pFl_Lthetas=zeros(2562,timepointsM2);
	m2pFl_Rthetas=zeros(2562,timepointsM2);
	% need to convert cartesian vectors (xyz) to circular: circle is tangent to the surface of the sphere
	% for each vertex
 	for v=1:2562
		% for each timepoint
		for t=1:timepointsBV
			xComp_L=bvopFl_L(v,t,1);
			yComp_L=bvopFl_L(v,t,2);
			zComp_L=bvopFl_L(v,t,3);
			xComp_R=bvopFl_R(v,t,1);
			yComp_R=bvopFl_R(v,t,2);
			zComp_R=bvopFl_R(v,t,3);
			% convert to spherical coord system
			vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(v),eld_L(v));
			vs_R=cart2sphvec(double([xComp_R;yComp_R;zComp_R]),azd_R(v),eld_R(v));
			% convert to thetas
			bvOpFl_Lthetas(v,t)=cart2pol(vs_L(1),vs_L(2));
			bvOpFl_Rthetas(v,t)=cart2pol(vs_R(1),vs_R(2));
		end
		% for placebo timepoints
		for t=1:timepointsPL
			xComp_L=popFl_L(v,t,1);
			yComp_L=popFl_L(v,t,2);
			zComp_L=popFl_L(v,t,3);
			xComp_R=popFl_R(v,t,1);
			yComp_R=popFl_R(v,t,2);
			zComp_R=popFl_R(v,t,3);
			% convert to spherical coord system
			vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(v),eld_L(v));
			vs_R=cart2sphvec(double([xComp_R;yComp_R;zComp_R]),azd_R(v),eld_R(v));
			% convert to thetas
			pOpFl_Lthetas(v,t)=cart2pol(vs_L(1),vs_L(2));
			pOpFl_Rthetas(v,t)=cart2pol(vs_R(1),vs_R(2));
		end
		% for 80 mg timepoints
		for t=1:timepointsM1
			xComp_L=m1pFl_L(v,t,1);
			yComp_L=m1pFl_L(v,t,2);
			zComp_L=m1pFl_L(v,t,3);
			xComp_R=m1pFl_R(v,t,1);
			yComp_R=m1pFl_R(v,t,2);
			zComp_R=m1pFl_R(v,t,3);
			% convert to spherical coord system
			vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(v),eld_L(v));
			vs_R=cart2sphvec(double([xComp_R;yComp_R;zComp_R]),azd_R(v),eld_R(v));
			% convert to thetas
			m1pFl_Lthetas(v,t)=cart2pol(vs_L(1),vs_L(2));
			m1pFl_Rthetas(v,t)=cart2pol(vs_R(1),vs_R(2));
		end
		% for 120 mg timepoints
		for t=1:timepointsM2
			xComp_L=m2pFl_L(v,t,1);
			yComp_L=m2pFl_L(v,t,2);
			zComp_L=m2pFl_L(v,t,3);
			xComp_R=m2pFl_R(v,t,1);
			yComp_R=m2pFl_R(v,t,2);
			zComp_R=m2pFl_R(v,t,3);
			% convert to spherical coord system
			vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(v),eld_L(v));
			vs_R=cart2sphvec(double([xComp_R;yComp_R;zComp_R]),azd_R(v),eld_R(v));
			% convert to thetas
			m2pFl_Lthetas(v,t)=cart2pol(vs_L(1),vs_L(2));
			m2pFl_Rthetas(v,t)=cart2pol(vs_R(1),vs_R(2));
		end
	% end for each vertex
	end
	%%% kuiper test each vertex
	% initialize 2562 by 4 arrays for k stats
	plv80_Lk=zeros(2562,1);
	plv80_Rk=zeros(2562,1);
	plv120_Lk=zeros(2562,1);
	plv120_Rk=zeros(2562,1);
	plvdrug_Lk=zeros(2562,1);
	plvdrug_Rk=zeros(2562,1);
	bvdrug_Lk=zeros(2562,1);
	bvdrug_Rk=zeros(2562,1);
	bvplac_Lk=zeros(2562,1);
	bvplac_Rk=zeros(2562,1);
	% loop over each vertex and kuiper test dat baby
	for v=1:2562

	% pl vs. 80 
	% pl vs 120
	% pl. vs drug
	% bv vs. drug
	% bv vs. plac
	% vert surface print .pngs
% end for each subject
% average k stat over each vertex
% print comparisons
% pl vs. 80 
% pl vs 120
% pl. vs drug
% bv vs. drug
% bv vs. plac

