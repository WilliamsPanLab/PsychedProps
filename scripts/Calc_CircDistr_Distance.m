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

% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')

% set common filepath
commonFP=['/scratch/users/apines/data/mdma/'];

% get number of used subjects (copied and pasted from below)
numSubjs=length([1 2 3 5 7 8 9 11 12 13 14 15 16 17]);

% initialize across-subject visualizations
gplv80_Lk=zeros(2562,numSubjs);
gplv80_Rk=zeros(2562,numSubjs);
gplv120_Lk=zeros(2562,numSubjs);
gplv120_Rk=zeros(2562,numSubjs);
gplvdrug_Lk=zeros(2562,numSubjs);
gplvdrug_Rk=zeros(2562,numSubjs);
gbvdrug_Lk=zeros(2562,numSubjs);
gbvdrug_Rk=zeros(2562,numSubjs);
gbvplac_Lk=zeros(2562,numSubjs);
gbvplac_Rk=zeros(2562,numSubjs);

% subj counter for inserting into the arrays initialized above
subjCounter=1;

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
	% need to add # or surviving TRs check to EVERY load
	survivingTrsFP=strjoin([commonFP '/' subj '/' seshInfo{1} '/' subj '_' seshInfo{1} '_task-' task '_ValidSegments_Trunc.txt'],'');
	survivingTrs=load(survivingTrsFP);
	% first condition needed because some truly desolate scans dont even have a second column for # of valid TRs, because no sequences of TRs meeting motion criteria exists
	if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;
		bvopFl_L1=load(strjoin(bvfpl,''));
		bvopFl_R1=load(strjoin(bvfpr,''));
	else
		bvopFl_L1=0;
		bvopFl_R1=0;
	end
	% placebo
	pfpl=[commonFP '/' subj '/' seshInfo{2} '/' subj '_' seshInfo{2} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	pfpr=[commonFP '/' subj '/' seshInfo{2} '/' subj '_' seshInfo{2} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	% need to add # or surviving TRs check to EVERY load
        survivingTrsFP=strjoin([commonFP '/' subj '/' seshInfo{2} '/' subj '_' seshInfo{2} '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;
		popFl_L1=load(strjoin(pfpl,''));
		popFl_R1=load(strjoin(pfpr,''));
	else
		popFl_L1=0;
		popFl_R1=0;
	end
	% 80 mg
	m1fpl=[commonFP '/' subj '/' seshInfo{3} '/' subj '_' seshInfo{3} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	m1fpr=[commonFP '/' subj '/' seshInfo{3} '/' subj '_' seshInfo{3} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	% need to add # or surviving TRs check to EVERY load
        survivingTrsFP=strjoin([commonFP '/' subj '/' seshInfo{3} '/' subj '_' seshInfo{3} '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;
		m1pFl_L1=load(strjoin(m1fpl,''));
		m1pFl_R1=load(strjoin(m1fpr,''));
	else
		m1pFl_L1=0;
		m1pFl_R1=0;
        end
	% 120 mg
	m2fpl=[commonFP '/' subj '/' seshInfo{4} '/' subj '_' seshInfo{4} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	m2fpr=[commonFP '/' subj '/' seshInfo{4} '/' subj '_' seshInfo{4} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	% need to add # or surviving TRs check to EVERY load
        survivingTrsFP=strjoin([commonFP '/' subj '/' seshInfo{4} '/' subj '_' seshInfo{4} '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        if size(survivingTrs, 2) > 1 && sum(survivingTrs(:,2))>250;	
		m2pFl_L1=load(strjoin(m2fpl,''));
		m2pFl_R1=load(strjoin(m2fpr,''));
	else
		m2pFl_L1=0;
		m2pFl_R1=0;
	end
	% repeat for rs2
	task='rs2';
	% load in disributions
	bvfpl=[commonFP '/' subj '/' seshInfo{1} '/' subj '_' seshInfo{1} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	bvfpr=[commonFP '/' subj '/' seshInfo{1} '/' subj '_' seshInfo{1} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	% need to add # or surviving TRs check to EVERY load
        survivingTrsFP=strjoin([commonFP '/' subj '/' seshInfo{1} '/' subj '_' seshInfo{1} '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        if sum(survivingTrs(:,2))>250;
		bvopFl_L2=load(strjoin(bvfpl,''));
		bvopFl_R2=load(strjoin(bvfpr,''));
	else
		bvopFl_L2=0;
		bvopFl_R2=0;
	end
	% placebo
	pfpl=[commonFP '/' subj '/' seshInfo{2} '/' subj '_' seshInfo{2} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	pfpr=[commonFP '/' subj '/' seshInfo{2} '/' subj '_' seshInfo{2} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	% need to add # or surviving TRs check to EVERY load
        survivingTrsFP=strjoin([commonFP '/' subj '/' seshInfo{1} '/' subj '_' seshInfo{1} '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        if sum(survivingTrs(:,2))>250;
		popFl_L2=load(strjoin(pfpl,''));
		popFl_R2=load(strjoin(pfpr,''));
	else
		popFl_L2=0;
		bvopFl_R2=0;
	end
	% 80 mg
	m1fpl=[commonFP '/' subj '/' seshInfo{3} '/' subj '_' seshInfo{3} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	m1fpr=[commonFP '/' subj '/' seshInfo{3} '/' subj '_' seshInfo{3} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	% need to add # or surviving TRs check to EVERY load
        survivingTrsFP=strjoin([commonFP '/' subj '/' seshInfo{3} '/' subj '_' seshInfo{3} '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        if sum(survivingTrs(:,2))>250;
		m1pFl_L2=load(strjoin(m1fpl,''));
		m1pFl_R2=load(strjoin(m1fpr,''));
	else
		m1pFl_L2=0;
		m1pFl_R2=0;
	end
	% 120 mg
	m2fpl=[commonFP '/' subj '/' seshInfo{4} '/' subj '_' seshInfo{4} '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
	m2fpr=[commonFP '/' subj '/' seshInfo{4} '/' subj '_' seshInfo{4} '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
	% need to add # or surviving TRs check to EVERY load
        survivingTrsFP=strjoin([commonFP '/' subj '/' seshInfo{4} '/' subj '_' seshInfo{4} '_task-' task '_ValidSegments_Trunc.txt'],'');
        survivingTrs=load(survivingTrsFP);
        if sum(survivingTrs(:,2))>250;
		m2pFl_L2=load(strjoin(m2fpl,''));
		m2pFl_R2=load(strjoin(m2fpr,''));
	else
		m2pFl_L2=0;
		m2pFl_R2=0;
	end
	% combine distributions across resting states
	% if else conditionals are set to only rope in scans with >250 TRs remaining
	% baseline
	if isstruct(bvopFl_L1) && isstruct(bvopFl_L2)
		bvopFl_L=cat(2,bvopFl_L1.vertWise_Vecs_l,bvopFl_L2.vertWise_Vecs_l);
		bvopFl_R=cat(2,bvopFl_R1.vertWise_Vecs_r,bvopFl_R2.vertWise_Vecs_r);
	elseif isstruct(bvopFl_L1) && ~isstruct(bvopFl_L2)
		bvopFl_L=bvopFl_L1.vertWise_Vecs_l;
                bvopFl_R=bvopFl_R1.vertWise_Vecs_r;
	elseif ~isstruct(bvopFl_L1) && isstruct(bvopFl_L2)
		bvopFl_L=bvopFl_L2.vertWise_Vecs_l;
                bvopFl_R=bvopFl_R2.vertWise_Vecs_r;
	else
		bvopFl_L=0;
		bvopFl_R=0;
	end
	% placebo
	if isstruct(popFl_L1) && isstruct(popFl_L2)
		popFl_L=cat(2,popFl_L1.vertWise_Vecs_l,popFl_L2.vertWise_Vecs_l);
		popFl_R=cat(2,popFl_R1.vertWise_Vecs_r,popFl_R2.vertWise_Vecs_r);
	elseif isstruct(popFl_L1) && ~isstruct(popFl_L2)
		popFl_L=popFl_L1.vertWise_Vecs_l;
		popFl_R=popFl_R1.vertWise_Vecs_r;
	elseif ~isstruct(popFl_L1) && isstruct(popFl_L2)
		popFl_L=popFl_L2.vertWise_Vecs_l;
                popFl_R=popFl_R2.vertWise_Vecs_r;
	else
		popFl_L=0;
                popFl_R=0;
	end
	% 80 mg
	if isstruct(m1pFl_L1) && isstruct(m1pFl_L2)
		m1pFl_L=cat(2,m1pFl_L1.vertWise_Vecs_l,m1pFl_L2.vertWise_Vecs_l);
		m1pFl_R=cat(2,m1pFl_R1.vertWise_Vecs_r,m1pFl_R2.vertWise_Vecs_r);
	elseif isstruct(m1pFl_L1) && ~isstruct(m1pFl_L2)
		m1pFl_L=m1pFl_L1.vertWise_Vecs_l;
                m1pFl_R=m1pFl_R1.vertWise_Vecs_r;
	elseif ~isstruct(m1pFl_L1) && isstruct(m1pFl_L2)
		m1pFl_L=m1pFl_L2.vertWise_Vecs_l;
                m1pFl_R=m1pFl_R2.vertWise_Vecs_r;
	else
		m1pFl_L=0;
		m1pFl_R=0;
	end
	% 120 mg
	if isstruct(m2pFl_L1) && isstruct(m2pFl_L2)
		m2pFl_L=cat(2,m2pFl_L1.vertWise_Vecs_l,m2pFl_L2.vertWise_Vecs_l);
		m2pFl_R=cat(2,m2pFl_R1.vertWise_Vecs_r,m2pFl_R2.vertWise_Vecs_r);
	elseif isstruct(m2pFl_L1) && ~isstruct(m2pFl_L2)
		m2pFl_L=m2pFl_L1.vertWise_Vecs_l;
                m2pFl_R=m2pFl_R1.vertWise_Vecs_r;
	elseif ~isstruct(m2pFl_L1) && isstruct(m2pFl_L2)
		m2pFl_L=m2pFl_L2.vertWise_Vecs_l;
                m2pFl_R=m2pFl_R2.vertWise_Vecs_r;
	else
		m2pFl_L=0;
		m2pFl_R=0;
	end
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
		% if this is a valid array
		if length(bvopFl_L) ~=1
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
		else
		end
		% if this is a valid (populated) array
		if length(popFl_L) ~=1
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
		else
		end
		% if m1 is a valid array
		if length(m1pFl_L) ~=1
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
		else
		end
		% if m2pFl is a valid array
		if length(m2pFl_L) ~=1	
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
		else
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
	% make "drug" distributions which are 80 and concatenated
	% but only if both vectors are valid
	if length(m1pFl_L) ~= 1 && length(m2pFl_L) ~=1
		drug_Lthetas=cat(2,m1pFl_Lthetas,m2pFl_Lthetas);
		drug_Rthetas=cat(2,m1pFl_Rthetas,m2pFl_Rthetas);
	elseif length(m1pFl_L)==1 && length(m2pFl_L) ~=1
		drug_Lthetas=m2pFl_Lthetas;
                drug_Rthetas=m2pFl_Rthetas;
	elseif length(m1pFl_L) ~= 1 && length(m2pFl_L)==1
		drug_Lthetas=m1pFl_Lthetas;
                drug_Rthetas=m1pFl_Rthetas;	
	end
	% loop over each vertex and kuiper test dat baby
	% get sizes for conditionals 
	sizePl=size(pOpFl_Lthetas);
	sizeBV=size(bvOpFl_Lthetas);
	sizem1=size(m1pFl_Lthetas);
	sizem2=size(m2pFl_Lthetas);
	sizeD=size(drug_Lthetas);
	for v=1:2562
		% more conditionals for missing data
		if sizePl(2) > 1 && sizem1(2) > 1
			[pval plv80_Lk(v) K] = circ_kuipertest(pOpFl_Lthetas(v,:),m1pFl_Lthetas(v,:));
			[pval plv80_Rk(v) K] = circ_kuipertest(pOpFl_Rthetas(v,:),m1pFl_Rthetas(v,:));
		else
		end
		if sizePl(2) > 1 && sizem2(2) > 1
			[pval plv120_Lk(v) K] = circ_kuipertest(pOpFl_Lthetas(v,:),m2pFl_Lthetas(v,:));
			[pval plv120_Rk(v) K] = circ_kuipertest(pOpFl_Rthetas(v,:),m2pFl_Rthetas(v,:));
		else
		end
		if sizePl(2) > 1 && sizeD(2) > 1
			[pval plvdrug_Lk(v) K] = circ_kuipertest(pOpFl_Lthetas(v,:),drug_Lthetas(v,:));
			[pval plvdrug_Rk(v) K] = circ_kuipertest(pOpFl_Rthetas(v,:),drug_Rthetas(v,:));
		else
		end
		if sizeBV(2) > 1 && sizeD(2) > 1
			[pval bvdrug_Lk(v) K] = circ_kuipertest(bvOpFl_Lthetas(v,:),drug_Lthetas(v,:));
        	        [pval bvdrug_Rk(v) K] = circ_kuipertest(bvOpFl_Rthetas(v,:),drug_Rthetas(v,:));
		else
		end
		if sizeBV(2) > 1 && sizePl(2) > 1
        	        [pval bvplac_Lk(v) K] = circ_kuipertest(bvOpFl_Lthetas(v,:),pOpFl_Lthetas(v,:));
        	        [pval bvplac_Rk(v) K] = circ_kuipertest(bvOpFl_Rthetas(v,:),pOpFl_Rthetas(v,:));
		else
		end
	end
	% these same fucking if else statements again
	if sizePl(2) > 1 && sizem1(2) > 1
		% placebo vs 80
		Vis_Vertvec(plv80_Lk,plv80_Rk,strjoin([commonFP subj '_Pl_v_80_k.png'],''))
	else
	end
	if sizePl(2) > 1 && sizem2(2) > 1
		% placebo vs 120
		Vis_Vertvec(plv120_Lk,plv120_Rk,strjoin([commonFP subj '_Pl_v_120_k.png'],''))
	else
	end
	if sizePl(2) > 1 && sizeD(2) > 1
		% placebo vs drug
		Vis_Vertvec(plvdrug_Lk,plvdrug_Rk,strjoin([commonFP subj '_Pl_v_drug_k.png'],''))	
	else
	end
	if sizeBV(2) > 1 && sizeD(2) > 1
		% baseline vs drug
		Vis_Vertvec(bvdrug_Lk,bvdrug_Rk,strjoin([commonFP subj '_Bv_v_drug_k.png'],''))
	else
	end
	if sizeBV(2) > 1 && sizePl(2) > 1
		% baseline vs placebo
        	Vis_Vertvec(bvplac_Lk,bvplac_Rk,strjoin([commonFP subj '_BV_v_pl_k.png'],''))
	else
	end
	% input to across-subject arrays
	gplv80_Lk(:,subjCounter)=plv80_Lk;
	gplv80_Rk(:,subjCounter)=plv80_Rk;
	gplv120_Lk(:,subjCounter)=plv120_Lk;
	gplv120_Rk(:,subjCounter)=plv120_Rk;
	gplvdrug_Lk(:,subjCounter)=plvdrug_Lk;
	gplvdrug_Rk(:,subjCounter)=plvdrug_Rk;
	gbvdrug_Lk(:,subjCounter)=bvdrug_Lk;
	gbvdrug_Rk(:,subjCounter)=bvdrug_Rk;
	gbvplac_Lk(:,subjCounter)=bvplac_Lk;
	gbvplac_Rk(:,subjCounter)=bvplac_Rk;
	% update subjcounter
	subjCounter=subjCounter+1
end
%%%%%%%%%%%%%%
% find where group-level vectors are set to dummy value, remove from arrays prior to averaging
%%%%%%%%%%%%%%%
% save out group structures as test
save('/scratch/users/apines/gplv80_Lk.mat', 'gplv80_Lk')
save('/scratch/users/apines/gplv80_Rk.mat', 'gplv80_Rk')
save('/scratch/users/apines/gplv120_Lk.mat', 'gplv120_Lk')
save('/scratch/users/apines/gplv120_Rk.mat', 'gplv120_Rk')
save('/scratch/users/apines/gplvdrug_Lk.mat', 'gplvdrug_Lk')
save('/scratch/users/apines/gplvdrug_Rk.mat', 'gplvdrug_Rk')
save('/scratch/users/apines/gbvdrug_Lk.mat', 'gbvdrug_Lk')
save('/scratch/users/apines/gbvdrug_Rk.mat', 'gbvdrug_Rk')
save('/scratch/users/apines/gbvplac_Lk.mat', 'gbvplac_Lk')
save('/scratch/users/apines/gbvplac_Rk.mat', 'gbvplac_Rk')
% Check and remove columns with all zero values
gplv80_Lk(:, all(gplv80_Lk == 0, 1)) = [];
gplv80_Rk(:, all(gplv80_Rk == 0, 1)) = [];
gplv120_Lk(:, all(gplv120_Lk == 0, 1)) = [];
gplv120_Rk(:, all(gplv120_Rk == 0, 1)) = [];
gplvdrug_Lk(:, all(gplvdrug_Lk == 0, 1)) = [];
gplvdrug_Rk(:, all(gplvdrug_Rk == 0, 1)) = [];
gbvdrug_Lk(:, all(gbvdrug_Lk == 0, 1)) = [];
gbvdrug_Rk(:, all(gbvdrug_Rk == 0, 1)) = [];
gbvplac_Lk(:, all(gbvplac_Lk == 0, 1)) = [];
gbvplac_Rk(:, all(gbvplac_Rk == 0, 1)) = [];

% size check
disp(size(gplv80_Lk))
disp(size(gplv120_Lk))
disp(size(gplvdrug_Lk))
disp(size(gbvdrug_Lk))
disp(size(gbvplac_Lk))

% get average values across subjects
gplv80_Lk=mean(gplv80_Lk,2);
gplv80_Rk=mean(gplv80_Rk,2);
gplv120_Lk=mean(gplv120_Lk,2);
gplv120_Rk=mean(gplv120_Rk,2);
gplvdrug_Lk=mean(gplvdrug_Lk,2);
gplvdrug_Rk=mean(gplvdrug_Rk,2);
gbvdrug_Lk=mean(gbvdrug_Lk,2);
gbvdrug_Rk=mean(gbvdrug_Rk,2);
gbvplac_Lk=mean(gbvplac_Lk,2);
gbvplac_Rk=mean(gbvplac_Rk,2);

% vert surface print .pngs
Vis_Vertvec(gplv80_Lk,plv80_Rk,'~/g_Pl_v_80_k.png')
Vis_Vertvec(gplv120_Lk,plv120_Rk,'~/g_Pl_v_120_k.png')
Vis_Vertvec(gplvdrug_Lk,plvdrug_Rk,'~/g_Pl_v_drug_k.png')
Vis_Vertvec(gbvdrug_Lk,bvdrug_Rk,'~/g_Bv_v_drug_k.png')
Vis_Vertvec(gbvplac_Lk,bvplac_Rk,'~/g_BV_v_pl_k.png')
