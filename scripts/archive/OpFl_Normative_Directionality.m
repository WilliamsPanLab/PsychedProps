% read in subject dosage correspondence, has to be before addpath for some silly reason
restoredefaultpath
subSeshDose=readtable('~/subjSeshDoseCorresp.csv');

% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfL = [SubjectsFolder '/lh.sphere'];
surfR = [SubjectsFolder '/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% get incenters of triangles
TR_L = TriRep(faces_l,vx_l);
P_L = TR_L.incenters;
TR_R = TriRep(faces_r,vx_r);
P_R = TR_R.incenters;

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

% normalize verts to unit sphere
% left
numV=length(vx_l);
vx_l(numV+1:end, :) = VecNormalize(vx_l(numV+1:end, :));
% right
numV=length(vx_r);
vx_r(numV+1:end, :) = VecNormalize(vx_r(numV+1:end, :));

% calculate normative directionality streamlines for this subject and session
n = size(faces_l, 1);
T = TriRep(faces_l, vx_l);
P = T.incenters./100;

% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));

% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);

% initialize baseline, placebo, 80 mg, and 120 mg vectors
bv_g_anglesx=zeros(17,5120);
p_g_anglesx=zeros(17,5120);
m1_g_anglesx=zeros(17,5120);
m2_g_anglesx=zeros(17,5120);
bv_g_anglesy=zeros(17,5120);
p_g_anglesy=zeros(17,5120);
m1_g_anglesy=zeros(17,5120);
m2_g_anglesy=zeros(17,5120);
bv_g_anglesz=zeros(17,5120);
p_g_anglesz=zeros(17,5120);
m1_g_anglesz=zeros(17,5120);
m2_g_anglesz=zeros(17,5120);
% within subject standard deviations
bv_g_std_th=zeros(17,5120);
p_g_std_th=zeros(17,5120);
m1_g_std_th=zeros(17,5120);
m2_g_std_th=zeros(17,5120);
% load each subject and get mean value within and across subjects at each vertex
% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')

% variable task name
task='rs2'

% initialize a vector to track length of OpFl, so I can threshold by same criteria used for main stats (>250 TRs remaining)
RemTrsVec_p=zeros(1,17);
RemTrsVec_m1=zeros(1,17);
RemTrsVec_m2=zeros(1,17);
% for each subject
for s=[1 2 3 5 7 8 9 11 12 13 14 15 16 17]
	s
	% get session info
        seshInfo=subSeshDose{s,2:5};
        %% ADD IN REST OF FILEPATHS
        bvFP=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{1} '/' subjList(s) '_' seshInfo{1} '_' task '_OpFl.mat'];
        bvFP=strjoin(bvFP,'');
	pFP=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{2} '/' subjList(s) '_' seshInfo{2} '_' task '_OpFl.mat'];
	pFP=strjoin(pFP,'');
	m1FP=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{3} '/' subjList(s) '_' seshInfo{3} '_' task '_OpFl.mat'];
        m1FP=strjoin(m1FP,'');
	m2FP=['/scratch/users/apines/data/mdma/' subjList(s) '/' seshInfo{4} '/' subjList(s) '_' seshInfo{4} '_' task '_OpFl.mat'];
        m2FP=strjoin(m2FP,'');

	% calculate normative directionality streamlines for this subject and session
	n = size(faces_l, 1);
	T = TriRep(faces_l, vx_l);
	P = T.incenters./100;

	% if bv exists
	if exist(bvFP,'file')
	bv=load(bvFP);
	%%% bv as left Hemi
	bv=bv.us.vf_left;
	NumTRs=size(bv);
	NumTRs=NumTRs(2);
	lenOpFl=NumTRs;
	% save thetas a rho
	Thetas_L=zeros(lenOpFl,5120);
	Rhos_L=zeros(lenOpFl,5120);
	% loop over each face
	for F=g_noMW_combined_L
		% initialize angles for this face over time
		x_Face=zeros(1,lenOpFl);
		y_Face=zeros(1,lenOpFl);
		z_Face=zeros(1,lenOpFl);
		% loop over each frame to extract from unfortunate cell structure
		for fr=1:lenOpFl
			% current vector field
	        	relVf_L=bv{fr};
			% xyz components
	        	x_Face(fr)=relVf_L(F,1);
	        	y_Face(fr)=relVf_L(F,2);
	        	z_Face(fr)=relVf_L(F,3);
			% convert to circular for circular STD
			vs_L=cart2sphvec(double([x_Face(fr);y_Face(fr);z_Face(fr)]),azd_L(F),eld_L(F));
			[Thetas_L(fr,F),Rhos_L(fr,F)]=cart2pol(vs_L(1),vs_L(2));
		end
		% get average angle over time
		bv_g_anglesx(s,F)=mean(x_Face);
		bv_g_anglesy(s,F)=mean(y_Face);
		bv_g_anglesz(s,F)=mean(z_Face);
		bv_g_std_th(s,F)=circ_std(Thetas_L(:,F));
	end
	else
	end

	%%% p as left Hemi (if it exists)
	if exist(pFP,'file')
	p=load(pFP);
	p=p.us.vf_left;
	NumTRs=size(p);
	NumTRs=NumTRs(2);
	lenOpFl=NumTRs;
	RemTrsVec_p(s)=lenOpFl;
	% rhos and thetas
	Thetas_L=zeros(lenOpFl,5120);
	Rhos_L=zeros(lenOpFl,5120);
	% loop over each face
	for F=g_noMW_combined_L
		% initialize angles for this face over time
                x_Face=zeros(1,lenOpFl);
                y_Face=zeros(1,lenOpFl);
                z_Face=zeros(1,lenOpFl);
		% loop over each frame to extract from unfortunate cell structure
		for fr=1:lenOpFl
			% current vector field
	        	relVf_L=p{fr};
			% xyz components
	        	x_Face(fr)=relVf_L(F,1);
	        	y_Face(fr)=relVf_L(F,2);
	        	z_Face(fr)=relVf_L(F,3);
			% convert to circular for circular STD
                        vs_L=cart2sphvec(double([x_Face(fr);y_Face(fr);z_Face(fr)]),azd_L(F),eld_L(F));
                        [Thetas_L(fr,F),Rhos_L(fr,F)]=cart2pol(vs_L(1),vs_L(2));
		end
		% get average angle over time
                p_g_anglesx(s,F)=mean(x_Face);
                p_g_anglesy(s,F)=mean(y_Face);
                p_g_anglesz(s,F)=mean(z_Face);
		p_g_std_th(s,F)=circ_std(Thetas_L(:,F));
	end
	else
	end
	
	%%% m1 as left Hemi
	if exist(m1FP,'file')
	m1=load(m1FP);
	m1=m1.us.vf_left;
	NumTRs=size(m1);
	NumTRs=NumTRs(2);
	lenOpFl=NumTRs;
	RemTrsVec_m1(s)=lenOpFl;
	% rhos and thetas
	Thetas_L=zeros(lenOpFl,5120);
	Rhos_L=zeros(lenOpFl,5120);
	% loop over each face
	for F=g_noMW_combined_L
		% initialize angles for this face over time
                x_Face=zeros(1,lenOpFl);
                y_Face=zeros(1,lenOpFl);
                z_Face=zeros(1,lenOpFl);
		% loop over each frame to extract from unfortunate cell structure
		for fr=1:lenOpFl
			% current vector field
	        	relVf_L=m1{fr};
			% xyz components
			x_Face(fr)=relVf_L(F,1);
                        y_Face(fr)=relVf_L(F,2);
                        z_Face(fr)=relVf_L(F,3);
			% convert to circular for circular STD
                        vs_L=cart2sphvec(double([x_Face(fr);y_Face(fr);z_Face(fr)]),azd_L(F),eld_L(F));
                        [Thetas_L(fr,F),Rhos_L(fr,F)]=cart2pol(vs_L(1),vs_L(2));
		end
		% get average angle over time
                m1_g_anglesx(s,F)=mean(x_Face);
                m1_g_anglesy(s,F)=mean(y_Face);
                m1_g_anglesz(s,F)=mean(z_Face);
		m1_g_std_th(s,F)=circ_std(Thetas_L(:,F));
	end
	else
	end
	%%% m2 as left Hemi
	if exist(m2FP,'file')
	m2=load(m2FP);
	m2=m2.us.vf_left;
	NumTRs=size(m2);
	NumTRs=NumTRs(2);
	lenOpFl=NumTRs;
	RemTrsVec_m2(s)=lenOpFl;
	% rhos and thetas
	Thetas_L=zeros(lenOpFl,5120);
	Rhos_L=zeros(lenOpFl,5120);
	% loop over each face
	for F=g_noMW_combined_L
		% initialize angles for this face over time
                x_Face=zeros(1,lenOpFl);
                y_Face=zeros(1,lenOpFl);
                z_Face=zeros(1,lenOpFl);
		% loop over each frame to extract from unfortunate cell structure
		for fr=1:lenOpFl
			% current vector field
	        	relVf_L=m2{fr};
			% xyz components
                        x_Face(fr)=relVf_L(F,1);
                        y_Face(fr)=relVf_L(F,2);
                        z_Face(fr)=relVf_L(F,3);
			% convert to circular for circular STD
                        vs_L=cart2sphvec(double([x_Face(fr);y_Face(fr);z_Face(fr)]),azd_L(F),eld_L(F));
                        [Thetas_L(fr,F),Rhos_L(fr,F)]=cart2pol(vs_L(1),vs_L(2));
		end
		% get average angle over time
                m2_g_anglesx(s,F)=mean(x_Face);
                m2_g_anglesy(s,F)=mean(y_Face);
                m2_g_anglesz(s,F)=mean(z_Face);
		m2_g_std_th(s,F)=circ_std(Thetas_L(:,F));
	end
	else
	end
end

% now threhsold by remaining TRs
ThreshVal=250;
passes_p=RemTrsVec_p>ThreshVal;
passes_m1=RemTrsVec_m1>ThreshVal;
passes_m2=RemTrsVec_m2>ThreshVal;
% apply threshold
p_g_anglesx=p_g_anglesx(passes_p,:);
p_g_anglesy=p_g_anglesy(passes_p,:);
p_g_anglesz=p_g_anglesz(passes_p,:);
m1_g_anglesx=m1_g_anglesx(passes_m1,:);
m1_g_anglesy=m1_g_anglesy(passes_m1,:);
m1_g_anglesz=m1_g_anglesz(passes_m1,:);
m2_g_anglesx=m2_g_anglesx(passes_m2,:);
m2_g_anglesy=m2_g_anglesy(passes_m2,:);
m2_g_anglesz=m2_g_anglesz(passes_m2,:);

% same to standard deviation
p_g_std_th=p_g_std_th(passes_p,:);
m1_g_std_th=m1_g_std_th(passes_m1,:);
m2_g_std_th=m2_g_std_th(passes_m2,:);

% remove blank rows (should be 4, 6, 10)
%bv_g_anglesx = bv_g_anglesx(~all(bv_g_anglesx == 0, 2), :);
%bv_g_anglesy = bv_g_anglesy(~all(bv_g_anglesy == 0, 2), :);
%bv_g_anglesz = bv_g_anglesz(~all(bv_g_anglesz == 0, 2), :);
%p_g_anglesx = p_g_anglesx(~all(p_g_anglesx == 0, 2), :);
%p_g_anglesy = p_g_anglesy(~all(p_g_anglesy == 0, 2), :);
%p_g_anglesz = p_g_anglesz(~all(p_g_anglesz == 0, 2), :);
%m1_g_anglesx = m1_g_anglesx(~all(m1_g_anglesx == 0, 2), :);
%m1_g_anglesy = m1_g_anglesy(~all(m1_g_anglesy == 0, 2), :);
%m1_g_anglesz = m1_g_anglesz(~all(m1_g_anglesz == 0, 2), :);
%m2_g_anglesx = m2_g_anglesx(~all(m2_g_anglesx == 0, 2), :);
%m2_g_anglesy = m2_g_anglesy(~all(m2_g_anglesy == 0, 2), :);
%m2_g_anglesz = m2_g_anglesz(~all(m2_g_anglesz == 0, 2), :);
% apply same to standard deviations
%bv_g_std_th = bv_g_std_th(~all(bv_g_std_th == 0, 2), :);
%p_g_std_th = p_g_std_th(~all(p_g_std_th == 0, 2), :);
%m1_g_std_th = m1_g_std_th(~all(m1_g_std_th == 0, 2), :);
%m2_g_std_th = m2_g_std_th(~all(m2_g_std_th == 0, 2), :);

% get means across angles across each participant
bv_g_angles_meanx=mean(bv_g_anglesx);
p_g_angles_meanx=mean(p_g_anglesx);
m1_g_angles_meanx=mean(m1_g_anglesx);
m2_g_angles_meanx=mean(m2_g_anglesx);
% y components
bv_g_angles_meany=mean(bv_g_anglesy);
p_g_angles_meany=mean(p_g_anglesy);
m1_g_angles_meany=mean(m1_g_anglesy);
m2_g_angles_meany=mean(m2_g_anglesy);
% z components
bv_g_angles_meanz=mean(bv_g_anglesz);
p_g_angles_meanz=mean(p_g_anglesz);
m1_g_angles_meanz=mean(m1_g_anglesz);
m2_g_angles_meanz=mean(m2_g_anglesz);

% initialize across-subject mean angle thetas (radians) by session
sizeP=size(p_g_anglesx);
sizem1=size(m1_g_anglesx);
sizem2=size(m2_g_anglesx);
cross_subjs_means_p=zeros(sizeP(1),5120);
cross_subjs_means_m1=zeros(sizem1(1),5120);
cross_subjs_means_m2=zeros(sizem2(1),5120);

% and output scalar
cross_subjs_std_p=zeros(1,5120);
cross_subjs_std_m1=zeros(1,5120);
cross_subjs_std_m2=zeros(1,5120);

% get across subject SD
for F=g_noMW_combined_L
	subjMeanX_p=p_g_anglesx(:,F);
	subjMeanY_p=p_g_anglesy(:,F);
	subjMeanZ_p=p_g_anglesz(:,F);
	subjMeanX_m1=m1_g_anglesx(:,F);
        subjMeanY_m1=m1_g_anglesy(:,F);
        subjMeanZ_m1=m1_g_anglesz(:,F);
	subjMeanX_m2=m2_g_anglesx(:,F);
        subjMeanY_m2=m2_g_anglesy(:,F);
        subjMeanZ_m2=m2_g_anglesz(:,F);
	% get size of remaining observations
	p_obs=size(subjMeanX_p);
	p_obs=p_obs(1);
	m1_obs=size(subjMeanX_m1);
        m1_obs=m1_obs(1);
	m2_obs=size(subjMeanX_m2);
        m2_obs=m2_obs(1);
	% for each remaining S
	for S=1:p_obs
		% convert each to radians
		vs_L_p=cart2sphvec(double([subjMeanX_p(S);subjMeanY_p(S);subjMeanZ_p(S)]),azd_L(F),eld_L(F));
		[cross_subjs_means_p(S,F), temp ]=cart2pol(vs_L_p(1),vs_L_p(2));
	end
	% 80 mg
	for S=1:m1_obs
                % convert each to radians
                vs_L_m1=cart2sphvec(double([subjMeanX_m1(S);subjMeanY_m1(S);subjMeanZ_m1(S)]),azd_L(F),eld_L(F));
                [cross_subjs_means_m1(S,F), temp ]=cart2pol(vs_L_m1(1),vs_L_m1(2));
        end
	% 120 mg
	for S=1:m2_obs
                % convert each to radians
                vs_L_m2=cart2sphvec(double([subjMeanX_m2(S);subjMeanY_m2(S);subjMeanZ_m2(S)]),azd_L(F),eld_L(F));
                [cross_subjs_means_m2(S,F), temp ]=cart2pol(vs_L_m2(1),vs_L_m2(2));
        end
	% get std of this face
	cross_subjs_std_p(F)=circ_std(cross_subjs_means_p(:,F));
	cross_subjs_std_m1(F)=circ_std(cross_subjs_means_m1(:,F));
	cross_subjs_std_m2(F)=circ_std(cross_subjs_means_m2(:,F));	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% orthogonalize to the surface of inflated
% inflated
surfLi = ['/oak/stanford/groups/leanew1/users/apines/surf/lh.inflated'];
surfRi = ['/oak/stanford/groups/leanew1/users/apines/surf/rh.inflated'];
% surface topography
[vx_li, faces_li] = read_surf(surfLi);
[vx_ri, faces_ri] = read_surf(surfRi);
% +1 the faces: begins indexing at 0
faces_li = faces_li + 1;
faces_ri = faces_ri + 1;
% Get incenters of triangles.
TriR_L = TriRep(faces_li, vx_li);
L_inc = TriR_L.incenters;
TriR_R = TriRep(faces_ri, vx_ri);
R_inc = TriR_R.incenters;
% initialize flattened vectors
flattenedVectors_L=zeros(2562,3);
flattenedVectors_L_m1=zeros(2562,3);
flattenedVectors_L_m2=zeros(2562,3);
flattenedSTD_Win_p_L=zeros(2562,1);
flattenedSTD_Btw_p_L=zeros(2562,1);
flattenedSTD_Win_m1_L=zeros(2562,1);
flattenedSTD_Btw_m1_L=zeros(2562,1);
flattenedSTD_Win_m2_L=zeros(2562,1);
flattenedSTD_Btw_m2_L=zeros(2562,1);
% change OG vecs for each condition. First = placebo
OGVecs_L=[p_g_angles_meanx; p_g_angles_meany; p_g_angles_meanz]';
OGVecs_L_m1=[m1_g_angles_meanx; m1_g_angles_meany; m1_g_angles_meanz]';
OGVecs_L_m2=[m2_g_angles_meanx; m2_g_angles_meany; m2_g_angles_meanz]';
% for each vertex (this one is for placebo and the vector STDs
for vertInd=1:2562
	% get faces intertwined with this vertex
	[InvolvedFaces_l,~]=find(faces_l==vertInd);
	% extract valid faces from this vector
	InvolvedFaces_l=intersect(InvolvedFaces_l,g_noMW_combined_L);	
	% get mean vector directionality at these faces (x y z)
	OGVec_L=mean(OGVecs_L(InvolvedFaces_l,:),1);
	% get mean SDs at these faces	
	flattenedSTD_Win_p_L(vertInd)=mean(mean(p_g_std_th(:,InvolvedFaces_l)));
	flattenedSTD_Win_m1_L(vertInd)=mean(mean(m1_g_std_th(:,InvolvedFaces_l)));
	flattenedSTD_Win_m2_L(vertInd)=mean(mean(m2_g_std_th(:,InvolvedFaces_l)));
	flattenedSTD_Btw_p_L(vertInd)=mean(cross_subjs_std_p(InvolvedFaces_l));
	flattenedSTD_Btw_m1_L(vertInd)=mean(cross_subjs_std_m1(InvolvedFaces_l));
	flattenedSTD_Btw_m2_L(vertInd)=mean(cross_subjs_std_m2(InvolvedFaces_l));	
	% initialize projection vector (directionality and magnitude w/r/t neighboring faces)
        projection_L=zeros(5,3);
	faceCounter=1;
        % construe vector in terms of how much it points to each neighboring face
        for neighborFace = InvolvedFaces_l'
        	% Calculate normal vector for the neighboring face
                normalVectors = cross(vx_l(faces_l(neighborFace, 2), :) - vx_l(faces_l(neighborFace, 1), :), ...
                                vx_l(faces_l(neighborFace, 3), :) - vx_l(faces_l(neighborFace, 1), :));
               	meanNormalVector = mean(normalVectors, 1);
                meanNormalVector = VecNormalize(meanNormalVector);

                % Project the original vector onto the normal vector
                projection_L(faceCounter, :) = dot(OGVec_L, meanNormalVector) * meanNormalVector;
                faceCounter = faceCounter + 1;
        end
	% Reconstruct the vector on the inflated surface (new surface)
        % Using vx_li and vx_ri for new surface vertices
        reconstructedVec_L = zeros(1, 3);
        for i = 1:length(InvolvedFaces_l)
        	neighborFace = InvolvedFaces_l(i);
        	newVertexCoords_L = vx_li(faces_li(neighborFace, :), :);
        	meanNewVertex_L = mean(newVertexCoords_L, 1); % Get the mean position of the new face vertices
        	reconstructedVec_L = reconstructedVec_L + projection_L(i, :) .* meanNewVertex_L; % Element-wise multiplication
        end	
	% Remove orthogonal component to flatten to the inflated surface
        involvedFaces_L = find(any(faces_li == vertInd, 2)); % Find faces involving this vertex
        normalVectors_L = cross(vx_li(faces_li(involvedFaces_L, 2), :) - vx_li(faces_li(involvedFaces_L, 1), :), ...
                        vx_li(faces_li(involvedFaces_L, 3), :) - vx_li(faces_li(involvedFaces_L, 1), :));
        meanNormalVector_L = mean(normalVectors_L, 1);
        meanNormalVector_L = VecNormalize(meanNormalVector_L);

        OGvecOrthogonal_L = dot(reconstructedVec_L, meanNormalVector_L) * meanNormalVector_L;
        modVec_L = reconstructedVec_L - OGvecOrthogonal_L;
        modVec_L = VecNormalize(modVec_L);
	flattenedVectors_L(vertInd, :) = modVec_L;
end
% for m1
% for each vertex (this one is for placebo and the vector STDs
for vertInd=1:2562
        % get faces intertwined with this vertex
        [InvolvedFaces_l,~]=find(faces_l==vertInd);
        % extract valid faces from this vector
        InvolvedFaces_l=intersect(InvolvedFaces_l,g_noMW_combined_L);
        % get mean vector directionality at these faces (x y z)
        OGVec_L=mean(OGVecs_L_m1(InvolvedFaces_l,:),1);
        % initialize projection vector (directionality and magnitude w/r/t neighboring faces)
        projection_L=zeros(5,3);
        faceCounter=1;
        % construe vector in terms of how much it points to each neighboring face
        for neighborFace = InvolvedFaces_l'
                % Calculate normal vector for the neighboring face
                normalVectors = cross(vx_l(faces_l(neighborFace, 2), :) - vx_l(faces_l(neighborFace, 1), :), ...
                                vx_l(faces_l(neighborFace, 3), :) - vx_l(faces_l(neighborFace, 1), :));
                meanNormalVector = mean(normalVectors, 1);
                meanNormalVector = VecNormalize(meanNormalVector);

                % Project the original vector onto the normal vector
                projection_L(faceCounter, :) = dot(OGVec_L, meanNormalVector) * meanNormalVector;
                faceCounter = faceCounter + 1;
        end
        % Reconstruct the vector on the inflated surface (new surface)
        % Using vx_li and vx_ri for new surface vertices
        reconstructedVec_L = zeros(1, 3);
        for i = 1:length(InvolvedFaces_l)
                neighborFace = InvolvedFaces_l(i);
                newVertexCoords_L = vx_li(faces_li(neighborFace, :), :);
                meanNewVertex_L = mean(newVertexCoords_L, 1); % Get the mean position of the new face vertices
                reconstructedVec_L = reconstructedVec_L + projection_L(i, :) .* meanNewVertex_L; % Element-wise multiplication
        end
        % Remove orthogonal component to flatten to the inflated surface
        involvedFaces_L = find(any(faces_li == vertInd, 2)); % Find faces involving this vertex
        normalVectors_L = cross(vx_li(faces_li(involvedFaces_L, 2), :) - vx_li(faces_li(involvedFaces_L, 1), :), ...
                        vx_li(faces_li(involvedFaces_L, 3), :) - vx_li(faces_li(involvedFaces_L, 1), :));
        meanNormalVector_L = mean(normalVectors_L, 1);
        meanNormalVector_L = VecNormalize(meanNormalVector_L);

        OGvecOrthogonal_L = dot(reconstructedVec_L, meanNormalVector_L) * meanNormalVector_L;
        modVec_L = reconstructedVec_L - OGvecOrthogonal_L;
        modVec_L = VecNormalize(modVec_L);
        flattenedVectors_L_m1(vertInd, :) = modVec_L;
end
% for m2
% for each vertex (this one is for placebo and the vector STDs
for vertInd=1:2562
        % get faces intertwined with this vertex
        [InvolvedFaces_l,~]=find(faces_l==vertInd);
        % extract valid faces from this vector
        InvolvedFaces_l=intersect(InvolvedFaces_l,g_noMW_combined_L);
        % get mean vector directionality at these faces (x y z)
        OGVec_L=mean(OGVecs_L_m2(InvolvedFaces_l,:),1);
        % initialize projection vector (directionality and magnitude w/r/t neighboring faces)
        projection_L=zeros(5,3);
        faceCounter=1;
        % construe vector in terms of how much it points to each neighboring face
        for neighborFace = InvolvedFaces_l'
                % Calculate normal vector for the neighboring face
                normalVectors = cross(vx_l(faces_l(neighborFace, 2), :) - vx_l(faces_l(neighborFace, 1), :), ...
                                vx_l(faces_l(neighborFace, 3), :) - vx_l(faces_l(neighborFace, 1), :));
                meanNormalVector = mean(normalVectors, 1);
                meanNormalVector = VecNormalize(meanNormalVector);

                % Project the original vector onto the normal vector
                projection_L(faceCounter, :) = dot(OGVec_L, meanNormalVector) * meanNormalVector;
                faceCounter = faceCounter + 1;
        end
        % Reconstruct the vector on the inflated surface (new surface)
        % Using vx_li and vx_ri for new surface vertices
        reconstructedVec_L = zeros(1, 3);
        for i = 1:length(InvolvedFaces_l)
                neighborFace = InvolvedFaces_l(i);
                newVertexCoords_L = vx_li(faces_li(neighborFace, :), :);
                meanNewVertex_L = mean(newVertexCoords_L, 1); % Get the mean position of the new face vertices
                reconstructedVec_L = reconstructedVec_L + projection_L(i, :) .* meanNewVertex_L; % Element-wise multiplication
        end
        % Remove orthogonal component to flatten to the inflated surface
        involvedFaces_L = find(any(faces_li == vertInd, 2)); % Find faces involving this vertex
        normalVectors_L = cross(vx_li(faces_li(involvedFaces_L, 2), :) - vx_li(faces_li(involvedFaces_L, 1), :), ...
                        vx_li(faces_li(involvedFaces_L, 3), :) - vx_li(faces_li(involvedFaces_L, 1), :));
        meanNormalVector_L = mean(normalVectors_L, 1);
        meanNormalVector_L = VecNormalize(meanNormalVector_L);

        OGvecOrthogonal_L = dot(reconstructedVec_L, meanNormalVector_L) * meanNormalVector_L;
        modVec_L = reconstructedVec_L - OGvecOrthogonal_L;
        modVec_L = VecNormalize(modVec_L);
        flattenedVectors_L_m2(vertInd, :) = modVec_L;
end

%%%%%%%%%%%%%%%%%%%%%%% plotting

% set right hemisphere vect rs as null
flattenedVectors_R=zeros(2562,3);
% read in DMN as background
networks=load(['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/gro_Nets_fs4.mat']);
nets_LH=networks.nets.Lnets(:,1);
nets_RH=networks.nets.Rnets(:,1);
% binarized DMN mask
% print em out: p m1 m2 by normative directionality color, then colored by within subj SD and between subj SD
fn='~/placeboVecs_directions.png'
Vis_Surf_n_Vecfield_NormativeDir(nets_LH,nets_RH,flattenedVectors_L,flattenedVectors_R,fn,'Directional');
% replace nans with 0s
flattenedSTD_Win_p_L(isnan(flattenedSTD_Win_p_L))=0;
fn='~/placeboVecs_WinSTD.png'
Vis_Surf_n_Vecfield_NormativeDir(flattenedSTD_Win_p_L,nets_RH,flattenedVectors_L,flattenedVectors_R,fn,'STD');
% replace nans with 0s
flattenedSTD_Btw_p_L(isnan(flattenedSTD_Btw_p_L))=0;
fn='~/placeboVecs_BWSTD.png'
Vis_Surf_n_Vecfield_NormativeDir(flattenedSTD_Btw_p_L,nets_RH,flattenedVectors_L,flattenedVectors_R,fn,'STD');
% m1
fn='~/m80_directions.png'
Vis_Surf_n_Vecfield_NormativeDir(nets_LH,nets_RH,flattenedVectors_L_m1,flattenedVectors_R,fn,'Directional');
% replace nans with 0s
flattenedSTD_Win_m1_L(isnan(flattenedSTD_Win_m1_L))=0;
fn='~/m80_WinSTD.png'
Vis_Surf_n_Vecfield_NormativeDir(flattenedSTD_Win_m1_L,nets_RH,flattenedVectors_L,flattenedVectors_R,fn,'STD');
% replace nans with 0s
flattenedSTD_Btw_m1_L(isnan(flattenedSTD_Btw_m1_L))=0;
fn='~/m80_BWSTD.png'
Vis_Surf_n_Vecfield_NormativeDir(flattenedSTD_Btw_m1_L,nets_RH,flattenedVectors_L,flattenedVectors_R,fn,'STD');
% m2
fn='~/m120_directions.png'
Vis_Surf_n_Vecfield_NormativeDir(nets_LH,nets_RH,flattenedVectors_L_m2,flattenedVectors_R,fn,'Directional');
% replace nans with 0s
flattenedSTD_Win_m2_L(isnan(flattenedSTD_Win_m2_L))=0;
fn='~/m120_WinSTD.png'
Vis_Surf_n_Vecfield_NormativeDir(flattenedSTD_Win_m2_L,nets_RH,flattenedVectors_L,flattenedVectors_R,fn,'STD');
% replace nans with 0s
flattenedSTD_Btw_m2_L(isnan(flattenedSTD_Btw_m2_L))=0;
fn='~/m120_BWSTD.png'
Vis_Surf_n_Vecfield_NormativeDir(flattenedSTD_Btw_m2_L,nets_RH,flattenedVectors_L,flattenedVectors_R,fn,'STD');
% done
disp('donezo')
~                                                                                                              
~                                                                                                              
