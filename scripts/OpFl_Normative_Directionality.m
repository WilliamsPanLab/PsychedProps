% read in subject dosage correspondence, has to be before addpath for some silly reason
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
% load each subject and get mean value within and across subjects at each vertex
% get subj list
subjPrefix=repmat('sub-MDMA0',17,1);
subjSuffix=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17"];
subjList=strcat(subjPrefix,subjSuffix')

% variable task name
task='rs1'

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
	Thetas_L=zeros(1,5120);
	Rhos_L=zeros(1,5120);
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
		end
		% get average angle over time
		bv_g_anglesx(s,F)=mean(x_Face);
		bv_g_anglesy(s,F)=mean(y_Face);
		bv_g_anglesz(s,F)=mean(z_Face);
		% might still use this for Standard deviation (circular)
		vs_L=cart2sphvec(double([mean(x_Face);mean(y_Face);mean(z_Face)]),azd_L(F),eld_L(F));
		[Thetas_L(F),Rhos_L(F)]=cart2pol(vs_L(1),vs_L(2));
	end
	% populate baseline dataframe leaving this out for now to keep aggregated angles as cartesian
	%[bv_g_anglesx(s,:),bv_g_anglesy(s,:)]=pol2cart(Thetas_L,Rhos_L);
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
	Thetas_L=zeros(1,5120);
	Rhos_L=zeros(1,5120);
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
		end
		% get average angle over time
                p_g_anglesx(s,F)=mean(x_Face);
                p_g_anglesy(s,F)=mean(y_Face);
                p_g_anglesz(s,F)=mean(z_Face);
		% convert to spherical coord system
	       	vs_L=cart2sphvec(double([mean(x_Face);mean(y_Face);mean(z_Face)]),azd_L(F),eld_L(F));
		% store in output vector (r is redundant across all vecs, only using az and el)
		[Thetas_L(F),Rhos_L(F)]=cart2pol(vs_L(1),vs_L(2));
	end
	% populate placebo dataframe - leaving out now in favor of cartesian storage
        %[p_g_anglesx(s,:),p_g_anglesy(s,:)]=pol2cart(Thetas_L,Rhos_L);

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
	Thetas_L=zeros(1,5120);
	Rhos_L=zeros(1,5120);
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
		end
		% get average angle over time
                m1_g_anglesx(s,F)=mean(x_Face);
                m1_g_anglesy(s,F)=mean(y_Face);
                m1_g_anglesz(s,F)=mean(z_Face);
		% convert to spherical coord system
	       	vs_L=cart2sphvec(double([mean(x_Face);mean(y_Face);mean(z_Face)]),azd_L(F),eld_L(F));
		% store in output vector (r is redundant across all vecs, only using az and el)
		[Thetas_L(F),Rhos_L(F)]=cart2pol(vs_L(1),vs_L(2));
	end
	% populate m1 dataframe
        % leaving out for now in favor of cartesian storage[m1_g_anglesx(s,:),m1_g_anglesy(s,:)]=pol2cart(Thetas_L,Rhos_L);
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
	Thetas_L=zeros(1,5120);
	Rhos_L=zeros(1,5120);
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
		end
		% get average angle over time
                m2_g_anglesx(s,F)=mean(x_Face);
                m2_g_anglesy(s,F)=mean(y_Face);
                m2_g_anglesz(s,F)=mean(z_Face);
		% convert to spherical coord system
	       	vs_L=cart2sphvec(double([mean(x_Face);mean(y_Face);mean(z_Face)]),azd_L(F),eld_L(F));
		% store in output vector (r is redundant across all vecs, only using az and el)
		[Thetas_L(F),Rhos_L(F)]=cart2pol(vs_L(1),vs_L(2));
	end
	% populate m2 dataframe commenting out for now in favor of cartesian storage
        %[m2_g_anglesx(s,:),m2_g_anglesy(s,:)]=pol2cart(Thetas_L,Rhos_L);
	else
	end
end

% remove blank rows (should be 4, 6, 10)
bv_g_anglesx = bv_g_anglesx(~all(bv_g_anglesx == 0, 2), :);
bv_g_anglesy = bv_g_anglesy(~all(bv_g_anglesy == 0, 2), :);
bv_g_anglesz = bv_g_anglesz(~all(bv_g_anglesz == 0, 2), :);
p_g_anglesx = p_g_anglesx(~all(p_g_anglesx == 0, 2), :);
p_g_anglesy = p_g_anglesy(~all(p_g_anglesy == 0, 2), :);
p_g_anglesz = p_g_anglesz(~all(p_g_anglesz == 0, 2), :);
m1_g_anglesx = m1_g_anglesx(~all(m1_g_anglesx == 0, 2), :);
m1_g_anglesy = m1_g_anglesy(~all(m1_g_anglesy == 0, 2), :);
m1_g_anglesz = m1_g_anglesz(~all(m1_g_anglesz == 0, 2), :);
m2_g_anglesx = m2_g_anglesx(~all(m2_g_anglesx == 0, 2), :);
m2_g_anglesy = m2_g_anglesy(~all(m2_g_anglesy == 0, 2), :);
m2_g_anglesz = m2_g_anglesz(~all(m2_g_anglesz == 0, 2), :);
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

% this would be a good spot to orthogonalize to the surface of inflated
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

% change OG vecs for each condition. First = placebo
OGVecs_L=[p_g_angles_meanx; p_g_angles_meany; p_g_angles_meanz]';
% get magnitude of each vector
magnitudes = sqrt(sum(OGVecs_L.^2, 2));
% for each vertex
for vertInd=1:2562
	% get faces intertwined with this vertex
	[InvolvedFaces_l,~]=find(faces_l==vertInd);
	% get mean vector directionality at these faces (x y z)
	OGVec_L=mean(OGVecs_L(InvolvedFaces_l,:));
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
% set right hemisphere vectors as null
flattenedVectors_R=zeros(2562,3);
% adjust flattenedVectors as function of original magnitude
flattenedVectors_L=flattenedVectors_L/magnitudes;
% read in DMN as background
networks=load(['/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/gro_Nets_fs4.mat']);
nets_LH=networks.nets.Lnets(:,1);
nets_RH=networks.nets.Rnets(:,1);
fn='~/testAggVecs.png'
Vis_Surf_n_Vecfield_NormativeDir(nets_LH,nets_RH,flattenedVectors_L,flattenedVectors_R,fn,'Directional');


% ADD IN CIRCULAR SD, CALC WITHIN SUBJECTS 17X5120, CAN THEN AVERAGE
% ADD IN CIRCULAR SD FOR BETWEEN SUBJECTS: 1X5120 (FROM 17X5120 NORMATIVE DIRECTIONALITY
% COLOR ARROWS BY SD: WILL NEED TO RUN SD VALUES NEXT TO INTERPOLATION TO VERTICES
