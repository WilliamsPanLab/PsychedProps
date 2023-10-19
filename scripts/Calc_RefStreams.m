% calculate reference streams

% addpaths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))

% Dorsal Stream: Calc Fissure and S1
% Ventral Stream: Calc Fissure and Temporal Pole
% Insular: Purely into insula
% Medial Posterior: Calc Fissure into sub parietal sulcus
% Medial Anterior: M1 to suborbital sulcus

% load in fs4 and labels
%%% Load in surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfL = [SubjectsFolder '/lh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
% faces_L
F_L=faces_l;
% vertices V
V_L=vx_l;

% read table
[v,label,ct]=read_annotation('/share/software/user/open/freesurfer/7.4.1/subjects/fsaverage4/label/lh.aparc.a2009s.annot');

% extract ROIs
CalcFisInd=ct.table(46,5);
CalcFisLocs=find(label==CalcFisInd);
M1Ind=ct.table(30,5);
M1Locs=find(label==M1Ind);
S1Ind=ct.table(29,5);
S1Locs=find(label==S1Ind);
InsulaInd=ct.table(50,5);
InsulaLocs=find(label==InsulaInd);
TempPoleInd=ct.table(45,5);
TempPoleLocs=find(label==TempPoleInd);
SubPariInd=ct.table(73,5);
SubPariLocs=find(label==SubPariInd);
SubOrbInd=ct.table(72,5);
SubOrbLocs=find(label==SubOrbInd);
% get centroid locations of each ROI
CalcFis_cent=mean(vx_l(CalcFisLocs,:));
M1_cent=mean(vx_l(M1Locs,:));
S1_cent=mean(vx_l(S1Locs,:));
Insula_cent=mean(vx_l(InsulaLocs,:));
TempPole_cent=mean(vx_l(TempPoleLocs,:));
SubPari_cent=mean(vx_l(SubPariLocs,:));
SubOrb_cent=mean(vx_l(SubOrbLocs,:));

% get lobe indices
[vl,labell,ctl]=read_annotation('/oak/stanford/groups/leanew1/users/apines/fsaverage4/label/lh.1.annot');
OccInd=ctl.table(5,5);
ParInd=ctl.table(7,5);
FrontInd=ctl.table(2,5);
InsLobeInd=ctl.table(8,5);
TempLobeInd=ctl.table(6,5);

% initialize a vector field for each stream (x5)
refStreams=zeros(2562,3,5);

%%%%%%%%%%% Dorsal Stream: Calc Fissure and S1
% mask by union of lobes: occipital and parietal 
OccPar = (labell == OccInd) | (labell == ParInd);
DStreamInds=find(OccPar);
% for each involved vertex, get distance from Calcarine and S1
for v=1:length(DStreamInds)
	% get real vertex index
	vertInd=DStreamInds(v);
	% if it's not a vertex already in calc fis. or S1
	if label(vertInd)~=CalcFisInd && label(vertInd) ~= S1Ind
		% get coordinates in euclidean space
		EucCoords=vx_l(vertInd,:);
		% get distance from calcarine
		CDist=EucCoords-CalcFis_cent;
		% convert 3 coordinate difference to euclidean distance
		CEuclideanDists = sqrt(sum(CDist.^2, 2));	
		% find nearest calcarine vertex
		[Cminimum,minInd]=min(CEuclideanDists);
		% get vector pointing directly away from CENTROID %%% (note away!)
		CAwayVector = VecNormalize(CDist(minInd, :));
		% get distance from S1
		SDist=EucCoords-S1_cent;
		% convert 3 coordinates to difference in euclidean distance
		SEuclideanDists = sqrt(sum(SDist.^2, 2));
		% find nearest S1 vertex
		[Sminimum,minInd]=min(SEuclideanDists);
		% get vector pointing directly to it (inverse of CAway vector)
		STowardsVector = VecNormalize(-SDist(minInd,:));
		% weight them both of them by both distances for aggregate vector
		CS1Ratio=Cminimum/Sminimum;
		S1CRatio=Sminimum/Cminimum;
		AggVec=((CAwayVector.*S1CRatio)+(STowardsVector.*CS1Ratio));
		% flatten to surface
		% find the faces involved in this vertex 
		[InvolvedFaces,~]=find(faces_l==vertInd);
		normalVectors = cross(vx_l(faces_l(InvolvedFaces, 2), :) - vx_l(faces_l(InvolvedFaces, 1), :), vx_l(faces_l(InvolvedFaces, 3), :) - vx_l(faces_l(InvolvedFaces, 1), :));
		meanNormalVector = mean(normalVectors, 1);
		% normalize normal vector
	        meanNormalVector=VecNormalize(meanNormalVector);
		% get dot product of orthogonal vector and original vector
		OGvecOrthogonal = dot(AggVec, meanNormalVector) * meanNormalVector;
		modVec = AggVec - OGvecOrthogonal;
		% convert to unit vector
		refStreams(vertInd,:,1)=VecNormalize(modVec);
	else
	end
end

%%%%%%%%% Ventral Stream: Calc Fissure and Temporal Pole
% mask by union of lobes: temporal and occipital
OccTemp = (labell == TempLobeInd) | (labell == OccInd);
VStreamInds=find(OccTemp);
% for each involved vertex, get distance from Calcarine and Temporal Pole
for v=1:length(VStreamInds)
	% get real vertex index
	vertInd=VStreamInds(v);
	% if it's not a vertex already in calc fis. or S1
	if label(vertInd)~=CalcFisInd && label(vertInd) ~=TempPoleInd
			% get coordinates in euclidean space
			EucCoords=vx_l(vertInd,:);
			% get distance from calcarine
			CDist=EucCoords-CalcFis_cent;
			% convert 3 coordinate difference to euclidean distance
			CEuclideanDists = sqrt(sum(CDist.^2, 2));
			% find nearest calcarine vertex
			[Cminimum,minInd]=min(CEuclideanDists);
			% get vector pointing directly away from it (note away!)
			CAwayVector = VecNormalize(CDist(minInd, :));
			% get distance from Temporal Pole
			TDist=EucCoords-TempPole_cent;
			% convert 3 coordinates to difference in euclidean distance
			TEuclideanDists = sqrt(sum(TDist.^2, 2));
			% find nearest Temporal Pole vertex
			[Tminimum,minInd]=min(TEuclideanDists);
			% get vector pointing directly to it (inverse of CAway vector)
			TTowardsVector = VecNormalize(-TDist(minInd,:));
			% weight them both of them by both distances for aggregate vector
			CTempRatio=Cminimum/Tminimum;
			TempCRatio=Tminimum/Cminimum;
			AggVec=((CAwayVector.*TempCRatio)+(TTowardsVector.*CTempRatio));
			% flatten to surface
			% find the faces involved in this vertex
			[InvolvedFaces,~]=find(faces_l==vertInd);
			normalVectors = cross(vx_l(faces_l(InvolvedFaces, 2), :) - vx_l(faces_l(InvolvedFaces, 1), :), vx_l(faces_l(InvolvedFaces, 3), :) - vx_l(faces_l(InvolvedFaces, 1), :));
			meanNormalVector = mean(normalVectors, 1);
			% normalize normal vector
			meanNormalVector=VecNormalize(meanNormalVector);
			% get dot product of orthogonal vector and original vector
			OGvecOrthogonal = dot(AggVec, meanNormalVector) * meanNormalVector;
			modVec = AggVec - OGvecOrthogonal;
			% convert to unit vector
			refStreams(vertInd,:,2)=VecNormalize(modVec);
	else
	end
end

%%%%% Insular Stream: Purely into insula
% mask by union of lobes: insula, temporal, frontal
InsFroTemp = (labell == InsLobeInd) | (labell == TempLobeInd) | (labell == FrontInd);
InsInds=find(InsFroTemp);
% for each involved vertex, get distance from Insula
for v=1:length(InsInds)
	% get real vertex index
	vertInd=InsInds(v);
	% if it's not a vertex already in insula ROI
	if label(vertInd)~=InsulaInd
		% get coordinates in euclidean space
		EucCoords=vx_l(vertInd,:);
		% get distance from insula
		IDist=EucCoords-Insula_cent;
		% convert 3 coordinate difference to euclidean distance
		IEuclideanDists = sqrt(sum(IDist.^2, 2));
		% find nearest insula vertex
		[Iminimum,minInd]=min(IEuclideanDists);
		% get vector pointing directly to it
		ITowardsVector = VecNormalize(-IDist(minInd, :));
		vecsL(vertInd,:)=ITowardsVector;
		% flatten to surface
		% find the faces involved in this vertex
		[InvolvedFaces,~]=find(faces_l==vertInd);
		normalVectors = cross(vx_l(faces_l(InvolvedFaces, 2), :) - vx_l(faces_l(InvolvedFaces, 1), :), vx_l(faces_l(InvolvedFaces, 3), :) - vx_l(faces_l(InvolvedFaces, 1), :));
		meanNormalVector = mean(normalVectors, 1);
		% normalize normal vector
		meanNormalVector=VecNormalize(meanNormalVector);
		% get dot product of orthogonal vector and original vector
		OGvecOrthogonal = dot(ITowardsVector, meanNormalVector) * meanNormalVector;
		modVec = ITowardsVector - OGvecOrthogonal;
		% convert to unit vector
		refStreams(vertInd,:,3)=VecNormalize(modVec);
	else
	end
end

%%% Medial Posterior: Calc Fissure into Sub Parietal Sulcus
% mask by lobe: occipital and parietal
OccPar = (labell == OccInd) | (labell == ParInd);
MPInds=find(OccPar);
% for each involved vertex, get distance from Calcarine and Sub Parietal Sulcus
for v=1:length(MPInds)
	% get real vertex index
	vertInd=MPInds(v);
	% if it's not a vertex already in calc fis. or SubPariSulc
	if label(vertInd)~=CalcFisInd && label(vertInd) ~=SubPariInd
		% get coordinates in euclidean space
		EucCoords=vx_l(vertInd,:);
		% get distance from calcarine
		CDist=EucCoords-CalcFis_cent;
		% convert 3 coordinate difference to euclidean distance
		CEuclideanDists = sqrt(sum(CDist.^2, 2));
		% find nearest calcarine vertex
		[Cminimum,minInd]=min(CEuclideanDists);
		% get vector pointing directly away from it (note away!)
		CAwayVector = VecNormalize(CDist(minInd, :));
		% get distance from Sub Parietal Sulcus
		SDist=EucCoords-SubPari_cent;
		% convert 3 coordinates to difference in euclidean distance
		SEuclideanDists = sqrt(sum(SDist.^2, 2));
		% find nearest Sub Parietal Sulcus vertex
		[Sminimum,minInd]=min(SEuclideanDists);
		% get vector pointing directly to it (inverse of CAway vector)
		STowardsVector = VecNormalize(-SDist(minInd,:));
		% weight them both of them by both distances for aggregate vector, more important that it is directly into subPari than directly out of calc (*3)
		CSubRatio=3*(Cminimum/Sminimum);
		SubCRatio=Sminimum/Cminimum;
		AggVec=((CAwayVector.*SubCRatio)+(STowardsVector.*CSubRatio));
		% flatten to surface
		% find the faces involved in this vertex
		[InvolvedFaces,~]=find(faces_l==vertInd);
		normalVectors = cross(vx_l(faces_l(InvolvedFaces, 2), :) - vx_l(faces_l(InvolvedFaces, 1), :), vx_l(faces_l(InvolvedFaces, 3), :) - vx_l(faces_l(InvolvedFaces, 1), :));
		meanNormalVector = mean(normalVectors, 1);
		% normalize normal vector
		meanNormalVector=VecNormalize(meanNormalVector);
		% get dot product of orthogonal vector and original vector
		OGvecOrthogonal = dot(AggVec, meanNormalVector) * meanNormalVector;
		modVec = AggVec - OGvecOrthogonal;
		% convert to unit vector
		refStreams(vertInd,:,4)=VecNormalize(modVec);
	else
	end
end

%%%%%%% Medial Anterior: M1 to suborbital
% mask by lobe: frontal
Front = (labell == FrontInd);
MAInds=find(Front);
% for each involved vertex, get distance from M1 and suborbital
for v=1:length(MAInds)
	% get real vertex index
	vertInd=MAInds(v);
	% if it's not a vertex already in M1 or suborbital
	if label(vertInd)~=M1Ind && label(vertInd) ~= SubOrbInd
		% get coordinates in euclidean space
		EucCoords=vx_l(vertInd,:);
		% get distance from M1 (centroid is too lateral, need to use closest vertex of M1)
		MDist=EucCoords-vx_l(M1Locs,:);
		% convert 3 coordinate difference to euclidean distance
		MEuclideanDists = sqrt(sum(MDist.^2, 2));
		% find nearest M1 vertex
		[Mminimum,minInd]=min(MEuclideanDists);
		% get vector pointing directly away from it (note away!)
		MAwayVector = VecNormalize(MDist(minInd, :));
		% get distance from suborbital
		SDist=EucCoords-SubOrb_cent;
		% convert 3 coordinates to difference in euclidean distance
		SEuclideanDists = sqrt(sum(SDist.^2, 2));
		% find nearest suborbital vertex
		[Sminimum,minInd]=min(SEuclideanDists);
		% get vector pointing directly to it (inverse of MAway vector)
		STowardsVector = VecNormalize(-SDist(minInd,:));
		% weight them both of them by both distances for aggregate vector
		MSubRatio=Mminimum/Sminimum;
		SubMRatio=Sminimum/Mminimum*2;
		AggVec=((MAwayVector.*SubMRatio)+(STowardsVector.*MSubRatio));
		% flatten to surface
		% find the faces involved in this vertex
		[InvolvedFaces,~]=find(faces_l==vertInd);
		normalVectors = cross(vx_l(faces_l(InvolvedFaces, 2), :) - vx_l(faces_l(InvolvedFaces, 1), :), vx_l(faces_l(InvolvedFaces, 3), :) - vx_l(faces_l(InvolvedFaces, 1), :));
		meanNormalVector = mean(normalVectors, 1);
		% normalize normal vector
		meanNormalVector=VecNormalize(meanNormalVector);
		% get dot product of orthogonal vector and original vector
		OGvecOrthogonal = dot(AggVec, meanNormalVector) * meanNormalVector;
		modVec = AggVec - OGvecOrthogonal;
		% convert to unit vector
		refStreams(vertInd,:,5)=VecNormalize(modVec);
	else
	end
end

% insert lobular corrections!
% real in pial for euclidean approximation
surfL = [SubjectsFolder '/lh.pial'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
% update centroids to pial surface
CalcFis_cent=mean(vx_l(CalcFisLocs,:));
M1_cent=mean(vx_l(M1Locs,:));
S1_cent=mean(vx_l(S1Locs,:));
Insula_cent=mean(vx_l(InsulaLocs,:));
TempPole_cent=mean(vx_l(TempPoleLocs,:));
SubPari_cent=mean(vx_l(SubPariLocs,:));
SubOrb_cent=mean(vx_l(SubOrbLocs,:));

% dorsal stream needs vertices anterior to S1 removed and inferior to V1 removed
% anterior to s1 removal purely by masking out central sulcus
CSInd=ct.table(47,5);
CSLocs=find(label==CSInd);
refStreams(CSLocs,:,1)=0;
% inferior to V1 removed, 3 for z dimension
DStreamInds_tooInf=DStreamInds(vx_l(DStreamInds,3)<(CalcFis_cent(3)+3));
refStreams(DStreamInds_tooInf,:,1)=0;
% lateral to S1 removed
%DStreamInds_tooLat=DStreamInds(vx_l(DStreamInds,1)<-55);
%refStreams(DStreamInds_tooLat,:,1)=0;
% remove supramarginal gyrus, planum temporale, and subcentral gyrus
SupraMargInd=ct.table(27,5);
SupraMargLocs=find(label==SupraMargInd);
PLTInd=ct.table(37,5);
PLTLocs=find(label==PLTInd);
SCGInd=ct.table(5,5);
SCGLocs=find(label==SCGInd);
% remove Lat_Fis-post and S_circular_insula_sup
LFPInd=ct.table(42,5);
LFPLocs=find(label==LFPInd);
SCIInd=ct.table(51,5);
SCILocs=find(label==SCIInd);
% mask out
refStreams(SupraMargLocs,:,1)=0;
refStreams(PLTLocs,:,1)=0;
refStreams(SCGLocs,:,1)=0;
refStreams(LFPLocs,:,1)=0;
refStreams(SCILocs,:,1)=0;

% ventral stream needs vertices superior to V1 removed
VStreamInds_tooSup=VStreamInds(vx_l(VStreamInds,3)>(CalcFis_cent(3)));
refStreams(VStreamInds_tooSup,:,2)=0;

% Insular needs medial vertices removed
InsInds_tooMed=InsInds(vx_l(InsInds,1)>-20);
refStreams(InsInds_tooMed,:,3)=0;

% medial posterior needs lateral vertices removed
MPInds_tooLat=MPInds(vx_l(MPInds,1)<-20);
refStreams(MPInds_tooLat,:,4)=0;
% and inferior to V1
MPInds_tooInf=MPInds(vx_l(MPInds,3)<(CalcFis_cent(3)+3));
refStreams(MPInds_tooInf,:,4)=0;
% and posterior to V1
MPInds_tooPost=MPInds(vx_l(MPInds,2)<(CalcFis_cent(2)));
refStreams(MPInds_tooPost,:,4)=0;
% and superior to sub parietal sulcus
MPInds_tooSup=MPInds(vx_l(MPInds,3)>(SubPari_cent(3)));
refStreams(MPInds_tooSup,:,4)=0;

% medial anterior needs lateral vertices removed
MAInds_tooLat=MAInds(vx_l(MAInds,1)<-20);
refStreams(MAInds_tooLat,:,5)=0;
% and vertices posterior to m1 removed
MAInds_tooPost=MAInds(vx_l(MAInds,2) < (M1_cent(2)));
refStreams(MAInds_tooPost,:,5)=0;

% make dummy surface values and visualize each stream
DsurfL=zeros(2562,1);
DsurfR=zeros(2562,1);
DVecsR=zeros(2562,3);

% vis dorsal stream
DsurfL_d=DsurfL;
DsurfL_d(CalcFisLocs)=1;
DsurfL_d(S1Locs)=2;
Vis_Surf_n_Vecfield(DsurfL_d,DsurfR,refStreams(:,:,1),refStreams(:,:,1),'~/Dorsal_Stream.png')
% vis ventral stream
DsurfL_v=DsurfL;
DsurfL_v(CalcFisLocs)=1;
DsurfL_v(TempPoleLocs)=2;
Vis_Surf_n_Vecfield(DsurfL_v,DsurfR,refStreams(:,:,2),refStreams(:,:,2),'~/Ventral_Stream.png')
% vis insular stream
DsurfL_i=DsurfL;
DsurfL_i(InsulaLocs)=1;
Vis_Surf_n_Vecfield(DsurfL_i,DsurfR,refStreams(:,:,3),refStreams(:,:,3),'~/Insular_Stream.png')
% vis medial posterior stream
DsurfL_mp=DsurfL;
DsurfL_mp(CalcFisLocs)=1;
DsurfL_mp(SubPariLocs)=2;
Vis_Surf_n_Vecfield(DsurfL_mp,DsurfR,refStreams(:,:,4),refStreams(:,:,4),'~/MedialPosterior_Stream.png')
% vis medial posterior stream
DsurfL_ma=DsurfL;
DsurfL_ma(M1Locs)=1;
DsurfL_ma(SubOrbLocs)=2;
Vis_Surf_n_Vecfield(DsurfL_ma,DsurfR,refStreams(:,:,5),refStreams(:,:,5),'~/MedialAnterior_Stream.png')

% all could benefit from a vector smooth
