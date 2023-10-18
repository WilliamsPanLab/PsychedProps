% calculate reference streams

% addpaths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))

% Dorsal Stream: Calc Fissure and S1
% Ventral Stream: Calc Fissure and Temporal Pole
% Insular: Purely into insula
% Medial Posterior: Calc Fissure into PCC
% Medial Anterior: M1 to mpfc

% load in fs4 and labels
%%% Load in surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
surfL = [SubjectsFolder '/lh.pial'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
% faces_L
F_L=faces_l;
% vertices V
V_L=vx_l;

%%%%%%% LOOK HERE
% load in spherical for vector norming
% change vx_l(CalcFisLocs,:) to centroid of CalcFis! and so on and so forth!

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
		CDist=EucCoords-vx_l(CalcFisLocs,:);
		% convert 3 coordinate difference to euclidean distance
		CEuclideanDists = sqrt(sum(CDist.^2, 2));	
		% find nearest calcarine vertex
		[Cminimum,minInd]=min(CEuclideanDists);
		% get vector pointing directly away from it (note away!)
		CAwayVector = VecNormalize(CDist(minInd, :));
		% get distance from S1
		SDist=EucCoords-vx_l(S1Locs,:);
		% convert 3 coordinates to difference in euclidean distance
		SEuclideanDists = sqrt(sum(SDist.^2, 2));
		% find nearest S1 vertex
		[Sminimum,minInd]=min(SEuclideanDists);
		% get vector pointing directly to it (inverse of CAway vector)
		STowardsVector = VecNormalize(-SDist(minInd,:));
		% weight them both of them by both distances for aggregate vector
		CS1Ratio=Cminimum/Sminimum;
		S1CRatio=Sminimum/Cminimum;
		AggVec=((CAwayVector.*CS1Ratio)+(STowardsVector.*S1CRatio))*2;
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
				CDist=EucCoords-vx_l(CalcFisLocs,:);
				% convert 3 coordinate difference to euclidean distance
				CEuclideanDists = sqrt(sum(CDist.^2, 2));
				% find nearest calcarine vertex
				[Cminimum,minInd]=min(CEuclideanDists);
				% get vector pointing directly away from it (note away!)
				CAwayVector = VecNormalize(CDist(minInd, :));
				% get distance from Temporal Pole
				TDist=EucCoords-vx_l(TempPoleLocs,:);
				% convert 3 coordinates to difference in euclidean distance
				TEuclideanDists = sqrt(sum(TDist.^2, 2));
				% find nearest Temporal Pole vertex
				[Tminimum,minInd]=min(TEuclideanDists);
				% get vector pointing directly to it (inverse of CAway vector)
				TTowardsVector = VecNormalize(-TDist(minInd,:));
				% weight them both of them by both distances for aggregate vector
				CTempRatio=Cminimum/Tminimum;
				TempCRatio=Tminimum/Cminimum;
				AggVec=((CAwayVector.*CTempRatio)+(TTowardsVector.*TempCRatio))*2;
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
		IDist=EucCoords-vx_l(InsulaLocs,:);
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
		CDist=EucCoords-vx_l(CalcFisLocs,:);
		% convert 3 coordinate difference to euclidean distance
		CEuclideanDists = sqrt(sum(CDist.^2, 2));
		% find nearest calcarine vertex
		[Cminimum,minInd]=min(CEuclideanDists);
		% get vector pointing directly away from it (note away!)
		CAwayVector = VecNormalize(CDist(minInd, :));
		% get distance from Sub Parietal Sulcus
		SDist=EucCoords-vx_l(SubPariLocs,:);
		% convert 3 coordinates to difference in euclidean distance
		SEuclideanDists = sqrt(sum(SDist.^2, 2));
		% find nearest Sub Parietal Sulcus vertex
		[Sminimum,minInd]=min(SEuclideanDists);
		% get vector pointing directly to it (inverse of CAway vector)
		STowardsVector = VecNormalize(-SDist(minInd,:));
		% weight them both of them by both distances for aggregate vector
		CSubRatio=Cminimum/Sminimum;
		SubCRatio=Sminimum/Cminimum;
		AggVec=((CAwayVector.*CSubRatio)+(STowardsVector.*SubCRatio))*2;
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
		% get distance from M1
		MDist=EucCoords-vx_l(M1Locs,:);
		% convert 3 coordinate difference to euclidean distance
		MEuclideanDists = sqrt(sum(MDist.^2, 2));
		% find nearest M1 vertex
		[Mminimum,minInd]=min(MEuclideanDists);
		% get vector pointing directly away from it (note away!)
		MAwayVector = VecNormalize(MDist(minInd, :));
		% get distance from suborbital
		SDist=EucCoords-vx_l(SubOrbLocs,:);
		% convert 3 coordinates to difference in euclidean distance
		SEuclideanDists = sqrt(sum(SDist.^2, 2));
		% find nearest suborbital vertex
		[Sminimum,minInd]=min(SEuclideanDists);
		% get vector pointing directly to it (inverse of MAway vector)
		STowardsVector = VecNormalize(-SDist(minInd,:));
		% weight them both of them by both distances for aggregate vector
		MSubRatio=Mminimum/Sminimum;
		SubMRatio=Sminimum/Mminimum;
		AggVec=((MAwayVector.*MSubRatio)+(STowardsVector.*SubMRatio))*2;
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
% dorsal stream needs vertices anterior to S1 removed, lateral to S1 removed, and inferior to V1 removed
% ventral stream needs vertices superior to V1 removed
% Insular needs medial vertices removed
% medial anterior needs lateral vertices removed
% medial posterior needs lateral vertices removed
% all could benefit from a vector smooth
