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
CalcFisLocs=v(CalcFisInd,:);
M1Ind=ct.table(30,5);
M1Locs=v(M1Ind,:);
S1Ind=ct.table(29,5);
S1Locs=v(S1Ind,:);
InsulaInd=ct.table(50,5);
InsulaLocs=v(InsulaInd,:);
TempPoleInd=ct.table(45,5);
TempPoleLocs=v(TempPoleInd,:);

% get lobe indices
[vl,labell,ctl]=read_annotation('/oak/stanford/groups/leanew1/users/apines/fsaverage4/label/lh.1.annot');
OccInd=ctl.table(5,5);
ParInd=ctl.table(7,5);
FrontInd=ctl.table(2,5);
InsLobeInd=ctl.table(8,5);
TempLobeInd=ctl.table(6,5);

% initialize a vector field for each stream (x5)
refStreams=zeros(2562,3,5);

% Dorsal Stream: Calc Fissure and S1
% mask by union of lobes: occipital and parietal 
OccPar = (labell == OccInd) | (labell == ParInd);
DStreamInds=find(OccPar);
% for each involved vertex, get distance from Calcarine and S1
for v=1:length(DStreamInds)
	% get real vertex index
	vertInd=DStreamInds(v);
	% get coordinates in euclidean space
	EucCoords=vx_l(vertInd,:);
	% get distance from calcarine
	Dist=EucCoords-CalcFisLocs
	% find nearest calcarine vertex
	% get vector pointing directly away from it (note away!)
	
	% get distance from S1
	% get vector pointing directly to it
	% weight the both of them by both distances for aggregate vector
	AggVec=
	%%% flatten to spherical surface (holllla divya ty)
	% get normal vector from nearby vertices
	% normvec step
	% 	
end



% Ventral Stream: Calc Fissure and Temporal Pole
% mask by union of lobes: temporal and occipital
OccTemp = (labell == TempLobeInd) | (labell == OccInd);

% Insular: Purely into insula
% mask by union of lobes: insula, temporal, frontal
InsFroTemp = (labell == InsLobeInd) | (labell == TempLobeInd) | (labell == FrontInd);

% Medial Posterior: Calc Fissure into PCC
% use same OccPar!

% Medial Anterior: M1 to mpfc
% mask by lobe: frontal
Front = (labell == FrontInd);




% will have to do something to convert verex-wise values to faces. Interpolate/average all vertices involved in the center of a face?



% per-face vector = sum of two components
% component 1 is directly away from Calc fissure * weighting (1/distance from calc fissure)
% component 2 is directly into S1 * weighting (1/distance from S1)
