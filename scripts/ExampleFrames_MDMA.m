function ExampleFrames_MDMA(subj,sesh,task) 
% add libraries
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% filepath to pull from
outFP=['/scratch/users/apines/data/mdma/' subj '/' sesh];
AngDistFP=[outFP '/' subj '_' sesh '_' task '_k1_AngDistMat.mat'];
AngDist=load(AngDistFP).AngDist;
% Load in Prop TS L
PropsL=AngDist.Left;
% Prop TS R
PropsR=AngDist.Right;
% combine them (Stack)
Props = [PropsL; PropsR];
% average at each timepoint
colMeans = mean(Props, 1);
% find columns with 3 highest average values (to be labeled prop 1, 2, 3 in output)
[~, idx] = sort(colMeans, 'descend'); % Sort column means in descending order
top3_cols = idx(1:3); % Get indices of the top 3 columns
% print maximum full-brain value for TD prop angle
disp('3 highest average values:');
disp(colMeans(top3_cols));
% pull in time series
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
% load in time series
fpL=[childfp '/' subj '_' sesh '_task-' task '_p2mm_masked_L.mgh'];
fpR=[childfp '/' subj '_' sesh '_task-' task '_p2mm_masked_R.mgh'];
dataL=MRIread(fpL).vol;
dataR=MRIread(fpR).vol;
% files to data
TRs_l=squeeze(dataL);
TRs_r=squeeze(dataR);
% pull in vector fields
ofFP=[childfp '/' subj '_' sesh '_' task '_OpFl.mat'];
OpFl=load(ofFP);
% for each prop
for prop=1:3
	% index signal
	% index vector fields
	% orthoganalize vectors to inflated surface
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
	% print it out: 4 views (lateral + medial *L/R)
end
