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
[~, idx] = sort(colMeans, 'ascend'); % Sort column means in descending order
top3_cols = idx(1:3); % Get indices of the top 3 columns
% print maximum full-brain value for TD prop angle
disp('3 lowest average values:');
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
VFs=OpFl.us;
vfs_L=VFs.vf_left;
vfs_R=VFs.vf_right;

% get length of VF time series (should be < OG signal time series)
lengthOF=length(VFs.vf_left);
% get number of faces
lengthFaces=length(VFs.vf_left{1});
% get length of vertices
lengthVertices=length(TRs_l);
%%% load in freesurfer surfaces
% Load in surface data - sphere
surfL = ['/oak/stanford/groups/leanew1/users/apines/surf/lh.sphere'];
surfR = ['/oak/stanford/groups/leanew1/users/apines/surf/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;

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
% make directories in scratch for saveouts
subjDirtoMake=['/scratch/users/apines/PropInstances/' subj];
system(['mkdir ', subjDirtoMake]);
sesDirtoMake=['/scratch/users/apines/PropInstances/' subj '/' sesh];
system(['mkdir ', sesDirtoMake]);
% for each prop
for prop=1:3
	% pull out average degrees from DMN
	PropMeans=colMeans(top3_cols);
	PropMean=PropMeans(prop);
	% pull out this propagation timepoint
	propInd=top3_cols(prop);
	% extrapolate to a range: 4 before and 4 after, but need to account for possibility of this range extended beyond time series
	rangeStart = max(propInd - 4, 1); % ensure the start is not less than 1
	rangeEnd = min(propInd + 4, lengthOF); % same for ending range, not > timeseries length
	propRange = rangeStart:rangeEnd;
	% index vector fields
	VFsOfInt_L=vfs_L{propRange};
	VFsOfInt_R=vfs_R{propRange};
	% for each frame
	for fr=propRange
		% index signal
		signalOfInt_L=TRs_l(:,fr);
		signalOfInt_R=TRs_r(:,fr);
		% initialzie flattened vectors array for this frame
		flattenedVectors_L=zeros(lengthVertices,3);
		flattenedVectors_R=zeros(lengthVertices,3);	
		% retrieve original vectors
                OGVecs_L=vfs_L{fr};
                OGVecs_R=vfs_R{fr};
		% for each point in space
		for vertInd=1:lengthVertices
			% get faces involved in this vertex
			[InvolvedFaces_l,~]=find(faces_l==vertInd);
			[InvolvedFaces_r,~]=find(faces_r==vertInd);
			% retrieve original vectors and average them across these faces
			OGVec_L=mean(OGVecs_L(InvolvedFaces_l,:));
			OGVec_R=mean(OGVecs_R(InvolvedFaces_r,:));
			% initialize projection vector (directionality and magnitude w/r/t neighboring faces)
			projection_L=zeros(5,3);
			projection_R=zeros(5,3);	
			faceCounter=1;
			% constriue vector in terms of how much it points to each neighboring face
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

			% Construe vector in terms of how much it points to each neighboring face for right hemisphere
            		faceCounter = 1;
            		for neighborFace = InvolvedFaces_r'
                		% Calculate normal vector for the neighboring face
                		normalVectors = cross(vx_r(faces_r(neighborFace, 2), :) - vx_r(faces_r(neighborFace, 1), :), ...
                        		              vx_r(faces_r(neighborFace, 3), :) - vx_r(faces_r(neighborFace, 1), :));
                		meanNormalVector = mean(normalVectors, 1);
                		meanNormalVector = VecNormalize(meanNormalVector);

                		% Project the original vector onto the normal vector
                		projection_R(faceCounter, :) = dot(OGVec_R, meanNormalVector) * meanNormalVector;
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
            		% right hemisphere
            		reconstructedVec_R = zeros(1, 3);
            		for i = 1:length(InvolvedFaces_r)
            		    neighborFace = InvolvedFaces_r(i);
            		    newVertexCoords_R = vx_ri(faces_ri(neighborFace, :));
            		    meanNewVertex_R = mean(newVertexCoords_R, 1); % Get the mean position of the new face vertices
            		    reconstructedVec_R = reconstructedVec_R + projection_R(i, :) .* meanNewVertex_R; % Element-wise multiplication
            		end
			
			% Normalize the reconstructed vectors
            		reconstructedVec_L = VecNormalize(reconstructedVec_L);
            		reconstructedVec_R = VecNormalize(reconstructedVec_R);	
			

			% Remove orthogonal component to flatten to the inflated surface
            		% Left hemisphere
            		involvedFaces_L = find(any(faces_li == vertInd, 2)); % Find faces involving this vertex
            		normalVectors_L = cross(vx_li(faces_li(involvedFaces_L, 2), :) - vx_li(faces_li(involvedFaces_L, 1), :), ...
                		                    vx_li(faces_li(involvedFaces_L, 3), :) - vx_li(faces_li(involvedFaces_L, 1), :));
            		meanNormalVector_L = mean(normalVectors_L, 1);
            		meanNormalVector_L = VecNormalize(meanNormalVector_L);
            
            		OGvecOrthogonal_L = dot(reconstructedVec_L, meanNormalVector_L) * meanNormalVector_L;
            		modVec_L = reconstructedVec_L - OGvecOrthogonal_L;
            		modVec_L = VecNormalize(modVec_L);

            		% Right hemisphere
            		involvedFaces_R = find(any(faces_ri == vertInd, 2)); % Find faces involving this vertex
            		normalVectors_R = cross(vx_ri(faces_ri(involvedFaces_R, 2), :) - vx_ri(faces_ri(involvedFaces_R, 1), :), ...
            		                        vx_ri(faces_ri(involvedFaces_R, 3), :) - vx_ri(faces_ri(involvedFaces_R, 1), :));
            		meanNormalVector_R = mean(normalVectors_R, 1);
            		meanNormalVector_R = VecNormalize(meanNormalVector_R);
            
            		OGvecOrthogonal_R = dot(reconstructedVec_R, meanNormalVector_R) * meanNormalVector_R;
            		modVec_R = reconstructedVec_R - OGvecOrthogonal_R;
            		modVec_R = VecNormalize(modVec_R);

            		% Store the flattened vectors for later use or visualization
            		flattenedVectors_L(vertInd, :) = modVec_L;
            		flattenedVectors_R(vertInd, :) = modVec_R;

			%	
			%%%% Visualization check: compare spherical to inflated surface projection
			%
			%

			%
			%%%% Visualization check: compare surface-flattened to inflated surface projection
			%
			%
		end
		% end for each point in space
		% plot it!
		fn=['/scratch/users/apines/PropInstances/' subj '/' sesh '/' num2str(PropMean) '_subj_' sesh '_Prop' num2str(prop) '_fr' num2str(fr) '_directions.png'];
		Vis_Surf_n_Vecfield(signalOfInt_L,signalOfInt_R,flattenedVectors_L,flattenedVectors_R,fn,'Directional');
		fn=['/scratch/users/apines/PropInstances/' subj '/' sesh '/' num2str(PropMean) '_subj_' sesh '_Prop' num2str(prop) '_fr' num2str(fr) '_BOLD.png'];
                Vis_Surf_n_Vecfield(signalOfInt_L,signalOfInt_R,flattenedVectors_L,flattenedVectors_R,fn,'BOLD');
	% end for this frame
	end
% end for this prop
end
disp('done plotting prop instances')
