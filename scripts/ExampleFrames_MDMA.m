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
VFs=OpFl.us;
vfs_L=VFs.vf_left;
vfs_R=VFs.vf_right;

% get length of VF time series (should be < OG signal time series)
lengthOF=length(VFs.vf_left);
% get number of faces
lengthFaces=length(VFs.vf_left{1});

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

% create colormap
% roy-big-bl palette imitation, inferno is just template
roybigbl_cm=inferno(16);
roybigbl_cm(1,:)=[255, 255, 0 ];
roybigbl_cm(2,:)=[255, 200, 0];
roybigbl_cm(3,:)=[255, 120, 0];
roybigbl_cm(4,:)=[255, 0, 0 ];
roybigbl_cm(5,:)=[200, 0, 0 ];
roybigbl_cm(6,:)=[150, 0, 0 ];
roybigbl_cm(7,:)=[100, 0, 0 ];
roybigbl_cm(8,:)=[60, 0, 0 ];
roybigbl_cm(9,:)=[0, 0, 80 ];
roybigbl_cm(10,:)=[0, 0, 170];
roybigbl_cm(11,:)=[75, 0, 125];
roybigbl_cm(12,:)=[125, 0, 160];
roybigbl_cm(13,:)=[75, 125, 0];
roybigbl_cm(14,:)=[0, 200, 0];
roybigbl_cm(15,:)=[0, 255, 0];
roybigbl_cm(16,:)=[0, 255, 255];
% scale to 1
roybigbl_cm=roybigbl_cm.*(1/255);
% interpolate color gradient
interpsteps=[0 0.0666 0.1333 0.2 0.2666 0.3333 0.4 0.4666 0.5333 0.6 0.6666 0.7333 0.8 0.86666 0.9333 1];
roybigbl_cm=interp1(interpsteps,roybigbl_cm,linspace(0,1,255));
% yellow as high
roybigbl_cm=flipud(roybigbl_cm);
% reduce just a little bit on the close-to-white coloring
roybigbl_cm=roybigbl_cm(15:240,:);

% for each prop
for prop=1:3
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
		flattenedVectors_L=zeros(lengthFaces,3);
		flattenedVectors_R=zeros(lengthFaces,3);	
		% retrieve original vectors
                OGVecs_L=vfs_L{fr};
                OGVecs_R=vfs_R{fr};
		% for each point in space
		for F=1:lengthFaces
			% retrieve original vector
			OGVec_L=OGVecs_L(F,:);
			OGVec_R=OGVecs_R(F,:);
			% find each neighboring face: faces that share two vertices (should be 3)
			currFace_L=faces_l(F,:);
			currFace_R=faces_r(F,:);
			% Initialize neighboring faces storage
            		neighboringFaces_L = [];
            		neighboringFaces_R = [];
			% Find neighboring faces for the left hemisphere
            		for f = 1:size(faces_l, 1)
                		if f ~= F
                		    commonVertices = intersect(currFace_L, faces_l(f, :));
                		    if length(commonVertices) == 2
                		        neighboringFaces_L = [neighboringFaces_L; f];
		   		    end
               		 	end
            		end
			% right hemi
            		for f = 1:size(faces_r, 1)
                                if f ~= F
                                    commonVertices = intersect(currFace_R, faces_r(f, :));
                                    if length(commonVertices) == 2
                                        neighboringFaces_R = [neighboringFaces_R; f];   
                                    end   
                                end
                        end
			% initialize projection vector (directionality and magnitude w/r/t neighboring faces)
			projection_L=zeros(3,3);	
			faceCounter=1;
			% construe vector in terms of how much it points to each neighboring face
            		for neighborFace = neighboringFaces_L'
                		% Calculate normal vector for the neighboring face
                		normalVectors = cross(vx_l(faces_l(neighborFace, 2), :) - vx_l(faces_l(neighborFace, 1), :), ...
                                      vx_l(faces_l(neighborFace, 3), :) - vx_l(faces_l(neighborFace, 1), :));
                		meanNormalVector = mean(normalVectors, 1);
                		meanNormalVector = VecNormalize(meanNormalVector);
                
                		% Project the original vector onto the normal vector
                		projection_L(faceCounter,:) = dot(OGVec_L, meanNormalVector) * meanNormalVector;
				faceCounter=faceCounter+1;
            		end
			% right hemi
			projection_R=zeros(3,3);
                        faceCounter=1;
                        % construe vector in terms of how much it points to each neighboring face
                        for neighborFace = neighboringFaces_R'
                                % Calculate normal vector for the neighboring face
                                normalVectors = cross(vx_r(faces_r(neighborFace, 2), :) - vx_r(faces_r(neighborFace, 1), :), ...
                                      vx_r(faces_r(neighborFace, 3), :) - vx_r(faces_r(neighborFace, 1), :));
                                meanNormalVector = mean(normalVectors, 1);
                                meanNormalVector = VecNormalize(meanNormalVector);

                                % Project the original vector onto the normal vector
                                projection_R(faceCounter,:) = dot(OGVec_R, meanNormalVector) * meanNormalVector;
                                faceCounter=faceCounter+1;
                        end
			
			% reconstruct on inflated surface in terms of how much it points to each neighboring face
			reconstructedVec_L = zeros(1, 3);
			for i = 1:length(neighboringFaces_L)
                		neighborFace = neighboringFaces_L(i);
                		newVertexCoords_L = vx_li(faces_li(neighborFace, :), :);
                		meanNewVertex_L = mean(newVertexCoords_L, 1); % Get the mean position of the new face vertices
                		reconstructedVec_L = reconstructedVec_L + projection_L(i, :) .* meanNewVertex_L;
            		end
			% right hemisphere
			reconstructedVec_R = zeros(1, 3);
                        for i = 1:length(neighboringFaces_R)
                                neighborFace = neighboringFaces_R(i);
                                newVertexCoords_R = vx_ri(faces_ri(neighborFace, :), :);
                                meanNewVertex_R = mean(newVertexCoords_R, 1); % Get the mean position of the new face vertices
                                reconstructedVec_R = reconstructedVec_R + projection_R(i, :) .* meanNewVertex_R;
                        end
			
			%	
			%%%% Visualization check: compare spherical to inflated surface projection
			%
			%

			% Remove orthogonal component to flatten to the inflated surface
            		% Left hemisphere
			involvedFace_L = F;
            		normalVectors_L = cross(vx_li(faces_li(involvedFace_L, 2), :) - vx_li(faces_li(involvedFace_L, 1), :), ...
            		                        vx_li(faces_li(involvedFace_L, 3), :) - vx_li(faces_li(involvedFace_L, 1), :));
            		meanNormalVector_L = mean(normalVectors_L, 1);
            		meanNormalVector_L = VecNormalize(meanNormalVector_L);
            
            		vecOrthogonal_L = dot(reconstructedVec_L, meanNormalVector_L) * meanNormalVector_L;
            		modVec_L = reconstructedVec_L - vecOrthogonal_L;
            		modVec_L = VecNormalize(modVec_L);

            		% Right hemisphere
			involvedFace_R = F;
            		normalVectors_R = cross(vx_ri(faces_ri(involvedFace_R, 2), :) - vx_ri(faces_ri(involvedFace_R, 1), :), ...
                		                    vx_ri(faces_ri(involvedFace_R, 3), :) - vx_ri(faces_ri(involvedFace_R, 1), :));
            		meanNormalVector_R = mean(normalVectors_R, 1);
            		meanNormalVector_R = VecNormalize(meanNormalVector_R);
            
            		vecOrthogonal_R = dot(reconstructedVec_R, meanNormalVector_R) * meanNormalVector_R;
            		modVec_R = reconstructedVec_R - vecOrthogonal_R;
            		modVec_R = VecNormalize(modVec_R);

            		% Store the flattened vectors for later use or visualization
            		flattenedVectors_L(F, :) = modVec_L;
            		flattenedVectors_R(F, :) = modVec_R;			

			% remove orthogonal component of vectors to flatten to inflated surface
			% find the faces involved in this vertex 
			%[InvolvedFaces,~]=find(faces_l==vertInd);
			%normalVectors = cross(vx_l(faces_l(InvolvedFaces, 2), :) - vx_l(faces_l(InvolvedFaces, 1), :), vx_l(faces_l(InvolvedFaces, 3), :) - vx_l(faces_l(InvolvedFaces, 1), :));
			%meanNormalVector = mean(normalVectors, 1);
			% normalize normal vector
		        %meanNormalVector=VecNormalize(meanNormalVector);
			% get dot product of orthogonal vector and original vector
			%OGvecOrthogonal = dot(AggVec, meanNormalVector) * meanNormalVector;
			%modVec = AggVec - OGvecOrthogonal;
			% convert to unit vector
			%refStreams(vertInd,:,1)=VecNormalize(modVec);

			%
			%%%% Visualization check: compare surface-flattened to inflated surface projection
			%
			%
		end
		% end for each point in space
		% plot it!
		figure('units','pixels','position',[0 0 2000 2000])
		trisurf(faces_li, vx_li(:, 1), vx_li(:, 2), vx_li(:, 3),signalOfInt_L , 'EdgeColor','none')
		hold on;
		quiver3d(L_inc(1:1345,1),L_inc(1:1345,2),L_inc(1:1345,3),flattenedVectors_L(1:1345,1),flattenedVectors_L(1:1345,2),flattenedVectors_L(1:1345,3),[0.2 0.2 0.8])
		view([90 0]);
		colormap(roybigbl_cm);
		daspect([1 1 1]);
		axis tight;
		axis vis3d off;
		lighting none;
		shading flat;
		camlight;
		colorbar
		fn=['~/boldvec_' subj '_' sesh '_Prop' num2str(prop) '_fr' num2str(fr) '.png'];
        	print(fn,'-dpng')
		% print it out: 4 views (lateral + medial *L/R)
	% end for this frame
	end
% end for this prop
end
disp('done plotting prop instances')
