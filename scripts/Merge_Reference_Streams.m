function Merge_Reference_Stream
% convert downsampled network .giis to .mats
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/lukaslang-ofd-614a2ffc50d6'))
% this participant's filepath
funcgiiFolder = ['/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/'];

% Left only for now

% get some freesurfer info
% Load in surface data
surfL = ['/oak/stanford/groups/leanew1/users/apines/surf/lh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
% faces_L
F_L=faces_l;
% vertices V
V_L=vx_l;

% load in manual reference streams
ManRef=load('~/vertexwise_RefVecs.mat');
ManRef=ManRef.ComDimVF;
ReferenceStreams_L=zeros(5120,3,12);

% Initialize arrays to store the face-wise values for x, y, and z components
F_X = zeros(size(F_L, 1), 12);
F_Y = zeros(size(F_L, 1), 12);
F_Z = zeros(size(F_L, 1), 12);

% Convert vertex-wise values to face-wise values for each reference stream
for referenceStream = 1:12
    for faceIdx = 1:size(F_L, 1)
        vertexIndices = F_L(faceIdx, :);  % Vertex indices for the current face
        % Compute the average x, y, and z components for the current face and reference stream
        avgX = mean(ManRef(vertexIndices, 1, referenceStream));
        avgY = mean(ManRef(vertexIndices, 2, referenceStream));
        avgZ = mean(ManRef(vertexIndices, 3, referenceStream));

        % Store the average values in the corresponding arrays
        F_X(faceIdx, referenceStream) = avgX;
        F_Y(faceIdx, referenceStream) = avgY;
        F_Z(faceIdx, referenceStream) = avgZ;
    end
end

% plug these values into referenceStreams_L
for referenceStream = 1:12
    % Assign the converted values to the referenceStreams_L array
    ReferenceStreams_L(:, 1, referenceStream) = F_X(:, referenceStream);
    ReferenceStreams_L(:, 2, referenceStream) = F_Y(:, referenceStream);
    ReferenceStreams_L(:, 3, referenceStream) = F_Z(:, referenceStream);
end

% get center of sphere for norming to surface
sphereCenter=mean(V_L);

% get mean directionality of each stream to anchor invidivual vectors
MeanVecs=zeros(12,3);
for referenceStream = 2:12
	RefVecs=ReferenceStreams_L(:,:,referenceStream);
	ExtantVecs=find(RefVecs(:,1));
	MeanVecs(referenceStream,:)=mean(RefVecs(ExtantVecs,:));
	% now that we have a normative direction for each stream, re-align each vector to be itself or its inverse, whichever is more aligned with normative direction
	for faceIdx = 1:size(F_L, 1)
		% if its a non-zero vector
		if sum(abs(ReferenceStreams_L(faceIdx,:,referenceStream)))>0
			% get face's vector
			fvec=ReferenceStreams_L(faceIdx,:,referenceStream);
			%%% parallelize to surface
			Curvertices = V_L(F_L(faceIdx, :), :);
			% Calculate the midpoint of the current face
			midpoint = mean(Curvertices);
			% vector to current position and sphere origin
			vectorToMidpoint = midpoint - sphereCenter;
			% Calculate the unit normal vector at the midpoint
			normalVector = vectorToMidpoint / norm(vectorToMidpoint);
			% Calculate the dot product between the vector and the normal vector
			dotProduct = dot(fvec, normalVector);
			% Calculate the orthogonal component of the vector
			orthogonalComponent = dotProduct * normalVector;
			% Calculate the vector parallel to the sphere's surface
			vectorParallelToSurface = fvec - orthogonalComponent;
			% get dot product of mean stream direction and current vector
%			dotProductSphereSurf=dot(vectorParallelToSurface,MeanVecs(referenceStream, :));
			% and inverse
%			ifvec=-1*vectorParallelToSurface;
			% see which is closer, sub it into OG
%			dotProductInverse = dot(ifvec, MeanVecs(referenceStream, :));
			% larger dot product indicates greater alignment
%			if dotProductSphereSurf < dotProductInverse
%				 ReferenceStreams_L(faceIdx, :, referenceStream) = ifvec;
%			else
				ReferenceStreams_L(faceIdx, :, referenceStream) = vectorParallelToSurface;
%			end
		else
		end
	end
	% make a face-wise vector just to fulfill script requirements
	dummyvec=zeros(5120,1);
	dummyvec(10)=1;
	% visualize
	%Vis_Surf_n_Vecfield_Faces(dummyvec,dummyvec,ReferenceStreams_L(:, :, referenceStream),ReferenceStreams_L(:, :, referenceStream),['~/MeanDir_aligned_' num2str(referenceStream) '.png']);
	%pause(10)
end
% saveout
save([funcgiiFolder 'goup_ReferenceStreams_L.mat'],'ReferenceStreams_L')
