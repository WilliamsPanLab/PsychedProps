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

% save out
save([funcgiiFolder 'group_ReferenceStreams_L.mat'],'ReferenceStreams_L')
