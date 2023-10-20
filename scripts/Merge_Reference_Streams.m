function Merge_Reference_Stream
% convert downsampled network .giis to .mats
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/lukaslang-ofd-614a2ffc50d6'))
% this participant's filepath
funcgiiFolder = ['/oak/stanford/groups/leanew1/users/apines/data/RobustInitialization/'];
netgiis_L=[funcgiiFolder 'group_L_AggNets_3k.func.gii'];
netgiis_R=[funcgiiFolder 'group_R_AggNets_3k.func.gii'];
% load in functional networks: Left
Lnets=gifti(netgiis_L);
% load in functional networks: Right
Rnets=gifti(netgiis_R);
% set to 4 networks
Lnets=Lnets.cdata(:,1:4);
Rnets=Rnets.cdata(:,1:4);
% select DMN, k=2
DMNL=Lnets(:,2);
% boolean out areas less than .3 DMN
DMNL(DMNL<0.3)=0;

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
% get gradient of DMN
DMNGrad = grad(F_L, V_L, DMNL);

% load in manual reference streams
ManRef=load('/scratch/users/apines/refStreams.mat');
ManRef=ManRef.refStreams;
ReferenceStreams_L=zeros(5120,3,6);

% DMN
ReferenceStreams_L(:,:,1)=DMNGrad;

% Initialize arrays to store the face-wise values for x, y, and z components
F_X = zeros(size(F_L, 1), size(ManRef, 3));
F_Y = zeros(size(F_L, 1), size(ManRef, 3));
F_Z = zeros(size(F_L, 1), size(ManRef, 3));

% Convert vertex-wise values to face-wise values for each reference stream
for referenceStream = 1:5
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
for referenceStream = 1:5
    % Assign the converted values to the referenceStreams_L array
    % +1 because slot 1 is dmn
    ReferenceStreams_L(:, 1, referenceStream+1) = F_X(:, referenceStream);
    ReferenceStreams_L(:, 2, referenceStream+1) = F_Y(:, referenceStream);
    ReferenceStreams_L(:, 3, referenceStream+1) = F_Z(:, referenceStream);
end

% save out
save([funcgiiFolder 'group_ReferenceStreams_L.mat'],'ReferenceStreams_L')
