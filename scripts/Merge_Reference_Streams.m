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

% ensure all reference vectors are parallel to surface of sphere
% for each stream
for referenceStream = 1:12
        % for each face
        for f=1:length(F_L)
                % retrieve original vector
                OGvec=ReferenceStreams_L(f,:,referenceStream);
                % if reference stream is >0
                if OGvec(1)>0
                        % find the three faces involved in this vertex
                        InvolvedFace=F_L(f);
                        % get normal vectors of each involved face
                        normalVectors = cross(vertices(faces(InvolvedFace, 2), :) - vertices(faces(InvolvedFace, 1), :), vertices(faces(InvolvedFace, 3), :) - vertices(faces(InvolvedFace, 1), :));
                        % find vector as close as possible to orthogonal from these three faces
                        meanNormalVector = mean(normalVectors, 1);
                        % normalize normal vector
                        meanNormalVector=VecNormalize(meanNormalVector);
                        % get dot product of orthogonal vector and original vector
                        OGvecOrthogonal = dot(OGvec, meanNormalVector) * meanNormalVector;
                        % subtract orthogonal component of original vector from original vector
                        modVec = OGvec - OGvecOrthogonal;;
                        % add modified vector to initialized matrix
                        ReferenceStreams_L(v,:,referemceStream)=VecNormalize(modVec);
                % if its a zero vector
                else
                end
        end
end
% test visualization
Vis_Surf_n

%%%%%%% Now, load in cortical landmarks to anchor within-stream directionality
% read table
[v,label,ct]=read_annotation('/share/software/user/open/freesurfer/7.4.1/subjects/fsaverage4/label/lh.aparc.a2009s.annot');

% extract ROIs
CalcFisInd=ct.table(46,5);
CalcFisLocs=find(label==CalcFisInd);
M1Ind=ct.table(30,5);
M1Locs=find(label==M1Ind);
% need PCC S_subparietal 73
SubParSulcInd=ct.table(73,5);
SubParLocs=find(label==SubParSulcInd);
% need Insula G_Ins_lg_and_S_cent_ins 18
InsInd=ct.table(18,5);
InsLocs=find(label==InsInd);
% need PEF S_precentral-inf-part 70
PEFInd=ct.table(70,5);
PEFLocs=find(label==PEFInd);
% need A1 Transverse temporal sulcus 76
A1Ind=ct.table(76,5);
A1Locs=find(label==A1Ind);
% get centroid locations of each ROI
CalcFis_cent=mean(vx_l(CalcFisLocs,:));
M1_cent=mean(vx_l(M1Locs,:));
SubPar_cent=mean(vx_l(SubParLocs,:));
Ins_cent=mean(vx_l(InsLocs,:));
PEF_cent=mean(vx_l(PEFLocs,:));
A1_cent=mean(vx_l(A1Locs,:));

% will need some locations in real space (pial) for anteriormost point establishment: load in surf
surfL_p = ['/oak/stanford/groups/leanew1/users/apines/surf/lh.pial'];
% surface topography
[vx_lp, faces_lp] = read_surf(surfL);

% extract anteriormost point of OFC stream
AntInd=find(ReferenceStreams_L(:,1,4));
% get x y z locations of each face

% get maximal (x?y?z?) location
[AnteriorMostInd AntMost]=max(ReferenceStreams_L(AntInd,1,4)

% extract anteriormost point of cingulate stream
% get x y z locations of each face
% get maximal (x?y?z?) location

% make an array of reference points, xyz coordinates for each reference stream

%%% set reference point for each stream
% 1 = null
% 2 = V1
% 3 = V1
% 4 = anteriormost point of ofc stream
% 5 M1
% 6 PCC
% 7 = Premotor Eye field
% 8 = insula
% 9 = insula
% 10 = A1
% 11 = anteriormost point of cingulate stream
% 12 = V1

%%% align vectors in each stream
% for each stream
for referenceStream = 1:12
% for each vector
        for f=1:length(F_L)
% get opposite-pointing vector
% see if vector or it's inverse points more towards hierarchical ascent standard

% save out
save([funcgiiFolder 'group_ReferenceStreams_L.mat'],'ReferenceStreams_L')

