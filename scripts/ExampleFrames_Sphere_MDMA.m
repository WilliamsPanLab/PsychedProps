function ExampleFrames_Sphere_MDMA(subj,sesh,task)
% add libraries
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% filepath to pull from
outFP=['/scratch/users/apines/data/mdma/' subj '/' sesh];
AngDistFPL=[outFP '/' subj '_' sesh '_' task '_k1_Prop_TS_dmn_L.csv'];
AngDistFPR=[outFP '/' subj '_' sesh '_' task '_k1_Prop_TS_dmn_R.csv'];
AngDistL=dlmread(AngDistFPL);
AngDistR=dlmread(AngDistFPR);
% combine them (Stack)
Props = [AngDistL; AngDistR];
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
	% for each frame
	for fr=propRange
                signalOfInt_L=TRs_l(:,fr);
                signalOfInt_R=TRs_r(:,fr);
		% init vec struct
		flattenedVectors_L=zeros(5120,3);
                flattenedVectors_R=zeros(5120,3);
		% pull vecs from cell struct
		flattenedVectors_L=vfs_L{fr};
		flattenedVectors_R=vfs_R{fr};
	        fn=['/scratch/users/apines/PropInstances/' subj '/' sesh '/' num2str(PropMean) '_subj_' sesh '_Prop' num2str(prop) '_fr' num2str(fr) '_sphere_directions.png'];
	         Vis_Surf_n_Vecfield_Sphere(signalOfInt_L,signalOfInt_R,flattenedVectors_L,flattenedVectors_R,fn,'Directional');
	         fn=['/scratch/users/apines/PropInstances/' subj '/' sesh '/' num2str(PropMean) '_subj_' sesh '_Prop' num2str(prop) '_fr' num2str(fr) '_sphere_BOLD.png'];
	         Vis_Surf_n_Vecfield_Sphere(signalOfInt_L,signalOfInt_R,flattenedVectors_L,flattenedVectors_R,fn,'BOLD');
	end
end

