function OpFlStreamlines_null_Left(subj,sesh,task)
% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
childfp=['/scratch/users/apines/'];
% load in optical flow output
data=load([childfp 'Simulated_OpFl_' subj '_' sesh '_' task '.mat']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uptake surface data
SubjectsFolder = '/oak/stanford/groups/leanew1/users/apines/surf';
% for surface data
surfL = [SubjectsFolder '/lh.sphere'];
surfR = [SubjectsFolder '/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% normalize verts to unit sphere
% left
numV=length(vx_l);
vx_l(numV+1:end, :) = VecNormalize(vx_l(numV+1:end, :));
% right
numV=length(vx_r);
vx_r(numV+1:end, :) = VecNormalize(vx_r(numV+1:end, :));

% calculate normative directionality streamlines for this subject and session
n = size(faces_l, 1);
T = TriRep(faces_l, vx_l);
P_L = T.incenters;
n = size(faces_r, 1);
T = TriRep(faces_r, vx_r);
P_R = T.incenters;
%%% L as left Hemi
L=data.us.vf_left;
%%% R as right
R=data.us.vf_right;
% use native freesurfer command for mw mask indices
surfML = [SubjectsFolder '/lh.Medial_wall.label'];
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = [SubjectsFolder '/rh.Medial_wall.label'];
mwIndVec_r = read_medial_wall_label(surfMR);
% make "isn't medial wall" indices
nonMW_L=setdiff(1:2562,mwIndVec_l);
nonMW_R=setdiff(1:2562,mwIndVec_r);
%%% get length of time serires
tsLength=length(L);

% initialize vector field to be plotted
VF_L=zeros(n,3,tsLength);
VF_R=zeros(n,3,tsLength);
% loop over to extract from unfortunate cell structure
for i=1:tsLength;
	VF_L(:,:,i)=L{i};
	VF_R(:,:,i)=R{i};
end
% correct spatial coordinates to be on scale of vector field
P_L=P_L./100;
P_R=P_R./100;
vx_l=vx_l./100;
vx_r=vx_r./100;
% initialize SC matrix equivalent
AdjMatrix_L=zeros(length(vx_l),length(vx_l));
AdjMatrix_R=zeros(length(vx_r),length(vx_r));


% this is annoying, but for extra clarity we'll be matching the simulated optical flow vector magnitude to the mean magnitude from real
% load in subj OpFl
childfp2=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
% load in optical flow output
data=load([childfp2 '/' subj '_' sesh '_' task '_OpFl.mat']);
Lreal=data.us.vf_left;
% get mean magnitude
VF_Lreal=zeros(n,3,tsLength);
for i=1:tsLength;
        VF_Lreal(:,:,i)=Lreal{i};
end
% get mean vector value
mVVreal=mean(mean(mean(VF_Lreal)));
% get ratio of mean vector values of real to sim
mVVr=mVVreal/(mean(mean(mean(VF_L))));
% correct null vectors
VF_L=VF_L.*mVVr;

% open up the pool
pool=parpool('local');

% get subject's motion mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load in continuous segment indices
childfp_subj=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
CSIfp=[childfp_subj '/' subj '_' sesh '_task-' task '_ValidSegments_Trunc.txt'];
CSI = importdata(CSIfp);
% assure that TR count is the same between time series and valid segments txt
SegNum=size(CSI);
SegNum=SegNum(1);
% number of trs in fmri ts
vf_size=size(VF_L);
mr_ts_trs=vf_size(3);
% trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
numTRsVS=CSI(SegNum,1)+CSI(SegNum,2)-1;
% - segnum because this is tr pairs in mr_ts_trs, not TRs
if numTRsVS ~= (mr_ts_trs + SegNum)
        disp('TRs from Valid Segments txt and mgh file do not match. Fix it.')
        return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for each vertex
parfor v=1:length(vx_l)
	% print v
	v
	% initialize row for this vertex to index into
	adjacency_row = zeros(1, length(vx_l));

	% for each continuous segment
        for seg=1:SegNum	
		 % just to print out count of current segment
                seg
                SegStart=CSI(seg,1);
                % minus 1 because it's between TR measures and discontinuity incurs missing between volume
                SegSpan=CSI(seg,2)-1;
                % get corresponding TRs from aggregate time series
                % note that we lose one "between frame" for each segment, as the bookend has no between calculation with the next book"begin"
                segTS_l=VF_L(:,:,(SegStart:((SegStart+SegSpan)-seg)));
		% loop over each TR-Pair: 1 fewer pair than number of TRs
                for TRP=1:(SegSpan)
			% plant a new seed
			CurrSeed=vx_l(v,:);
			% set endpoint of tracking: if t + 85 > length of time series, cant track beyond ts
                        % allowing them to run for up to one minute of scan time (current tr = .71)
                        endpoint=85;
                        if (endpoint+TRP)>SegSpan
                                endpoint=SegSpan-TRP;
                        else
                        end
			% for a max of one minute
			for t2=1:endpoint
				% find the three nearest faces
				[nearFaces,nearDistances]=find_nearest_faces(CurrSeed,P_L);
				% get minimum distance for weighting
				minDist=min(nearDistances);
				% get the 3 vector fields, original TR + 30 timepoints (-1 because 1 + 1 is iter 1)
				threeVFs=segTS_l(nearFaces,:,(TRP+t2-1));
				% get a weighting vector
				Wvec1=minDist/nearDistances(1);
				Wvec2=minDist/nearDistances(2);
				Wvec3=minDist/nearDistances(3);
				CorrectedVec1=threeVFs(1,:)*Wvec1;
				CorrectedVec2=threeVFs(2,:)*Wvec2;
				CorrectedVec3=threeVFs(3,:)*Wvec3;	
				% get a weighted average directionality across the 3 faces
				CorrectedAvg=(CorrectedVec1+CorrectedVec2+CorrectedVec3)/(Wvec1+Wvec2+Wvec3);
				% propagate ongoing seeds one step
				CurrSeed=CurrSeed+CorrectedAvg
				% find which vertex the resultant propagating seed is at
				nearest_vertex=find_nearest_vertex(CurrSeed,vx_l)
				% add them to adjacency row at (v,site of ongoing seed)
				adjacency_row(nearest_vertex)=adjacency_row(nearest_vertex)+1;
			end
		end	
	end
	% add adjacency row to adjacency matrix
	AdjMatrix_L(v,:)=adjacency_row;
end
delete(pool);
% save out matrices for this participant
fn=[childfp 'SimStreams/' subj '_' sesh '_' task '_streamConnectivity_L.mat'];
save(fn,'AdjMatrix_L','-v7.3');
