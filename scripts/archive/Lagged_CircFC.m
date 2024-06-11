function CircFC(subj,sesh,task)
%%%% Derive Circular Correlation Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
% Load in fsav4 opflow calc
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
datafp=[childfp '/' subj '_' sesh '_' task '_OpFl.mat'];
data=load(datafp);
%%% L as left Hemi
L=data.us.vf_left;
%%% R as right
R=data.us.vf_right;
%%% get length of time serires
tsLength=length(L);
% Load in surface data
surfL = ['/oak/stanford/groups/leanew1/users/apines/surf/lh.sphere'];
surfR = ['/oak/stanford/groups/leanew1/users/apines/surf/rh.sphere'];
% surface topography
[vx_l, faces_l] = read_surf(surfL);
[vx_r, faces_r] = read_surf(surfR);
% +1 the faces: begins indexing at 0
faces_l = faces_l + 1;
faces_r = faces_r + 1;
% faces_L
F_L=faces_l;
% vertices V
V_L=vx_l;
% vector fields
vfl=data.us.vf_left;
% faces_R
F_R=faces_r;
% vertices V
V_R=vx_r;
% vector fields
vfr=data.us.vf_right;
% get incenters of triangles
TR_L = TriRep(F_L,V_L);
P_L = TR_L.incenters;
TR_R = TriRep(F_R,V_R);
P_R = TR_R.incenters;
% get temporal length
lenOpFl=length(vfl);
% translate xyz spherical coordinates to az/el/r
[az_L,el_L,r_L]=cart2sph(P_L(:,1),P_L(:,2),P_L(:,3));
[az_R,el_R,r_R]=cart2sph(P_R(:,1),P_R(:,2),P_R(:,3));
% convert from radians to degrees
azd_L=rad2deg(az_L);
eld_L=rad2deg(el_L);
azd_R=rad2deg(az_R);
eld_R=rad2deg(el_R);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use native freesurfer command for mw mask indices
surfML = '/oak/stanford/groups/leanew1/users/apines/surf/lh.Medial_wall.label';
mwIndVec_l = read_medial_wall_label(surfML);
surfMR = '/oak/stanford/groups/leanew1/users/apines/surf/rh.Medial_wall.label';
mwIndVec_r = read_medial_wall_label(surfMR);
% make binary "is medial wall" vector for vertices
mw_L=zeros(1,2562);
mw_L(mwIndVec_l)=1;
mw_R=zeros(1,2562);
mw_R(mwIndVec_r)=1;
% convert to faces
% convert to faces
F_MW_L=sum(mw_L(faces_l),2)./3;
F_MW_R=sum(mw_R(faces_r),2)./3;
% convert "partial" medial wall to medial wall
F_MW_L=ceil(F_MW_L);
F_MW_R=ceil(F_MW_R);
% face mask indices
fmwIndVec_l=find(F_MW_L);
fmwIndVec_r=find(F_MW_R);
% make medial wall vector
g_noMW_combined_L=setdiff([1:5120],fmwIndVec_l);
g_noMW_combined_R=setdiff([1:5120],fmwIndVec_r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize output arrays
AnglesL=zeros(length(g_noMW_combined_L),lenOpFl);
AnglesR=zeros(length(g_noMW_combined_R),lenOpFl);
CFC_L=zeros(5120);
CFC_R=zeros(5120);
% master CFC, across segments
mCFC_L=zeros(5120);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load in continuous segment indices
CSIfp=[childfp '/' subj '_' sesh '_task-' task '_ValidSegments_Trunc.txt'];
CSI = importdata(CSIfp);
% assure that TR count is the same between time series and valid segments txt
SegNum=size(CSI);
SegNum=SegNum(1);
% number of trs in fmri ts
of_ts_trs=length(L);
% trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
numTRsVS=CSI(SegNum,1)+CSI(SegNum,2)-1;
% note we lose a "between" volume for each segment
if numTRsVS ~= (of_ts_trs+SegNum)
        disp('TRs from Valid Segments txt and derived vectors do not match. Fix it.')
        return
end
% initialize full-length counter
fullLength=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note loop opening without indentation
% loop over each segment and find delay assoc. with max cor, max corr %%%%%@@@@@@@@@@@@@@@2 
for seg=1:SegNum;
	SegStart=CSI(seg,1);
	SegSpan=CSI(seg,2);
	% need to limit to segments over 14 for this script
	if SegSpan>14
	% extract optical flow for this segment
	% indexing is a little complicated, because we are transforming from TRwise to betweenTRwise. -1 for each segment, but starts at 1
	vfl=L(((SegStart+1)-seg):(((SegStart)-seg)+SegSpan-1));
	% get length of optical flow in this segment
	lenOpFl=SegSpan-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FOR ALL LEFT FACES %%%%%%
%Get az/el, derive circ corrs%
% I. Az/el
% note we are loopign over the non-mw-face indices, skipping medial wall faces but leaving 0's in their stead
for F=g_noMW_combined_L
        % translate xyz vector fields from opfl to az/el/r to polar angles
        % note, by not saving rho (just theta), we are discarding magnitude information at this point
        for fr=1:lenOpFl
                % current vector field
                relVf_L=vfl{fr};
                % xyz components
                xComp_L=relVf_L(F,1);
                yComp_L=relVf_L(F,2);
                zComp_L=relVf_L(F,3);
                % convert to spherical coord system
                vs_L=cart2sphvec(double([xComp_L;yComp_L;zComp_L]),azd_L(F),eld_L(F));
                % store in output vector (r is redundant across all vecs, only using az and el)
                AnglesL(F,fr)=cart2pol(vs_L(1),vs_L(2));
        end
end	
% II. derive circ corrs
% get triu, loop over triu coordinates
[i,j]=meshgrid(g_noMW_combined_L,g_noMW_combined_L);
UpperCoords_x=i(i<j);
UpperCoords_y=j(i<j);
% loop over upper triangle of adjacency matrix
for F=1:length(UpperCoords_x);
	% initialize CFC_candidate
	top_CFC=0;
	% get x and y coord
	Xcoord=UpperCoords_x(F);
	Ycoord=UpperCoords_y(F);
	% inquire across possible delays 
	for delay=-4:4
		% Determine the start and end indices for the overlapping portion
		start_index = max(1, 1 - delay);
		end_index = min(lenOpFl, lenOpFl - delay);
		% Calculate the overlap between the two time series after applying the delay.
		% For 'truncated_series_x', we start at the maximum of the beginning of the series
		% and the beginning of the delayed series, and end at the minimum of the series
		% length and the length of the delayed series.
		truncated_series_x = AnglesL(Xcoord,start_index:end_index);
		% To calculate the corresponding portion of 'truncated_series_y', we start at the
		% same 'start_index', but since 'truncated_series_y' is delayed by 'delay' time steps,
		% we adjust the end_index by adding 'delay' to it.
		truncated_series_y = AnglesL(Ycoord,start_index + delay:end_index + delay);
		[CFC_candidate pval] = circ_corrcc(truncated_series_x,truncated_series_y);
		% if it's higher than max coord, take it!
		if abs(CFC_candidate)>top_CFC;
			top_CFC=CFC_candidate;
			CFC_L(Xcoord,Ycoord)=top_CFC;
	 	else
		end
	end		
end	
% end contingency of if segment is >14 TRs
else
end
% multiply by number of between TRs in this segment
CFC_L=CFC_L.*(SegSpan-1);
% add to master CFC matrix
mCFC_L=mCFC_L+CFC_L;
% update full length counter
fullLength=fullLength+(SegSpan-1);
% end for each segment
end
% divide by total number of TRs
mCFC_L=mCFC_L./fullLength;
% mirror
mCFC_L=triu(mCFC_L)+triu(mCFC_L,1)';
CFC_R=triu(CFC_R)+triu(CFC_R,1)';
% mask mw 0s out
mCFC_L=mCFC_L(g_noMW_combined_L,g_noMW_combined_L);
%CFC_R=CFC_R(g_noMW_combined_R,g_noMW_combined_R);
% saveout
%adjMats=struct('L',{CFC_L},'R',{CFC_R});
adjMats=struct('L',{mCFC_L});
fn=['/scratch/users/apines/data/mdma/' subj '/' subj '_' sesh '_' task '_CircFC.mat'];
save(fn,'adjMats')
