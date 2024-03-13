%000000000 Need to adapt all paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

%%%OOOOOO load in placebo subj/sesh vs. drug subj/sesh instead of rs 
Subjs=readtable('~/PWs/rs_subs.csv')

% initialize output array

% 18 bins for angular distances because 0 and 180 are endpoints (19 numbers)
OutDf=zeros(1,height(Subjs));
% 17 bins for modes because there is no 0th bin - 1 and 18 are endpoints (18 numbers)
% edges to apply in discretize
edges = 0:10:180;

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

%%% use TSNR mask, includes medial wall
mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
mw_L=zeros(1,2562);
mw_L(mwAndTSNR_L==1)=1;
mw_R=zeros(1,2562);
mw_R(mwAndTSNR_R==1)=1;
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

% loop over each subj
for s = 1:height(Subjs)
	Subj=table2array(Subjs(s,2));
	fp=['/cbica/projects/pinesParcels/results/PWs/Proced/' Subj{:}];
	% if file exists, then
	if isfile([fp '/' Subj{:} '_AngDistMat4.mat']);
		s
		% load in subj's distr

		% potentially faceMatrix,['/scratch/users/apines/gp/PropFeats/' subj '_' sesh '_' task '_k' num2str(k) '_faceMatrix_gro.csv'] ?

		Angs=load([fp '/' Subj{:} '_AngDistMat4.mat']);
		AngsL=Angs.AngDist.gLeft;
		AngsR=Angs.AngDist.gRight;
		% initialize facewise vectors for 1-18 mode for 
		facePropL=zeros(1,length(g_noMW_combined_L));
		facePropR=zeros(1,length(g_noMW_combined_R));
		% for each face: get mode over TRs
		for l=1:length(g_noMW_combined_L);
			binnedFaceDirs=histcounts(AngsL(:,g_noMW_combined_L(l)),edges);
			BUs=sum(binnedFaceDirs(1:9));
			TDs=sum(binnedFaceDirs(10:18));
			facePropL(l)=BUs/(BUs+TDs);	
		end
		% for right
		for r=1:length(g_noMW_combined_R);
                        binnedFaceDirs=histcounts(AngsR(:,g_noMW_combined_R(r)),edges);
                        BUs=sum(binnedFaceDirs(1:9));
                        TDs=sum(binnedFaceDirs(10:18));
                        facePropR(r)=BUs/(BUs+TDs);
		end
		% count of prop BU for this subj
		OutDf(s)=mean([facePropL facePropR]);
	end
end

% maybe write as csv: as table ya know
fn=['/cbica/projects/pinesParcels/results/PWs/rs_subsBuProp.csv'];
writetable(table(OutDf),fn)

