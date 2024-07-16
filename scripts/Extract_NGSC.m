function Extract_NGSC(subj,sesh,task)
% set parent directory
parentfp=['/scratch/groups/leanew1/xcpd_outP50_36p_bp/xcp_d/' subj '/' sesh '/func'];

% define some paths 
Paths{1} = '/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(Paths{1}))

% re-adjust for rsfmri naming conventions
if string(task)=="rs1"
	fp=[parentfp '/' subj '_' sesh '_task-rs_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii'];
	C=ft_read_cifti_mod(fp);
elseif string(task)=="rs2"
	fp=[parentfp '/' subj '_' sesh '_task-rs_acq-mb_dir-pe1_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii'];
	C=ft_read_cifti_mod(fp);
else
% read in time series
C=ft_read_cifti_mod([parentfp '/' subj '_' sesh '_task-' task '_acq-mb_dir-pe0_run-0_space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii']);
end

% extract time series
C_timeseries=C.data;

% subject to SNR Mask
% ALREADY BAKED INTO NMF COMPONENTS
%mwAndTSNR_L='/oak/stanford/groups/leanew1/users/apines/fs4surf/lh.Mask_SNR.func.gii';
%mwAndTSNR_R='/oak/stanford/groups/leanew1/users/apines/fs4surf/rh.Mask_SNR.func.gii';
%mwAndTSNR_L=gifti(mwAndTSNR_L).cdata(:,1);
%mwAndTSNR_R=gifti(mwAndTSNR_R).cdata(:,1);
% 0's where invalid, 1 where valid
%mw_L=ones(1,2562);
%mw_L(mwAndTSNR_L==1)=0;
%mw_R=ones(1,2562);
%mw_R(mwAndTSNR_R==1)=0;

% load in DMN
DMN=ft_read_cifti_mod('/oak/stanford/groups/leanew1/users/apines/maps/Network1_fslr.dscalar.nii');
%DMN_L=gifti('/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/Group_lh_Network_1_32k.func.gii');
%DMN_R=gifti('/oak/stanford/groups/leanew1/users/apines/data/Atlas_Visualize/Group_rh_Network_1_32k.func.gii');
%% 1st component to select DMN.
%Dnet_LH=DMN_L.cdata;
%Dnet_RH=DMN_R.cdata;

% load into gifti'd gordon parcels
%Gord_L=gifti('/oak/stanford/groups/leanew1/users/apines/maps/Gordon_333_left.label.gii').cdata;
%Gord_R=gifti('/oak/stanford/groups/leanew1/users/apines/maps/Gordon_333_right.label.gii').cdata;

% get average dmn value for each parcel
%pDMN=zeros(1,334);
%for p=1:334
	% if it's a left hemi parcel
%	if sum(Gord_L==p)>0;
%		pDMN(p)=mean(Dnet_LH(Gord_L==p));
%	elseif sum(Gord_R==p)>0;
%		pDMN(p)=mean(Dnet_RH(Gord_R==p));
%	end
%end

% get gordon parcels where average DMN value is .3 or greater
%DMNParcels=find(pDMN>.3);
%cxDMN=[];
DMNInds=find(DMN.data>.3);

% load in temporal mask
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh];
tmaskfp=[childfp '/' subj '_' sesh '_task-' task '_AllSegments.txt'];
tmask=load(tmaskfp);
% just get ts length from fmriprep outputs to be safe/consistent
confFilepath1=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subj '/' sesh '/func/' subj '_' sesh '_task-' task '_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv'];
% adapt if it is resting state
if string(task)=="rs1"
	confFilepath1=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subj '/' sesh '/func/' subj '_' sesh '_task-rs_acq-mb_dir-pe0_run-0_desc-confounds_timeseries.tsv'];
elseif string(task)=="rs2"
	confFilepath1=['/oak/stanford/groups/leanew1/SHARED_DATASETS/private/p50/bids/data/derivatives/fmriprep-20.2.3/fmriprep/' subj '/' sesh '/func/' subj '_' sesh '_task-rs_acq-mb_dir-pe1_run-0_desc-confounds_timeseries.tsv'];
else
end
conf1=readtable(confFilepath1,"FileType","text",'Delimiter', '\t');
% extract FD columns
FD=table2array(conf1(:,'framewise_displacement'));

% make binary mask for continuous segments
TRwise_mask_cont=zeros(1,length(FD));
% Loop through each row in Absolut
for row = 1:size(tmask, 1)
	if tmask(row, 3) == 1
        	% Extract the start and end values from the current row
                startValue = tmask(row, 1);
                endValue = tmask(row, 2);
                % change TRwise_mask_cont to 1 where this sequence of continuous good TRs occurs
                TRwise_mask_cont(startValue:endValue)=1;
        else
        end
end

% get complexity of full time series
dmn_ts=C_timeseries(DMNInds,logical(TRwise_mask_cont));
% explicitly using Josh's code, note scrubbing mask is TRwise_mask_cont
[~,~,~,~,EXPLAINED]=pca(dmn_ts);
EXPLAINED=EXPLAINED/100;
nGSC=-sum(EXPLAINED .* log(EXPLAINED))/log(length(EXPLAINED));
cxDMN = nGSC;

% save out normalized entropy for dmn
avComplexity=cxDMN;

T=table(avComplexity,'RowNames',"Row1");
outFP=['/scratch/users/apines/data/mdma/' subj '/' sesh];
writetable(T,[outFP '/' subj '_' sesh '_' task '_Complexity_gro.csv'],'WriteRowNames',true)


