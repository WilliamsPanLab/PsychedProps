function Extract_DMNMag_Spun_mice(subj,sesh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take optical flow results, get a bottom-up and top-down resultant vector in x,y coords for each pixel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));
% set number of spins
numSpins=2000;
% Load in flatmouse opflow calc
% note this is for LSD only! Adapt recording date for other drugs if needed
childfp='/oak/stanford/groups/leanew1/users/apines/p50_mice/proc2/proc/20200228/'
datafp=[childfp subj '_vf_out_' num2str(sesh) '.mat']
% adding in if it exists: ends all the way at the end of the script
if exist(datafp)
data=load(datafp)
% vector fields
vf=data.outputs.velocityFields;
% and actual signal
signalGrid=data.outputs.filteredSignal;
% load in spun DMN values
networks_OG=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
networks=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_permutations_10k.mat');
% initialize cross-spin vector
SpunMags=zeros(1,numSpins);
% for each spin
for s=1:numSpins
	Dnet=networks.DMN_perms(:,:,s);
        Mask=networks_OG.mask;
        % need to re-boolean from the interpolation in the downsample
        Mask=Mask==1;
        % need to transpose mask to match python-matlab difference: to be confirmed via visualization
        Mask=Mask';
        % get size of mask and indices to include
        MaskSize=sum(sum(Mask));
        MaskInds=find(Mask);
        % extract size of time series
        NumTRs=size(vf);
        NumTRs=NumTRs(3);
        lenOpFl=NumTRs;
        % and an "over time" version of the mask to apply to the time series
        fullMask = repmat(Mask, [1, 1, lenOpFl]);

        % initialize out dataframes
        Propvec=[];
        stringVec={};

        % and azez and els for opflow vectors
        OpF_x=zeros(size(vf,1),size(vf,2),lenOpFl);
        OpF_y=zeros(size(vf,1),size(vf,2),lenOpFl);
        for tp=1:lenOpFl
                OpF_x(:,:,tp)=real(vf(:,:,tp));
                OpF_y(:,:,tp)=imag(vf(:,:,tp));
        end
        % mask out angles
        OpF_x(~fullMask)=0;
        OpF_y(~fullMask)=0;
        % create pixel-wise network mask
        DMN_bool=Dnet;
        DMN_bool(DMN_bool>.6)=1;
        DMN_bool(DMN_bool<.6)=0;
        % adding mask into this step 7/1/24
        DMN_bool(Mask==0)=0;
        DMN_bool=logical(DMN_bool);
        % initialize matrix for each pixel to saveout to scratch
        faceMatrix=zeros(sum(sum(DMN_bool)));
        % network of interest
        net=Dnet	

        % initialize angular distance vector for each network (l and r) above
        Mags=zeros(sum(sum(DMN_bool)),lenOpFl);
	% get angular distance for each face for each timepoint
        for F=1:sum(sum(DMN_bool));
                % loop over each tp
                for fr=1:lenOpFl
			% pull out optical flow vectors for this pixel (denoted by F, because it represents faces in human workflow)
			curOpF_x=OpF_x(:,:,fr);
			curOpF_y=OpF_y(:,:,fr);
			curOpF_x_DMN=curOpF_x(DMN_bool);
			curOpF_y_DMN=curOpF_y(DMN_bool);
                        % get optical flow vector
                        OpFlVec=[curOpF_x_DMN(F) curOpF_y_DMN(F)];
			% calculate magnitude!
                        Mags(F,fr)=sqrt(sum(OpFlVec .^2));

                % end tp loop
                end
	% end each face loop
        end
	avgD=mean(mean(Mags));
	SpunMags(s)=avgD;
end
stringVec = compose("Spin%d", 1:2000);
% calc outFP
outFP=['/scratch/users/apines/data/mouse/'];
% add label
T=table(SpunMags','RowNames',stringVec);
% write out
writetable(T,[outFP subj '_' num2str(sesh) '_Prop_Feats_Mag_Spun.csv'],'WriteRowNames',true)
else
	disp('file not found')
end
