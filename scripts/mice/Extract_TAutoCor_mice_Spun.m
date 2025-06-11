function Extract_TAutoCor_mice(subj,run)
% set parent directoryP load in data: pre LSD
basefp='/oak/stanford/groups/leanew1/users/apines/p50_mice/proc2/proc/20200228/'
% load in specified scan
if run==1
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_preLSD0p3mgkg_1/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end
%% add if/else
% post 1
if run==2
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_0/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end
% post 2
if run==3
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_5/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end
% post 3
if run==4
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_10/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end
% post 4
if run==5
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_15/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end
% post 5
if run==6
        fn = [basefp 'thy1gc6s_0p3mgkg_' subj '_postLSD0p3mgkg_20/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
        if exist(fn)
                data=h5read(fn, '/processed_data');
        else
                disp('no run found')
        end
end

% pull in data as "C", refers to cifti in human scripts
C_timeseries=data;
lenOpFl=size(C_timeseries,3);

% set number of spins
numSpins=2000;
% load in spun DMN values
networks_OG=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_components.mat');
networks=load('/oak/stanford/groups/leanew1/users/apines/p50_mice/Mouse_DMN_permutations_10k.mat');
% initialize cross-spin vector
SpunTAs=zeros(1,numSpins);
% for each spin
for s=1:numSpins
	% load in mouse DMN
	Dnet=networks.DMN_perms(:,:,s);
	Mask=networks_OG.mask;
	% ensure boolean
	Mask=Mask==1;
	% transpose
	Mask=Mask';
	% Look in spun DMN
	Mask=Mask(Dnet>.6);
	% initialize DMN time series
	DMN_mat = zeros(sum(Mask(:)), lenOpFl);
	% use for loop to straightforwardly extract time series in pixels of interest
	for t = 1:lenOpFl
	    ts_slice = C_timeseries(:, :, t);
	    DMN_mat(:, t) = ts_slice(Mask);
	end

	% initialize an autocor array
	ACarray=zeros(1,sum(Mask(:)));
	% number of pixels
	NumPix=sum(Mask(:));

	% loop over each pixel within DMN
	for p = 1:NumPix;
		% extract data
		dataInPix=DMN_mat(p,:);
		% OG time series (-1)
		OGTS=dataInPix(1:end-1);
        	ShiftTS=dataInPix(2:end);
		% calculate autocor (t-1 vs. t)
		correlationOfInterest=corrcoef(OGTS,ShiftTS);
        	% note this returns a correlation matrix: only interested in the between-timeseries measurement (off-diag)
        	ACarray(p)=correlationOfInterest(1,2);
	% end pixel loop
	end

	% get average
	av_AutoCor=mean(ACarray);
	SpunTAs(s)=av_AutoCor;
end
stringVec = compose("Spin%d", 1:2000);
T=table(SpunTAs','RowNames',stringVec);
% makedir
outFP=['/scratch/users/apines/data/mouse/' subj '/'];
% makedir
writetable(T,[outFP '/' subj '_' num2str(run) '_av_AutoCor_Spun.csv'],'WriteRowNames',true)


