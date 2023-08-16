function partialC_and_BP_and_TDif_ISPOT(subj)
% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))
% define filter
% sampling frequency per second
F=1/.8
% nyquist
Ny=F/2
% lower and upper limit
BW_Hz = [.009 .08];
bp_order=2;
BW_N = BW_Hz/Ny;
[b, aa] = butter(ceil(bp_order/2),BW_N);
% task list
tlist=["gng","noncon","con"];
% init avg. valuevec
LH_avgValueVec=zeros(2562,3);
RH_avgValueVec=zeros(2562,3);
% for each task
for t=1:3
	task=tlist(t)
	%%% LEFT
	% input image filepath
	imagefp=char(strjoin(['/scratch/users/apines/data/ispot/' string(subj) '/' task '_lh.mgh'],''));
	imagefp_out=char(strjoin(['/scratch/users/apines/data/ispot/' string(subj) '/' task '_lh_pcent.mgh'],''));
	% load image
	image=MRIread(imagefp);
	% convert to true dimensionality
	series=squeeze(image.vol);
	% commented out for task contrast	
	% normalize each vertex  to be mean centered
	% series=normalize(series,2);
	% normalize so each timepoint is mean centered
	normseries=normalize(series,1);
	% need to transpose to filter over time rather than arbitrary spatial dimension
	Y = filtfilt(b,aa,double(normseries'));
	image.vol(1,:,1,:)=Y';
	LH_avgValueVec(:,t)=mean(Y',2);   
	MRIwrite(image,imagefp_out)
	%%% RIGHT
	% input image filepath
        imagefp=char(strjoin(['/scratch/users/apines/data/ispot/' string(subj) '/' task '_rh.mgh'],''));
        imagefp_out=char(strjoin(['/scratch/users/apines/data/ispot/' string(subj) '/' task '_rh_pcent.mgh'],''));
	% load image
        image=MRIread(imagefp);
        % convert to true dimensionality
        series=squeeze(image.vol);
        % commented out for task contrast
	% normalize each vertex  to be mean centered
        % series=normalize(series,2);
        % normalize so each timepoint is mean centered
        normseries=normalize(series,1);
        % need to transpose to filter over time rather than arbitrary spatial dimension
        Y = filtfilt(b,aa,double(normseries'));
        image.vol(1,:,1,:)=Y';
	RH_avgValueVec(:,t)=mean(Y',2);
        MRIwrite(image,imagefp_out)
end

% just take simple difference
LH_GNG_NCF_dif=LH_avgValueVec(:,1)-LH_avgValueVec(:,2);
RH_GNG_NCF_dif=RH_avgValueVec(:,1)-RH_avgValueVec(:,2);

% convert to csv
csv_L=table(LH_GNG_NCF_dif(:,:));
csv_R=table(RH_GNG_NCF_dif(:,:));
writetable(csv_L,['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/' subj '/GNG-NCF_L.csv'])
writetable(csv_R,['/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/ispot/' subj '/GNG-NCF_R.csv'])
