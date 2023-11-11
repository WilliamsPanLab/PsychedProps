function C_and_BP_ISPOT(subj)
% add paths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs'))
% define filter based on acquisition speed
LegacyTR_i=2.5;
LegacyTR_r=2;
MBTR=.71;
% sampling frequency per second
F=1/LegacyTR_i;
% nyquist
Ny=F/2;
% lower and upper limit
BW_Hz = [.009 .08];
bp_order=2;
BW_N = BW_Hz/Ny;
[b, aa] = butter(ceil(bp_order/2),BW_N);
% task list
tlist=["gng","noncon","con"];
% for each task
for t=1:3
	task=tlist(t)
	%%% LEFT
	% input image filepath
	imagefp=char(strjoin(['/scratch/users/apines/data/ispot/' string(subj) '/' task '_lh.mgh'],''));
	imagefp_out=char(strjoin(['/scratch/users/apines/data/ispot/' string(subj) '/' task '_lh_cent.mgh'],''));
	% load image
	image=MRIread(imagefp);
	% convert to true dimensionality
	series=squeeze(image.vol);
	% normalize each vertex  to be mean centered
	series=normalize(series,2);
	% normalize so each timepoint is mean centered
	normseries=normalize(series,1);
	% need to transpose to filter over time rather than arbitrary spatial dimension
	Y = filtfilt(b,aa,double(normseries'));
	image.vol(1,:,1,:)=Y';   
	MRIwrite(image,imagefp_out)
	%%% RIGHT
	% input image filepath
        imagefp=char(strjoin(['/scratch/users/apines/data/ispot/' string(subj) '/' task '_rh.mgh'],''));
        imagefp_out=char(strjoin(['/scratch/users/apines/data/ispot/' string(subj) '/' task '_rh_cent.mgh'],''));
	% load image
        image=MRIread(imagefp);
        % convert to true dimensionality
        series=squeeze(image.vol);
        % normalize each vertex  to be mean centered
        series=normalize(series,2);
        % normalize so each timepoint is mean centered
        normseries=normalize(series,1);
        % need to transpose to filter over time rather than arbitrary spatial dimension
        Y = filtfilt(b,aa,double(normseries'));
        image.vol(1,:,1,:)=Y';
        MRIwrite(image,imagefp_out)
end
