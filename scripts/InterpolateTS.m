function InterpolateTS(subj,sesh)
% this function is to interpolate the functional time series to within-segment and between-TR timepoints. Should be v stable in band-passed fMR signal
% load in time series (downsampled)
childfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];
% load in data
fpL=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_L.mgh'];
fpR=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_R.mgh'];
dataL=MRIread(fpL);
dataR=MRIread(fpR);
% squeeze to get rid of extra dimensions
TRs_l_g=squeeze(dataL.vol);
TRs_r_g=squeeze(dataR.vol);
% get sizes
sizefMR=size(TRs_l_g);
lengthTS=sizefMR(2);
% load in continuous segment information
CSIfp=[childfp '/' subj '_' sesh '_task-rs_ValidSegments_Trunc.txt'];
CSI = importdata(CSIfp);
% assure that TR count is the same between time series and valid segments txt
SegNum=size(CSI);
SegNum=SegNum(1);
% trailing -1 is because the count column (,2) is inclusive of the start TR (,1)
numTRsVS=CSI(SegNum,1)+CSI(SegNum,2)-1;
if numTRsVS ~= lengthTS
        disp('TRs from Valid Segments txt and mgh file do not match. Fix it.')
        return
end
% store original timepoints as sequence
originalTimepoints = 1:lengthTS; 
%%% derive segment-ending TRs from CSI to factor them out
% get the number that each segment starts with
SEtrs=CSI(:,1);
% take out first starting point (i.e., first TR) from this exclusion list. Only seeking to omit interpolation between book-ends of segments
SEtrs=SEtrs(2:end);
% and finally, -1 to get 1 prior to starting point
SEtrs=SEtrs-1;
% now omit these TRs from originalTimepoints in terms of what TRs we want to interpolate between
% AKA remove segment-ending TRs from this sequence to not interpolate between-segment TRs
timepointsToInterpolate = setdiff(originalTimepoints, SEtrs);
% Calculate halfway timepoints
halfwayTimepoints = timepointsToInterpolate(1:end-1) + 0.5;
% initialize interpolated data
InterpData_L=zeros(sizefMR(1),length(halfwayTimepoints));
InterpData_R=zeros(sizefMR(1),length(halfwayTimepoints));
% Interpolate data at halfway timepoints
for v=1:10242
	interpolatedData_L(v,:) = interp1(originalTimepoints, TRs_l_g(v,:), halfwayTimepoints);
	interpolatedData_R(v,:) = interp1(originalTimepoints, TRs_r_g(v,:), halfwayTimepoints);
end
% save it out
ofpl=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_interp_L.mgh'];
ofpr=[childfp '/' subj '_' sesh '_task-rs_p2mm_masked_interp_R.mgh'];
% saveout
dataL.vol=interpolatedData_L;
dataR.vol=interpolatedData_R;
MRIwrite(dataL,ofpl)
MRIwrite(dataR,ofpr)
