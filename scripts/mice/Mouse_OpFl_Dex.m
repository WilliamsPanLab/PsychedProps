function Mouse_OpFl_Dex(subj,run)
% all credit to NeuroPattToolbox: https://github.com/BrainDynamicsUSYD/NeuroPattToolbox, at least Rory Townsend, Xian Long, and Pulin Gong
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/NeuroPattToolbox'))
% Most of code is from:
% Rory Townsend, Aug 2018
% rory.townsend@sydney.edu.au
% set sampling frequency
Fs=15;
params=struct;
params.zscoreChannels=0;
params.params.subtractBaseline=0;
% enforcement of smoothness of vector field
params.opAlpha = 0.5;
startTime=datetime
% https://github.com/rorygt/NeuroPattToolbox/blob/master/setParams.m claims opBeta=0.01 is default, but it's not clear we can run the code with that low of a beta and 1 is what is actually set in the code. Trying 0.1 to afford some nonlinearity while being computationally tractable
% enforcement of linearity of vector field (higher = more linear)
params.opBeta = 0.1;
% note we don't want to analyze phase directly in case propagations of interest are aperiodic
params.useAmplitude = true;

% AP load in data: pre dex
basefp='/oak/stanford/groups/leanew1/users/apines/p50_mice/proc2/proc/drugs/processed_data/20180404/'
% different date for m1
if subj==string('m1')
	basefp='/oak/stanford/groups/leanew1/users/apines/p50_mice/proc2/proc/drugs/processed_data/20180405/'
end
% load in specified scan
if run==1	
	fn = [basefp 'thy1gc6s_140mgkg_' subj '_preDex_1/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
	if exist(fn)
		data=h5read(fn, '/processed_data');
		formask=h5read(fn, '/mask');
	else
		disp('no run found')
	end	
end
%% add if/else
% NO POST_0 FOR DEX
% post 1
%if run==2
%	fn = [basefp 'thy1gc6s_140mgkg_' subj '_postDex_0/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
%	if exist(fn)
%		data=h5read(fn, '/processed_data');
%		formask=h5read(fn, '/mask');
%	else
%		disp('no run found')
%	end
%end
% post 2
if run==3	
	fn = [basefp 'thy1gc6s_140mgkg_' subj '_postDex_5/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
	if exist(fn)
		data=h5read(fn, '/processed_data');
		formask=h5read(fn, '/mask');
	else    
                disp('no run found')
        end
end 
% post 3
if run==4
	fn = [basefp 'thy1gc6s_140mgkg_' subj '_postDex_10/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
	if exist(fn)
		data=h5read(fn, '/processed_data');
		formask=h5read(fn, '/mask');
	else
		disp('no run found')
	end
end
% post 4
if run==5
	fn = [basefp 'thy1gc6s_140mgkg_' subj '_postDex_15/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
	if exist(fn)
		data=h5read(fn, '/processed_data');
		formask=h5read(fn, '/mask');
	else
		disp('no run found')
	end
end
% post 5
if run==6
	fn = [basefp 'thy1gc6s_140mgkg_' subj '_postDex_20/masked_dff_Gro_Masked_Sml_BP_Smoothed_Sml.h5']
	if exist(fn)
		data=h5read(fn, '/processed_data');
		formask=h5read(fn, '/mask');
	else	
		disp('no run found')
	end
end
% AP - commenting out and replacing with mask derived from atlas
mask = cellfun(@(x) strcmp(x, 'TRUE'), formask); 
badChannels=mask;
% Find any invalid channels that must be interpolated over
% Any channels with any NaN values
%nanChans = any(isnan(data(:,:,:)),3);
% Any channels that never change over time
%zeroChans = all(data(:,:,:)==0, 3);
%badChannels = find(nanChans | zeroChans);
% AP - replacing preproc with manual preproc (in previous script)
timeDim = 3;
ntimesteps = size(data, timeDim) - 1;
wvcfs=data;

% Pre-allocate velocity field variables
vfs = zeros(size(wvcfs));
vfs = vfs(:,:,1:end-1,:);
meanCSteps = zeros(size(wvcfs,4), 1);

% Calculate velocity fields for every trial: consider feeding in pre drug as trial one, and 4 post durgs as trials 1-5!
% steps to converge
% AP - Only one trial in this version
for itrial = 1
%for itrial = 1:size(wvcfs,4)
    [vx, vy, csteps] = opticalFlow2(wvcfs(:,:,:,itrial), badChannels, ...
        params.opAlpha, params.opBeta, ~params.useAmplitude);
    % alternatively, you can use Matlab built-in function 'opticalflowHS', 
    % which does not provide the penalty term Beta but is optimised for speed.
    % alpha = 1 ;
    % opticFlow = opticalFlowHS('Smoothness',alpha,'MaxIteration',100) ;
    % Vx = zeros(size(wvcfs,1),size(wvcfs,2),size(wvcfs,3)-1) ;
    % Vy = zeros(size(wvcfs,1),size(wvcfs,2),size(wvcfs,3)-1) ;
    %  for iTime = 1:downR:size(sigIn,3)
    %        flow = estimateFlow(opticFlow,wvcfs(:,:,iTime,itrial));
    %        Vx(:,:,iTime) = flow.Vx ;
    %        Vy(:,:,iTime) = flow.Vy ;
    %  end
    % not sure why we would want vfs - AP, maybe conversion to circular coordinate system is used for further derivatives?
    vfs(:,:,:,itrial) = vx + 1i*vy;

    meanCSteps(itrial) = mean(csteps);
    fprintf('Processed trial %i\n', itrial)
end
%
%%% AP: skipping SVD could be useful at some point
%
% this is their pattern estimation, looks cool
tic
% Set up pattern structures
allPatts = cell(1, size(wvcfs,4));
allLocs = allPatts;
% Loop over all repetitions to find patterns
for itrial = 1:size(wvcfs,4)
    thisvf = vfs(:,:,:,itrial);
    [patterns, patternTypes, patternColNames, pattLocs] = ...
        findAllPatterns(real(thisvf), imag(thisvf), params, ...
        angle(wvcfs(:,:,:,itrial)));
    allPatts{itrial} = patterns;
    allLocs{itrial} = pattLocs;
end

% evaluating pattern evolution, also looks cool
%% Analyse evolution between patterns
% Number of time steps before and after a pattern ends to search for other
% patterns
nafter = round(0.05*Fs);
nbefore = round(0.01*Fs);

pattTypeStr = '';
for itype = 1:length(patternTypes)
    pattTypeStr = sprintf('%s%i.%s ', pattTypeStr,itype,patternTypes{itype});
end

[nobs, nexp] = pattEvolution(allPatts, ntimesteps, nafter, nbefore);
rateDiff = (nobs - nexp) / ntimesteps * Fs;
%disp('Observed minus expected pattern transitions/sec')
disp(pattTypeStr)
%disp(nanmean(rateDiff,3))
%disp(median(nobs,3) - median(nexp,3))

disp('Fractional change between observed and expected')
disp(nanmean((nobs-nexp)./nexp, 3));


% AP - omitting tests across trials: runs too large to fit into one run anyway. Would have to run post-hoc after loading in derived vector fields

% Set all outputs
% AP omitting onlyPatterns: we want dem vecta fieldz
outputs.filteredSignal = wvcfs;
outputs.velocityFields = vfs;
outputs.badChannels = badChannels;
outputs.nTimeSteps = ntimesteps;
outputs.patternTypes = patternTypes;
outputs.patternResultColumns = patternColNames;
outputs.patterns = allPatts;
outputs.patternLocs = allLocs;
outputs.params = params;
outputs.pattTransitionsObs = nobs;
outputs.pattTransitionsExp = nexp;
outputs.Fs = Fs;
outputs.processTime = datetime - startTime;

% save outputs to mouse dir
outfp=[basefp subj '_vf_out_' num2str(run) '.mat'];
save(outfp,'outputs');
% visualize
useAmplitude= true;
outFoldName=['/scratch/users/apines/mouseViz/' subj '_' num2str(run)];
system(['mkdir ' outFoldName]);
vidName=[outFoldName '/vecField_'];
vidFps=15;
resizeScale=1;
vfScale=1;
saveVelocityFieldVideo(wvcfs, vfs, vidName, vidFps, ...
    Fs, resizeScale, vfScale, useAmplitude,'_Gen')
