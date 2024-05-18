function TestPatt(FB)
% all credit to NeuroPattToolbox: https://github.com/BrainDynamicsUSYD/NeuroPattToolbox, at least Rory Townsend, Xian Long, and Pulin Gong
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/scripts/NeuroPattToolbox'))
% Most of code is from:
% Rory Townsend, Aug 2018
% rory.townsend@sydney.edu.au

% load in patterns for this frequency band
mList={'m2000','m7520','m7522','m7589','m7594'};
% set params
params=struct;
params.zscoreChannels=0;
params.params.subtractBaseline=0;
% enforcement of smoothness of vector field
params.opAlpha = 0.5;
startTime=datetime
% https://github.com/rorygt/NeuroPattToolbox/blob/master/setParams.m claims opBeta=0.01 is default, but it's not clear we can run the code with that low of a beta and 1 is what is actually set in the code. Trying 0.1 to afford some nonlinearity while being computationally tractable
% enforcement of linearity of vector field (higher = more linear)
params.opBeta = 0.1;
% note we don't want to analyze phase directly incase propagations of interest are aperiodic
params.useAmplitude = true;
% initialize number of observations total and expected total
nobs_tot=zeros(7,6);
nexp_tot=zeros(7,6);
% base filepath for LSD
basefp='/scratch/users/apines/p50_mice/proc/20200228/';
% for each mouse
for i=1:length(mList)
	% initialize struct to store vector fields
	m_vfs=zeros(67,70,892,6);
	m=mList{i}
	% load in each run
	for run=1:6
		outfp=[basefp m '_vf_out_' FB '_' num2str(run) '.mat'];
		OFdata=load(outfp);
		% extract vector fields
		vfs=OFdata.outputs.velocityFields;
		m_vfs(:,:,:,run)=vfs;
	end
	% initialize patterns struct
	allPatts = cell(1, size(m_vfs,4));
	allLocs = allPatts;
	% Loop over all repetitions to find patterns
	for itrial = 1:size(m_vfs,4)
	    thisvf = m_vfs(:,:,:,itrial);
	    [patterns, patternTypes, patternColNames, pattLocs] = ...
	        findAllPatterns(real(thisvf), imag(thisvf), params, ...
	        angle(m_vfs(:,:,:,itrial)));
	    allPatts{itrial} = patterns;
	    allLocs{itrial} = pattLocs;
	    % Initialize pattFreq as a 7x1 vector with zeros
	    pattFreq = zeros(7, 1);
	    % Fill in the detected patterns
        	if ~isempty(patterns)
            		% Tabulate the patterns
            		tempPattFreq = tabulate(patterns(:,1))
            		% Update pattFreq with the detected patterns
            		for j = 1:size(tempPattFreq,1)
				currPatts = find(tempPattFreq(:,1)==j);
				freq=tempPattFreq(currPatts,2)
				% insert into cross-mouse df
				nobs_tot(j, itrial) = nobs_tot(j, itrial) + freq;
        		end
		end
	end
	%% Analyse evolution between patterns
	% Number of time steps before and after a pattern ends to search for other
	% patterns
	Fs=15;
	nafter = round(0.05*Fs);
	nbefore = round(0.01*Fs);

	pattTypeStr = '';
	for itype = 1:length(patternTypes)
	    pattTypeStr = sprintf('%s%i.%s ', pattTypeStr,itype,patternTypes{itype});
	end
	% set timesteps
	timeDim = 3;
	ntimesteps = size(m_vfs,timeDim);
	[nobs, nexp] = pattEvolution(allPatts, ntimesteps, nafter, nbefore);
	rateDiff = (nobs - nexp) / ntimesteps * Fs;
	%disp('Observed minus expected pattern transitions/sec')
	disp(pattTypeStr)
	disp(nanmean(rateDiff,3))
	disp(median(nobs,3) - median(nexp,3))

	disp('Fractional change between observed and expected')
	disp(nanmean((nobs-nexp)./nexp, 3));

	% Test differences between observed and expected if multiple trials are
	% present
	if size(m_vfs, 4) > 1

	    pvals = zeros(size(nobs,1));
    	for initPatt = 1:size(nobs,1)
    	    for nextPatt = 1:size(nobs,2)
    	        thisObs = nobs(initPatt, nextPatt, :);
    	        thisExp = nexp(initPatt, nextPatt, :);
    	        [h, p] = ttest(thisObs(:),  thisExp(:));
    	        pvals(initPatt, nextPatt) = p;
    	    end
    	end
    	disp('Paired t-test p-values')
    	fprintf('Bonferroni correction factor = %i\n', numel(pvals))
    	disp(pvals)
	end
	nobs_tot
end
disp('all mice ran')
nobs_tot
