function IntraParcelSynchrony_ot_left(subj,sesh,task)
% addpaths
addpath(genpath('/oak/stanford/groups/leanew1/users/apines/libs/'))
parentfp=['/scratch/users/apines/data/mdma/' subj '/' sesh ];

% load in vertex-wise angles
fpl=[parentfp '/' subj '_' sesh '_task-' task '_p2mm_masked_Vert_Angles_L.mat'];
fpr=[parentfp '/' subj '_' sesh '_task-' task '_p2mm_masked_Vert_Angles_R.mat'];
vertAngles_L=load(fpl).vertWise_Vecs_l;
vertAngles_R=load(fpr).vertWise_Vecs_r;

% load in interppolated time series (interpolated to b/w timepoints, like the angles)
fpl=[parentfp '/' subj '_' sesh '_task-' task '_p2mm_masked_interp_L.mgh'];
fpr=[parentfp '/' subj '_' sesh '_task-' task '_p2mm_masked_interp_R.mgh'];
bold_L=MRIread(fpl).vol;
bold_R=MRIread(fpr).vol;
% get length of scan
lenOpFl=size(bold_L);
lenOpFl=lenOpFl(2);

% Parcellation
schaef100_L=gifti('/oak/stanford/groups/leanew1/users/apines/maps/Gro_L_Schaef100_3k.func.gii');
schaef100_R=gifti('/oak/stanford/groups/leanew1/users/apines/maps/Gro_R_Schaef100_3k.func.gii');

% initialize output df: subj sesh task parcel TR starting point t1 value t2 value ... t25 value
outDF = struct('subj', {}, 'sesh', {}, 'task', {}, 'parcel', {}, 'startTR', {}, 'v1', {}, 'v2', {}, 'v3', {}, 'v4', {}, 'v5', {}, 'v6', {}, 'v7', {}, 'v8', {}, ...
    'v9', {}, 'v10', {}, 'v11', {}, 'v12', {}, 'v13', {}, 'v14', {}, 'v15', {}, 'v16', {}, 'v17', {}, 'v18', {}, ...
    'v19', {}, 'v20', {}, 'v21', {}, 'v22', {}, 'v23', {}, 'v24', {}, 'v25', {}, 't1', {}, 't2', {}, 't3', {}, 't4', {}, 't5', {}, 't6', {}, 't7', {}, 't8', {}, ...
    't9', {}, 't10', {}, 't11', {}, 't12', {}, 't13', {}, 't14', {}, 't15', {}, 't16', {}, 't17', {}, 't18', {}, ...
    't19', {}, 't20', {}, 't21', {}, 't22', {}, 't23', {}, 't24', {}, 't25', {});


% define event threshold
eventThresh=20; % ... SD Of angles < 20? print out number of times this comes up in a parcel time series to tune
eventThresh=1.5 % SD of signal increases and direction changes

% get number of parcels
parcelTable=tabulate(schaef100_L.cdata);
numParcels=length(parcelTable); % will have to subtract 1 because 0 is considered a parcel in this scheme
numParcels=numParcels-1;
% initialize row counter for struct
RowsCounter=0;
% for each parcel: left
for p=1:numParcels;
	% get parcel indices
	parcInds=find(schaef100_L.cdata==p);
	% pull out signal in this parcel
	parcSig=bold_L(parcInds,:);
	% get signal change TS
	parcSigDelta=diff(parcSig,1,2);
	parcSigDelta_avg=mean(parcSigDelta);
	% pull out vectors in this parcel
	parcVecs=vertAngles_L(parcInds,:,:);
	[num_vectors, num_timepoints, ~] = size(parcVecs);
	% Normalize the vectors
	normVecs = parcVecs ./ vecnorm(parcVecs, 2, 3);
	% Shift the normalized vectors to align consecutive time points
	normVecs1 = normVecs(:, 1:end-1, :); % Vectors at time t
	normVecs2 = normVecs(:, 2:end, :);   % Vectors at time t+1
	% Compute the Euclidean distance between consecutive normalized vectors
	euclideanDist = sqrt(sum((normVecs2 - normVecs1).^2, 3));	
	% get direction change TS
	aggregateChangeOverTime = sum(euclideanDist, 1); % Sum of Euclidean distances for each time point
        %%% get signal/direction change events
	sigChange_SD=std(parcSigDelta_avg);
	vecChange_SD=std(aggregateChangeOverTime);
	% ID > thresh SD signal change TRs
	sigChange_OvThr=find(parcSigDelta_avg>(sigChange_SD*eventThresh));
	% ID > thresh SD vector change TRs
	vecChange_OvThr=find(aggregateChangeOverTime>(vecChange_SD*eventThresh));
	% ID overlap
	overlap_OvThr=intersect(sigChange_OvThr,vecChange_OvThr);
	%%% ID best starting point in continuous sequences
	% within overlap_OvThr, find continuous sequences
	% find the start of each new continuous sequence
	isStartOfSequence = [true, diff(overlap_OvThr) ~= 1];
	% find the end of each continuous sequence
	isEndOfSequence = [diff(overlap_OvThr) ~= 1, true];
	% get the start and end indices of each sequence
	startIndices = find(isStartOfSequence);
	endIndices = find(isEndOfSequence);
	% if there are instances
	if length(overlap_OvThr)>0;
	% extract the continuous sequences using array indexing
	continuousSequences = arrayfun(@(s, e) overlap_OvThr(s:e), startIndices, endIndices, 'UniformOutput', false);
	for c=1:length(continuousSequences)
		continuousSequence=continuousSequences{c};
		if length(continuousSequence>1);
			% get signal within this window
			realIndices=continuousSequence+1;
			realSig=parcSig(realIndices);
			% get peak
			[maximum index]=max(realSig);
			% convert back to difference indices
			differenceIndices=realIndices-1;
			% get nonpeak TRs
			nonPeaks=differenceIndices(differenceIndices~=differenceIndices(index));
			% remove from overlap_OvThr
			overlap_OvThr=overlap_OvThr(~ismember(overlap_OvThr, nonPeaks));
		end
	end
	% convert pack to real indices rather than difference indices
	overlap_OvThr=overlap_OvThr+1;
	% last item is to remove indices <5 to avoid scanner startup artifacts and increase simplicity
	overlap_OvThr=overlap_OvThr(overlap_OvThr>5);
	% should actually apply same criterion to end-of-scan
	overlap_OvThr=overlap_OvThr(overlap_OvThr<(lenOpFl-20));
	% get number of events
	numEvents=length(overlap_OvThr);
	% intialize rows = # of events at end of df
	eventTSs_BOLD=zeros(numEvents,25);
	eventTSs_Vecs=zeros(numEvents,25);
	% for each event
	for E=1:numEvents
		% get indices of interest: 5 TRs before peak and 20 after
		indicesStart=overlap_OvThr(E)-5;
		indicesEnd=overlap_OvThr(E)+20;
		indices=indicesStart:indicesEnd;
		indexCounter=1;
		for i=indices	
			%%% extract synchrony: make distance matrix based on cross-sectional signal
			% fancy way of saying standard deviation
			% signal for this tp
			tpSig=parcSig(:,i);
			tpSynchrony=std(tpSig);
			% insert into eventTSs
			eventTSs_BOLD(E,indexCounter)=tpSynchrony;
			% extract vector uniformity: going to lean into euclidean distance 
			tpVecs=parcVecs(:,i,:);
			tpVecSync=(mean(std(tpVecs)));
			eventTSs_Vecs(E,indexCounter)=tpVecSync;
			% udpate index counter
			indexCounter=indexCounter+1;
		end
	end
	% this is dumb but trying to leverage the struct mixed data benefits
	% Dynamically add new elements to the structure array
	if numEvents>0;
		for i = 1:numEvents;
		    idx = RowsCounter + 1;
		    RowsCounter=RowsCounter+1;
		    outDF(idx).subj = subj;
		    for j = 1:25
		        outDF(idx).(['v' num2str(j)]) = eventTSs_BOLD(i,j);
		        outDF(idx).(['t' num2str(j)]) = eventTSs_Vecs(i,j);  % Initialize t fields as needed
		    end
		    outDF(idx).sesh = sesh;
		    outDF(idx).task = task;
		    outDF(idx).parcel = p;
		    outDF(idx).startTR = overlap_OvThr(i);
		end
	% end if there are no overlapping over thresholds
	end
	% normal end for each parcel
	end
end
% save out to scratch... r-friendly format
save([parentfp '/' task '_' 'Bold_Vec_Sync_L.mat'],'outDF');
% note parcel size is also stored via parcelTable: might need as covariate later
% also note: with removal of 1/std, we are now capturing asynchrony, not synchrony
