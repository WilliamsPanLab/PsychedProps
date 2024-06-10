% need to devise approach to aggregate each synchrony struct and print out time series in r friendly format. For a mixed effect model, so we'll need it in long format
% i.e., each row is subject sessions task parcel # and starting TR, as well as 25 timepoints in parcel synchrony and 25 timepoints in vector synchrony
% add paths
ToolFolder='/oak/stanford/groups/leanew1/users/apines/scripts/PersonalCircuits/scripts/code_nmf_cifti/tool_folder';
addpath(genpath(ToolFolder));

% get list of every struct
structList=g_ls(['/scratch/users/apines/data/mdma/*/*Bold_Vec_Sync_L.mat']);

% initialize a table
aggregatedTable = table();

lengthStructs=length(structList);

% load in accordingly
for s=1:lengthStructs
	s
	% load struct
	loadedStruct=load(structList{s}).outDF;
	% get size of struct to add to a table
	lengthStruct=size(loadedStruct);
	lengthStruct=lengthStruct(2);

	% iterate over each entry in the struct
	for i = 1:length(loadedStruct)
        	% Pull relevant fields
        	subj = loadedStruct(i).subj;
        	sesh = loadedStruct(i).sesh;
        	task = loadedStruct(i).task;
        	parcel = loadedStruct(i).parcel;
        	startTR = loadedStruct(i).startTR;
        
        	% Pull BOLD time series
        	boldTimeSeries = arrayfun(@(x) loadedStruct(i).(['t' num2str(x)]), 1:25, 'UniformOutput', false);
        	boldTimeSeries = cell2mat(boldTimeSeries);  % Convert cell array to a matrix
        
        	% Pull vector time series
        	vectorTimeSeries = arrayfun(@(x) loadedStruct(i).(['v' num2str(x)]), 1:25, 'UniformOutput', false);
        	vectorTimeSeries = cell2mat(vectorTimeSeries);  % Convert cell array to a matrix
        
        	% Create a temporary table for the current entry
       		tempTable = table({subj}, {sesh}, {task}, {parcel}, startTR, ...
                          boldTimeSeries, vectorTimeSeries, ...
                          'VariableNames', {'Subject', 'Session', 'Task', 'Parcel', 'StartTR', 'BoldSynchrony', 'VectorSynchrony'});
 
       		% Append the temporary table to the aggregated table
        	aggregatedTable = [aggregatedTable; tempTable];
    	end
end
% Convert cell columns to regular columns for CSV export
aggregatedTable.Subject = categorical(aggregatedTable.Subject);
aggregatedTable.Session = categorical(aggregatedTable.Session);
aggregatedTable.Task = categorical(aggregatedTable.Task);
aggregatedTable.Parcel = cell2mat(aggregatedTable.Parcel);

% Save the aggregated table to a CSV file
outputFileName = '/scratch/users/apines/data/mdma/aggregated_IntraParcel_ot_L.csv';
writetable(aggregatedTable, outputFileName);
