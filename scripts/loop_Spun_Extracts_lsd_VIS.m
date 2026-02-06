% loop over spins to produce equivalent outputs as true measures
numSpins = 2000;
startN = 1;

% Define the output directory and file pattern
outputDir = '/scratch/users/apines/data/';
filePattern = 'lsd_propsMerged_Spin_%d_VIS.csv';

% Find which N values are missing
missingN = [];
for n = startN:numSpins
    outputFile = fullfile(outputDir, sprintf(filePattern, n));
    if ~exist(outputFile, 'file')
        missingN = [missingN, n];
    end
end

% Display info about missing files
fprintf('Total spins to process: %d\n', numSpins - startN + 1);
fprintf('Already completed: %d\n', (numSpins - startN + 1) - length(missingN));
fprintf('Missing/To process: %d\n', length(missingN));

% this runs for mdma
for i=1:length(missingN)
	tic
	n = missingN(i)
	% new networks for revisions
        Extract_ang_dif_lsd_Spun_VIS(n);
	toc
        %Extract_VISMag_dif_lsd_Spun(n);
        %Extract_VISSeg_dif_lsd_Spun(n);
        %Extract_TA_dif_lsd_Spun_VIS(n);
end
