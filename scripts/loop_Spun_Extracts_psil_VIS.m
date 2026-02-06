% loop over spins to produce equivalent outputs as true measures
numSpins = 2000;
startN = 1;

% Define the output directory and file pattern
outputDir = '/scratch/users/apines/data/';
filePattern = 'rs1_Psil_VISTemporalAutoCor_Merged_Spin_%d.csv';

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
for i = 1:length(missingN)
	n = missingN(i)
	% new networks for revisions
        Extract_ang_dif_psil_Spun_VIS(n);
        Extract_VISMag_dif_psil_Spun(n);
        Extract_VISSeg_dif_psil_Spun(n);
        Extract_TA_dif_psil_Spun_VIS(n);
end
