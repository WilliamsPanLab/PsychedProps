% loop over spins to produce equivalent outputs as true measures
numSpins = 2000;
startN = 1;

% Define the output directory and file pattern
outputDir = '/scratch/users/apines/data/';
filePattern = 'rs1_Psil_MOTTemporalAutoCor_Merged_Spin_%d.csv';

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
	n = missingN(i)
        Extract_ang_dif_psil_Spun_MOT(n);
        Extract_MOTMag_dif_psil_Spun(n);
        Extract_MOTSeg_dif_psil_Spun(n);
	Extract_TA_dif_psil_Spun_MOT(n);
end
