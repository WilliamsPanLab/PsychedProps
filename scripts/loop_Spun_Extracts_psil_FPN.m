% loop over spins to produce equivalent outputs as true measures
numSpins = 2000;
startN = 1;

% Define the output directory and file pattern
outputDir = '/oak/stanford/groups/leanew1/users/apines/data/';
filePattern = 'rs1_Psil_FPNTemporalAutoCor_Merged_Spin_%d.csv';

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

% Only run for missing N values
for i = 1:length(missingN)
    n = missingN(i);
    fprintf('Processing spin %d (%d of %d remaining)...\n', n, i, length(missingN));    
    Extract_ang_dif_psil_Spun_FPN(n);
    Extract_FPNMag_dif_psil_Spun(n);
    Extract_FPNSeg_dif_psil_Spun(n);
    Extract_TA_dif_psil_Spun_FPN(n);
end
