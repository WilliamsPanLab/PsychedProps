
% loop over spins to produce equivalent outputs as true measures
numSpins = 2000;
startN = 1;

% Define the output directory and file pattern
outputDir = '/scratch/users/apines/data/';
filePattern = 'rs1_TAutoCorMerged_Spin_%d_FPN.csv';

% Find which N values are missing
%missingN = [];
%for n = startN:numSpins
%    outputFile = fullfile(outputDir, sprintf(filePattern, n));
%    if ~exist(outputFile, 'file')
%        missingN = [missingN, n];
%    end
%end

% Display info about missing files
%fprintf('Total spins to process: %d\n', numSpins - startN + 1);
%fprintf('Already completed: %d\n', (numSpins - startN + 1) - length(missingN));
%fprintf('Missing/To process: %d\n', length(missingN));

% this runs for mdma
%for i=1:length(missingN)
%	n = missingN(i)
for n=startN:numSpins
	Extract_ang_dif_Spun_FPN(n);
        Extract_FPNMag_dif_Spun(n);
        Extract_FPNSeg_dif_Spun(n);
        Extract_TAutocor_dif_Spun_FPN(n);
end
