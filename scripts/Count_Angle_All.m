% define the three study directories
study1_dir = '/scratch/users/apines/LSD_ICL/rest_proc/';
study2_dir = '/scratch/users/apines/data/psil/';
study3_dir = '/scratch/users/apines/data/mdma/';

% Define histogram bins (0-10, 10-20, ..., 170-180)
bin_edges = 0:10:180;
bin_counts = zeros(1, length(bin_edges)-1);

fprintf('Starting to process files...\n\n');

%% Process Study 1: LSD_ICL
fprintf('Processing Study 1 (LSD_ICL)...\n');
study1_files = dir(fullfile(study1_dir, '**', '*k1_Prop_TS_dmn*.csv'));
fprintf('Found %d files in Study 1\n', length(study1_files));

for i = 1:length(study1_files)
    filepath = fullfile(study1_files(i).folder, study1_files(i).name);
    data = csvread(filepath);
    data = data(:); % Flatten to vector
    bin_counts = bin_counts + histcounts(data, bin_edges);
    if mod(i, 50) == 0
        fprintf('  Processed %d/%d files\n', i, length(study1_files));
    end
end
fprintf('Study 1 complete.\n\n');

%% Process Study 2: Psilocybin
fprintf('Processing Study 2 (Psilocybin)...\n');
study2_files = dir(fullfile(study2_dir, '**', '*k1_Prop_TS_dmn*.csv'));
fprintf('Found %d files in Study 2\n', length(study2_files));

for i = 1:length(study2_files)
    filepath = fullfile(study2_files(i).folder, study2_files(i).name);
    data = csvread(filepath);
    data = data(:); % Flatten to vector
    bin_counts = bin_counts + histcounts(data, bin_edges);
    if mod(i, 50) == 0
        fprintf('  Processed %d/%d files\n', i, length(study2_files));
    end
end
fprintf('Study 2 complete.\n\n');

%% Process Study 3: MDMA
fprintf('Processing Study 3 (MDMA)...\n');
study3_files = dir(fullfile(study3_dir, '**', '*k1_Prop_TS_dmn*.csv'));
fprintf('Found %d files in Study 3\n', length(study3_files));

for i = 1:length(study3_files)
    filepath = fullfile(study3_files(i).folder, study3_files(i).name);
    data = csvread(filepath);
    data = data(:); % Flatten to vector
    bin_counts = bin_counts + histcounts(data, bin_edges);
    if mod(i, 50) == 0
        fprintf('  Processed %d/%d files\n', i, length(study3_files));
    end
end
fprintf('Study 3 complete.\n\n');

%% Create output table
% Create bin labels (e.g., "0-10", "10-20", etc.)
bin_labels = cell(length(bin_edges)-1, 1);
for i = 1:length(bin_edges)-1
    bin_labels{i} = sprintf('%d-%d', bin_edges(i), bin_edges(i+1));
end

% Create table with bin labels and counts
output_table = table(bin_labels, bin_counts', 'VariableNames', {'Bin', 'Count'});

% Display summary
fprintf('=== SUMMARY ===\n');
fprintf('Total files processed: %d\n', length(study1_files) + length(study2_files) + length(study3_files));
fprintf('Total values aggregated: %d\n', sum(bin_counts));
fprintf('\nHistogram counts per bin:\n');
disp(output_table);

%% Save to CSV
output_file = 'dmn_histogram_counts.csv';
writetable(output_table, output_file);
fprintf('\nResults saved to: %s\n', output_file);

%% Optional: Create visualization
figure;
bar(1:length(bin_counts), bin_counts);
set(gca, 'XTick', 1:length(bin_counts), 'XTickLabel', bin_labels, 'XTickLabelRotation', 45);
xlabel('Angle (degrees)');
ylabel('Count across all scans');
title('Optical flow vectors relative the the gradient of the DMN');
ylim([1e7, max(bin_counts)*1.01]);
grid on;
saveas(gcf, '~/dmn_histogram.png');
fprintf('Histogram plot saved to: dmn_histogram.png\n');
