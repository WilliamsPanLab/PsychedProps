function Viz_Streams(mat_fp,out_fp) 
% sort adjacency matrix for each subject
test=load(mat_fp)

% Load the adjacency matrix from the saved image
adjacency_matrix = test.AdjMatrix_L;

% Perform PCA on the adjacency matrix
[coeff, score, latent] = pca(adjacency_matrix);

% Set the diagonal elements to zero using the diag function
adjacency_matrix = adjacency_matrix - diag(diag(adjacency_matrix));

% Sort the vertices by the first principal component
[~, sorted_indices] = sort(score(:, 1));

% Reorder the adjacency matrix based on the sorting
sorted_adjacency_matrix = adjacency_matrix(sorted_indices, sorted_indices);

% Plot the sorted adjacency matrix
imagesc(log(sorted_adjacency_matrix));
colormap('jet');
colorbar;

% Save the sorted adjacency matrix as an image file
saveas(gcf,out_fp)
