% function for finding nearest faces
function [nearest_faces, nearest_distances] = find_nearest_faces(seed_point, faces)
    % Initialize arrays to store nearest faces and distances
    nearest_faces = zeros(1, 3);
    nearest_distances = zeros(1, 3);
    % get coordinate differences
    differences = seed_point - faces;
    % Calculate distances between the seed point and all faces
    distances = sqrt(sum(abs(differences).^2, 2));  % Take the absolute value before squaring
    % Sort distances and get indices of the nearest faces
    [~, sorted_indices] = sort(distances);
    nearest_indices = sorted_indices(1:3);
    % Retrieve the nearest faces and distances
    nearest_faces = nearest_indices;
    nearest_distances = distances(nearest_indices);
end
