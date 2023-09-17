% function for finding nearest vertices
function nearest_vertex = find_nearest_faces(seed_point, vertices)
    % get coordinate differences
    differences = seed_point - vertices;
    % Calculate distances between the seed point and all faces
    distances = sqrt(sum(differences.^2, 2));
    % Sort distances and get indices of the nearest faces
    [~, sorted_indices] = sort(distances);
    nearest_indices = sorted_indices(1);
    % Retrieve the nearest faces and distances
    nearest_vertex = nearest_indices;
end
