function path = geodesic_path(start_idx, end_idx, vertices, faces)
    % Computes the geodesic path along the cortical mesh
    % Inputs:
    %   start_idx - index of the starting vertex
    %   end_idx - index of the ending vertex
    %   vertices - Nx3 matrix of vertex coordinates
    %   faces - Mx3 matrix of triangle indices
    % Output:
    %   path - Nx3 matrix containing geodesic path coordinates

    % Convert faces into an edge list
    edges = [faces(:, [1,2]); faces(:, [2,3]); faces(:, [3,1])];  % Convert triangles into edges

    % Compute edge weights (Euclidean distances)
    edge_weights = vecnorm(vertices(edges(:,1), :) - vertices(edges(:,2), :), 2, 2);

    % Create a weighted graph
    G = graph(edges(:,1), edges(:,2), edge_weights);

    % Compute shortest path using geodesic distances
    path_indices = shortestpath(G, start_idx, end_idx);

    % Extract path coordinates
    path = vertices(path_indices, :);
end

% thanks to matlab devs for shortestpath super useful

