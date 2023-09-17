function distances = calculate_distances_to_boundary(network_coordinates, all_coordinates)
    % get number of vertices in this network and broadly
    num_network_vertices = size(network_coordinates, 1);
    num_all_vertices = size(all_coordinates, 1);
    
    % init
    distances = zeros(num_network_vertices, 1);
    
   % for each in-network vertex
    for i = 1:num_network_vertices
        network_vertex = network_coordinates(i, :);
        
        % Calculate distances between the network vertex and all vertices
        differences = all_coordinates - network_vertex;
        distances_to_all_vertices = sqrt(sum(differences.^2, 2));  % Euclidean distance
        
        % Find the closest vertex outside the network
        min_distance = 1000;
        for j = 1:num_all_vertices
            if is_inside_network(all_coordinates(j, :), network_coordinates)
                continue;  % Skip vertices within the network
            end
            
	    % king of the hill at out-of-network vertex distance here, each one gets a chance to be smallest distance
            if distances_to_all_vertices(j) < min_distance
                min_distance = distances_to_all_vertices(j);
            end
        end
        
        distances(i) = min_distance;
    end
end

function inside = is_inside_network(vertex, network_coordinates)
    % Check if a vertex is inside the network
    % You may need to define your criteria for what constitutes being "inside"
    % For example, you could check if the vertex is within a certain radius
    % of the network center or if it belongs to a specific community.
    % Modify this function accordingly.
    
    % Here, we assume that a vertex is inside the network if it is identical
    % to one of the network coordinates.
    inside = ismember(vertex, network_coordinates, 'rows');
end

