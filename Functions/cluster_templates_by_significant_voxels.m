function [similarity_matrix, sorted_indices] = cluster_templates_by_significant_voxels(t_matrix, threshold, num_clusters)
    % Cluster templates based on significant voxel overlap.
    %
    % Input:
    %   t_matrix: A matrix of size (voxels x templates) with t-scores
    %   threshold: T-score threshold for significance (e.g., 3)
    %   num_clusters: Number of clusters to form (e.g., 4)
    %
    % Output:
    %   similarity_matrix: Pairwise similarity matrix of templates
    %   sorted_indices: Indices of templates sorted by clustering

    % Step 1: Threshold the T-score matrix
    binary_matrix = abs(t_matrix) > threshold; % Binary matrix (1 for significant, 0 otherwise)

    % Step 2: Compute pairwise Jaccard similarity
    num_templates = size(binary_matrix, 2);
    similarity_matrix = zeros(num_templates, num_templates);

    for i = 1:num_templates
        for j = i:num_templates
            % Compute Jaccard similarity
            intersection = sum(binary_matrix(:, i) & binary_matrix(:, j));
            union = sum(binary_matrix(:, i) | binary_matrix(:, j));
            similarity_matrix(i, j) = intersection / union;
            similarity_matrix(j, i) = similarity_matrix(i, j); % Symmetric
        end
    end

    % Step 3: Convert similarity to distance
    distance_matrix = 1 - similarity_matrix;

    % Step 4: Cluster templates
    % Perform hierarchical clustering (Ward's method)
    Z = linkage(distance_matrix, 'ward');
    cluster_labels = cluster(Z, 'maxclust', num_clusters);

    % Step 5: Sort templates by cluster
    [~, sorted_indices] = sort(cluster_labels);

    % Display results (optional visualization)
    figure;
    img = similarity_matrix(sorted_indices, sorted_indices);
    img(logical(eye(size(img))))=nan;
    imagesc(img);
    colorbar;
    title('Similarity Matrix Sorted by Clustering');
    xlabel('Templates');
    ylabel('Templates');
end
