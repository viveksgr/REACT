function [similarity_matrix, sorted_indices] = cluster_templates_by_normalized_significance(t_matrix, threshold, num_clusters)
    % Cluster templates based on normalized voxel significance (positive/negative).
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
    positive_binary = t_matrix > threshold;   % Positive significance
    negative_binary = t_matrix < -threshold; % Negative significance

    % Step 2: Compute pairwise normalized signed similarity
    num_templates = size(t_matrix, 2);
    similarity_matrix = zeros(num_templates, num_templates);

    for i = 1:num_templates
        for j = i+1:num_templates
            % Compute total significant voxels for both templates
            total_significant_voxels = sum(positive_binary(:, i) | positive_binary(:, j)) + ...
                                       sum(negative_binary(:, i) | negative_binary(:, j));

            % Positive overlap (normalized by total positive voxels)
            pos_intersection = sum(positive_binary(:, i) & positive_binary(:, j));
            pos_similarity = pos_intersection / (total_significant_voxels + eps);

            % Negative overlap (normalized by total negative voxels)
            neg_intersection = sum(negative_binary(:, i) & negative_binary(:, j));
            neg_similarity = neg_intersection / (total_significant_voxels + eps);

            % Penalize mismatch between positive and negative significance
            pos_neg_mismatch = sum((positive_binary(:, i) & negative_binary(:, j)) | ...
                                    (negative_binary(:, i) & positive_binary(:, j)));
            mismatch_penalty = pos_neg_mismatch / (total_significant_voxels + eps);

            % Combined similarity
            similarity_matrix(i, j) = pos_similarity + neg_similarity - mismatch_penalty;
            similarity_matrix(j, i) = similarity_matrix(i, j); % Symmetry
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
    imagesc(similarity_matrix(sorted_indices, sorted_indices));
    colorbar;
    title('Normalized Similarity Matrix Sorted by Clustering');
    xlabel('Templates');
    ylabel('Templates');
end
