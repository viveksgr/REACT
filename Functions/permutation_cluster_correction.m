function permutation_cluster_correction(tmap_path)
    % Load the t-map using SPM
    tmap_vol = spm_vol(tmap_path);
    tmap = spm_read_vols(tmap_vol);
    tmap(isnan(tmap))=0;
    % Define parameters
    voxel_thresh = 2.3; % Uncorrected threshold for individual voxels
    num_permutations = 1000; % Number of permutations
    connectivity = 26; % 3D connectivity for clusters
    
    % Identify significant voxels (positive and negative tails)
    pos_voxels = tmap > voxel_thresh;
    neg_voxels = tmap < -voxel_thresh;
    
    % Find connected clusters
    cc_pos = bwconncomp(pos_voxels, connectivity);
    cc_neg = bwconncomp(neg_voxels, connectivity);
    
    % Calculate observed cluster sizes
    cluster_sizes_pos = cellfun(@numel, cc_pos.PixelIdxList);
    cluster_sizes_neg = cellfun(@numel, cc_neg.PixelIdxList);
    
    % Permutation testing to generate null distribution of cluster sizes
    max_cluster_sizes_pos = zeros(num_permutations, 1);
    max_cluster_sizes_neg = zeros(num_permutations, 1);

    fprintf('\n')
    for p = 1:num_permutations
        if mod(p,10)==0
            fprintf('.')
        end
    end

    fprintf('\n')
    for p = 1:num_permutations
        if mod(p,10)==0
            fprintf('|')
        end
        % Permute the sign of the t-values randomly
        % permuted_tmap = tmap .* (rand(size(tmap)) > 0.5) * 2 - 1;

        % Permutation of voxels
        % Flatten t-map
        tmap_vector = tmap(:);

        % Shuffle t-values randomly across all voxels
        shuffled_tmap_vector = tmap_vector(randperm(numel(tmap_vector)));

        % Reshape back to original image size
        permuted_tmap = reshape(shuffled_tmap_vector, size(tmap));

        % Identify positive and negative clusters separately in permuted map
        perm_pos_voxels = permuted_tmap > voxel_thresh;
        perm_neg_voxels = permuted_tmap < -voxel_thresh;

        % Find connected components separately for positive and negative clusters
        perm_cc_pos = bwconncomp(perm_pos_voxels, connectivity);
        perm_cc_neg = bwconncomp(perm_neg_voxels, connectivity);

        % Compute largest cluster sizes separately
        if perm_cc_pos.NumObjects > 0
            max_cluster_sizes_pos(p) = max(cellfun(@numel, perm_cc_pos.PixelIdxList));
        else
            max_cluster_sizes_pos(p) = 0;
        end

        if perm_cc_neg.NumObjects > 0
            max_cluster_sizes_neg(p) = max(cellfun(@numel, perm_cc_neg.PixelIdxList));
        else
            max_cluster_sizes_neg(p) = 0;
        end
    end
    fprintf('\n')
    % Determine cluster size threshold (95th percentile of permutation distribution)
    cluster_size_thresh_pos = prctile(max_cluster_sizes_pos, 95);
    cluster_size_thresh_neg = prctile(max_cluster_sizes_neg, 95);

    % Identify significant clusters
    significant_clusters_pos = cluster_sizes_pos >= cluster_size_thresh_pos;
    significant_clusters_neg = cluster_sizes_neg >= cluster_size_thresh_neg;

    sum(significant_clusters_neg)
    sum(significant_clusters_pos)
    
    % % Create corrected t-map (retain only significant clusters)
    % corrected_tmap = zeros(size(tmap));
    % 
    % % Keep positive clusters that pass the threshold
    % for i = 1:cc_pos.NumObjects
    %     if significant_clusters_pos(i)
    %         corrected_tmap(cc_pos.PixelIdxList{i}) = tmap(cc_pos.PixelIdxList{i});
    %     end
    % end
    % 
    % % Keep negative clusters that pass the threshold
    % for i = 1:cc_neg.NumObjects
    %     if significant_clusters_neg(i)
    %         corrected_tmap(cc_neg.PixelIdxList{i}) = tmap(cc_neg.PixelIdxList{i});
    %     end
    % end
    % 
    % % Save corrected t-map as a new NIfTI file
    % corrected_vol = tmap_vol; % Copy metadata
    % corrected_vol.fname = fullfile(fileparts(tmap_path), 'corrected_tmap.nii');
    % spm_write_vol(corrected_vol, corrected_tmap);
    % 
    % % Display results
    % fprintf('Permutation-based cluster correction completed.\n');
    % fprintf('Cluster size threshold: %d voxels\n', cluster_size_thresh);
    % fprintf('Corrected t-map saved to: %s\n', corrected_vol.fname);
end
