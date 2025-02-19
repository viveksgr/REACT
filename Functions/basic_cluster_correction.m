function basic_cluster_correction(save_dir,fname)
    % Perform group-level analysis on REACT beta files for neurotransmitters.
    %
    % Input:
    %   root_dir - root directory containing subdirectories for each neurotransmitter.

    % opt.fname = 'rtmap_pos.nii'
    
    % List all neurotransmitter directories
    neurotransmitter_dirs = dir(fullfile(save_dir, 'react_mask_Normalized*'));
    neurotransmitter_dirs = neurotransmitter_dirs([neurotransmitter_dirs.isdir]);
    
    % Loop through each neurotransmitter
    for n_idx = 1:length(neurotransmitter_dirs)
        neurotransmitter_dir = fullfile(save_dir, neurotransmitter_dirs(n_idx).name);
        fprintf('Processing neurotransmitter: %s\n', neurotransmitter_dirs(n_idx).name);
        
        tmap_path = fullfile( neurotransmitter_dir ,fname );       
        basic_cluster_correction_main(tmap_path,2.3)
    end
end


function basic_cluster_correction_main(tmap_path,thr)

    % Load the t-map using SPM
    tmap_vol = spm_vol(tmap_path);
    tmap = spm_read_vols(tmap_vol);
    tmap(isnan(tmap))=0;
    connectivity = 26;
    % Define parameters
    voxel_thresh = 2.3; % Uncorrected threshold for individual voxels
    
    % Identify significant voxels (positive and negative tails)
    pos_voxels = tmap > voxel_thresh;
    neg_voxels = tmap < -voxel_thresh;
    
    % Find connected clusters
    cc_pos = bwconncomp(pos_voxels, connectivity);
    cc_neg = bwconncomp(neg_voxels, connectivity);
    
    % Calculate observed cluster sizes
    cluster_sizes_pos = cellfun(@numel, cc_pos.PixelIdxList);
    cluster_sizes_neg = cellfun(@numel, cc_neg.PixelIdxList);
    
    % Identify significant clusters
    significant_clusters_pos = cluster_sizes_pos >= thr;
    significant_clusters_neg = cluster_sizes_neg >= thr;
    
    % Create corrected t-map (retain only significant clusters)
    corrected_tmap = zeros(size(tmap));

    % Keep positive clusters that pass the threshold
    for i = 1:cc_pos.NumObjects
        if significant_clusters_pos(i)
            corrected_tmap(cc_pos.PixelIdxList{i}) = tmap(cc_pos.PixelIdxList{i});
        end
    end

    % Keep negative clusters that pass the threshold
    for i = 1:cc_neg.NumObjects
        if significant_clusters_neg(i)
            corrected_tmap(cc_neg.PixelIdxList{i}) = tmap(cc_neg.PixelIdxList{i});
        end
    end

    % save_statistical_maps_neutral(tmap_vol, fileparts(tmap_path), corrected_tmap, 'corrected_tmap.nii')
    % Save corrected t-map as a new NIfTI file
    corrected_vol = tmap_vol; % Copy metadata
    corrected_vol.fname = fullfile(fileparts(tmap_path), 'corrected_tmap.nii');
    spm_write_vol(corrected_vol, corrected_tmap);

end
