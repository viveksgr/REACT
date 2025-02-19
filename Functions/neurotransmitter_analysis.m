function [results, groupstats,corrvox] = neurotransmitter_analysis(root_dir)
    % Perform group-level analysis on REACT beta files for neurotransmitters.
    %
    % Input:
    %   root_dir - root directory containing subdirectories for each neurotransmitter.
    
    % List all neurotransmitter directories
    neurotransmitter_dirs = dir(fullfile(root_dir, 'react_mask_Normalized*'));
    neurotransmitter_dirs = neurotransmitter_dirs([neurotransmitter_dirs.isdir]);
    
    % Initialize results
    results = struct();
    groupstats = struct();

    map_score = zeros(1,length(neurotransmitter_dirs));
    map_labels = cell(1,length(neurotransmitter_dirs));
    map_pval = zeros(1,length(neurotransmitter_dirs));
    stage2_mask_file = fullfile(root_dir,'React_mask',neurotransmitter_dirs(1).name, 'mask_stage2.nii');
    stage2_mask_vol = spm_vol(stage2_mask_file);
    stage2_mask = spm_read_vols(stage2_mask_vol) > 0;
    num_voxels = nnz(stage2_mask);

    corrvox.m = zeros(num_voxels,length(neurotransmitter_dirs));
    corrvox.n = zeros(num_voxels,length(neurotransmitter_dirs));
    corrvox.p = zeros(num_voxels,length(neurotransmitter_dirs));
    
    % Loop through each neurotransmitter
    for n_idx = 1:length(neurotransmitter_dirs)
        neurotransmitter_dir = fullfile(root_dir, neurotransmitter_dirs(n_idx).name);
        fprintf('Processing neurotransmitter: %s\n', neurotransmitter_dirs(n_idx).name);
        
        % Load subject files for nasal and mouth conditions
        nasal_files = dir(fullfile(neurotransmitter_dir, '*_nasal_react_stage2_map*.nii'));
        filtered_nasal_files = filter_image_files( nasal_files,'nasal',13);
        mouth_files = dir(fullfile(neurotransmitter_dir, '*_mouth_react_stage2_map*.nii'));
        filtered_mouth_files = filter_image_files( mouth_files,'mouth',13);

        if length(filtered_nasal_files) ~= length(filtered_mouth_files)
            error('Mismatch in the number of nasal and mouth files for %s.', neurotransmitter_dirs(n_idx).name);
        end
        
        num_subjects = length(mouth_files);
        
        % Load stage 2 mask (assumes it is common for all subjects and conditions)
        stage2_mask_file = fullfile(root_dir,'React_mask',neurotransmitter_dirs(n_idx).name, 'mask_stage2.nii');
        stage2_mask_vol = spm_vol(stage2_mask_file);
        stage2_mask = spm_read_vols(stage2_mask_vol) > 0;
        num_voxels = nnz(stage2_mask);
        
        % Initialize data matrices
        nasal_data = zeros(num_subjects, num_voxels);
        mouth_data = zeros(num_subjects, num_voxels);
        
        % Extract beta values for each subject
        for s_idx = 1:num_subjects
            % Load and mask nasal condition
            nasal_vol = spm_vol(fullfile(neurotransmitter_dir, nasal_files(s_idx).name));
            nasal_beta = spm_read_vols(nasal_vol);
            nasal_data(s_idx, :) = nasal_beta(stage2_mask);
            
            % Load and mask mouth condition
            mouth_vol = spm_vol(fullfile(neurotransmitter_dir, mouth_files(s_idx).name));
            mouth_beta = spm_read_vols(mouth_vol);
            mouth_data(s_idx, :) = mouth_beta(stage2_mask);
        end

        corrvox.m(:,n_idx)=mean(mouth_data)';
        corrvox.n(:,n_idx)=mean(nasal_data)';
        corrvox.p(:,n_idx)=mean(mouth_data - nasal_data)';
        
        % Perform statistical test (e.g., paired t-test)
        [~, p_values, ~, stats] = ttest(nasal_data, mouth_data);
        corrvox.p(:,n_idx)=stats.tstat';
        
        % Store results
        results(n_idx).neurotransmitter = neurotransmitter_dirs(n_idx).name;
        results(n_idx).t_statistic = stats.tstat; % T-statistic map      
        results(n_idx).p_values = p_values;      % P-value map
        

        [p_thr,p_values_adj] = fdr_benjhoc(p_values);

        % Save statistical maps
        save_statistical_maps_neutral(stage2_mask_vol, neurotransmitter_dir, stats.tstat, 'tmap_pos.nii')
        % save_statistical_maps(stage2_mask_vol, neurotransmitter_dir, stats.tstat, 'tmap.nii')
        save_statistical_maps_neutral(stage2_mask_vol, neurotransmitter_dir, -stats.tstat, 'tmap_neg.nii')
        save_statistical_maps_neutral(stage2_mask_vol, neurotransmitter_dir, 1-p_values_adj, 'pmap_neg.nii')
        % save_statistical_maps(stage2_mask_vol, stats.tstat, p_values_adj, neurotransmitter_dir);
        
        % fprintf('Finished processing %s\n', neurotransmitter_dirs(n_idx).name);
         map_score(n_idx) = sum(p_values_adj<0.05)/length(p_values_adj);
         % map_score(n_idx) = sum(p_values<0.05)/length(p_values);
         map_labels{n_idx}=neurotransmitter_dirs(n_idx).name;
         map_pval(n_idx)=tinv(p_thr,length(p_values_adj));
    end
    groupstats.map_score = map_score;
    groupstats.map_labels = map_labels;
    groupstats.map_pval = map_pval;
    % % Display summary of results
    % disp('Summary of group-level analysis:');
    % for n_idx = 1:length(results)
    %     fprintf('Neurotransmitter: %s\n', results(n_idx).neurotransmitter);
    % end
end

