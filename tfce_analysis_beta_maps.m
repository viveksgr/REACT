function tfce_analysis_beta_maps(root_dir, nperm, H, E, C, dh)
    % TFCE analysis on subject-wise beta maps for neurotransmitters
    %
    % Inputs:
    %   root_dir - root directory containing subdirectories for each neurotransmitter.
    %   nperm - number of permutations for TFCE
    %   H - TFCE height exponent
    %   E - TFCE extent exponent
    %   C - TFCE connectivity
    %   dh - TFCE step size

    % List all neurotransmitter directories
    neurotransmitter_dirs = dir(fullfile(root_dir, 'react_mask_Normalized*'));
    neurotransmitter_dirs = neurotransmitter_dirs([neurotransmitter_dirs.isdir]);
    
    % Loop through each neurotransmitter
    for n_idx = 1:length(neurotransmitter_dirs)
        neurotransmitter_dir = fullfile(root_dir, neurotransmitter_dirs(n_idx).name);
        fprintf('Processing neurotransmitter: %s\n', neurotransmitter_dirs(n_idx).name);
        
        % Load subject files for nasal and mouth conditions
        nasal_files = dir(fullfile(neurotransmitter_dir, '*_nasal_react_stage2_map1.nii'));
        mouth_files = dir(fullfile(neurotransmitter_dir, '*_mouth_react_stage2_map1.nii'));
        
        if length(nasal_files) ~= length(mouth_files)
            error('Mismatch in the number of nasal and mouth files for %s.', neurotransmitter_dirs(n_idx).name);
        end
        
        num_subjects = length(nasal_files);
        
        % Initialize 4D arrays for beta maps
        stage2_mask_file = fullfile(root_dir,'React_mask',neurotransmitter_dirs(n_idx).name, 'mask_stage2.nii');
        stage2_mask_vol = spm_vol(stage2_mask_file);
        stage2_mask = spm_read_vols(stage2_mask_vol) > 0;

        beta_nasal = zeros([size(stage2_mask), num_subjects]);
        beta_mouth = zeros([size(stage2_mask), num_subjects]);
        
        % Load beta maps for each subject and condition
        for s_idx = 1:num_subjects
            nasal_vol = spm_vol(fullfile(neurotransmitter_dir, nasal_files(s_idx).name));
            mouth_vol = spm_vol(fullfile(neurotransmitter_dir, mouth_files(s_idx).name));
            
            nasal_beta = spm_read_vols(nasal_vol);
            mouth_beta = spm_read_vols(mouth_vol);
            
            % Store masked beta values in 4D arrays
            beta_nasal(:, :, :, s_idx) = nasal_beta;
            beta_mouth(:, :, :, s_idx) = mouth_beta;
        end
        
        % Compute voxel-wise differences (nasal - mouth)
        beta_diff = beta_nasal - beta_mouth;

        % Apply TFCE
        num_subjects = size(beta_diff, 4);
        dummy_covariate = ones(num_subjects, 1);
        [pcorr_pos, pcorr_neg] = matlab_tfce_correlation(beta_diff, dummy_covariate, 2, 5000, 2, 0.5, 26, 0.1);

        % Save corrected p-value maps
        save_statistical_maps(stage2_mask_vol, neurotransmitter_dir, 1 - pcorr_pos, 'tfce_pmap_pos.nii');
        save_statistical_maps(stage2_mask_vol, neurotransmitter_dir, 1 - pcorr_neg, 'tfce_pmap_neg.nii');

    end
end
