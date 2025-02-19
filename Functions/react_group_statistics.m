function react_group_statistics(root_dir, test_type)
    % Group-level statistics on REACT files across PET templates
    %
    % Inputs:
    %   root_dir - the root directory containing subdirectories for each PET template
    %   test_type - 'ttest' for a voxel-wise t-test or 'anova' for an ANOVA analysis
    
    % Get list of all subdirectories (one per PET template)
    pet_template_dirs = dir(fullfile(root_dir, 'react_mask_Normalized*'));
    pet_template_dirs = pet_template_dirs([pet_template_dirs.isdir]);
    
    % Loop through each PET template subdirectory
    for d = 1:length(pet_template_dirs)
        pet_dir = fullfile(root_dir, pet_template_dirs(d).name);
        fprintf('Processing PET template: %s\n', pet_template_dirs(d).name);
        
        % Load stage 2 mask
        stage2_mask_file = fullfile(pet_dir, 'mask_stage2.nii');
        stage2_mask_vol = spm_vol(stage2_mask_file);
        stage2_mask = spm_read_vols(stage2_mask_vol) > 0;

        % Collect subject files for conditions "n" and "m"
        condition_n_files = dir(fullfile(pet_dir, '*_n_react*.nii'));
        condition_m_files = dir(fullfile(pet_dir, '*_m_react*.nii'));
        
        if length(condition_n_files) ~= length(condition_m_files)
            error('Mismatch in number of subject files between conditions "n" and "m".');
        end
        
        % Initialize matrices to store voxel values for each condition
        num_voxels = nnz(stage2_mask);  % Number of voxels in the stage 2 mask
        num_subjects = length(condition_n_files);
        
        % Preallocate arrays for voxel values in masked region for each condition
        condition_n_data = zeros(num_subjects, num_voxels);
        condition_m_data = zeros(num_subjects, num_voxels);
        
        % Load voxel data for each subject under condition "n"
        for s = 1:num_subjects
            % Load and mask subject data for condition "n"
            n_vol = spm_vol(fullfile(pet_dir, condition_n_files(s).name));
            n_data = spm_read_vols(n_vol);
            condition_n_data(s, :) = n_data(stage2_mask);
            
            % Load and mask subject data for condition "m"
            m_vol = spm_vol(fullfile(pet_dir, condition_m_files(s).name));
            m_data = spm_read_vols(m_vol);
            condition_m_data(s, :) = m_data(stage2_mask);
        end
        
        % Perform voxel-wise statistics
        if strcmpi(test_type, 'ttest')
            % T-test between conditions "n" and "m" for each voxel
            [~, p_values, ~, stats] = ttest(condition_n_data, condition_m_data);
            t_values = stats.tstat;
        elseif strcmpi(test_type, 'anova')
            % ANOVA across conditions for each voxel
            all_data = [condition_n_data; condition_m_data];
            group = [ones(num_subjects, 1); 2 * ones(num_subjects, 1)];
            p_values = zeros(1, num_voxels);
            f_values = zeros(1, num_voxels);
            for v = 1:num_voxels
                p = anova1(all_data(:, v), group, 'off');
                p_values(v) = p;
                [~, f_stats] = vartestn(all_data(:, v), group, 'TestType', 'LeveneQuadratic', 'Display', 'off');
                f_values(v) = f_stats.fstat;
            end
        else
            error('Unknown test type. Use "ttest" or "anova".');
        end

        % Prepare the output file for saving statistical results
        output_vol = stage2_mask_vol;
        output_vol.fname = fullfile(pet_dir, sprintf('group_stats_%s.nii', test_type));
        
        % Fill in the results into the masked region of the output volume
        output_data = zeros(size(stage2_mask));
        if strcmpi(test_type, 'ttest')
            output_data(stage2_mask) = t_values;  % Save t-values from t-test
        elseif strcmpi(test_type, 'anova')
            output_data(stage2_mask) = f_values;  % Save F-values from ANOVA
        end
        
        % Write the statistical map to a NIfTI file
        spm_write_vol(output_vol, output_data);
        
        % Save p-values as a separate file
        p_output_vol = stage2_mask_vol;
        p_output_vol.fname = fullfile(pet_dir, sprintf('group_stats_pvals_%s.nii', test_type));
        p_output_data = zeros(size(stage2_mask));
        p_output_data(stage2_mask) = p_values;
        spm_write_vol(p_output_vol, p_output_data);
        
        fprintf('Finished processing PET template: %s\n', pet_template_dirs(d).name);
    end
end
