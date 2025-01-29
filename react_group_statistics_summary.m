function react_group_statistics_summary(root_dir)
    % Group-level statistics on averaged REACT values for each PET template
    %
    % Inputs:
    %   root_dir - the root directory containing subdirectories for each PET template
    
    % Get list of all subdirectories (one per PET template)
    pet_template_dirs = dir(fullfile(root_dir, 'react_mask_Normalized*'));
    pet_template_dirs = pet_template_dirs([pet_template_dirs.isdir]);
    
    % Initialize results storage
    template_results = struct();
    
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
        
        num_subjects = length(condition_n_files);
        
        % Initialize arrays to store mean values for each condition and subject
        mean_values_n = zeros(num_subjects, 1);
        mean_values_m = zeros(num_subjects, 1);
        
        % Load and mask data for each subject
        for s = 1:num_subjects
            % Load and mask subject data for condition "n"
            n_vol = spm_vol(fullfile(pet_dir, condition_n_files(s).name));
            n_data = spm_read_vols(n_vol);
            mean_values_n(s) = mean(n_data(stage2_mask));  % Average value within mask for condition "n"
            
            % Load and mask subject data for condition "m"
            m_vol = spm_vol(fullfile(pet_dir, condition_m_files(s).name));
            m_data = spm_read_vols(m_vol);
            mean_values_m(s) = mean(m_data(stage2_mask));  % Average value within mask for condition "m"
        end
        
        % Perform a t-test on the mean values for this template
        [~, p_value, ~, stats] = ttest(mean_values_n, mean_values_m);
        
        % Store the results for this template
        template_results(d).template_name = pet_template_dirs(d).name;
        template_results(d).t_statistic = stats.tstat;
        template_results(d).p_value = p_value;
        template_results(d).mean_n = mean(mean_values_n);
        template_results(d).mean_m = mean(mean_values_m);
        
        fprintf('Finished processing PET template: %s\n', pet_template_dirs(d).name);
        fprintf('  T-statistic: %.4f, P-value: %.4f\n', stats.tstat, p_value);
    end
    
    % Display final results for all templates
    disp('Group-level t-test results for each PET template:');
    for d = 1:length(template_results)
        fprintf('Template: %s\n', template_results(d).template_name);
        fprintf('  T-statistic: %.4f, P-value: %.4f\n', template_results(d).t_statistic, template_results(d).p_value);
        fprintf('  Mean (Condition n): %.4f, Mean (Condition m): %.4f\n', ...
            template_results(d).mean_n, template_results(d).mean_m);
    end
end
