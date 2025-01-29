function react_group_statistics_anova(root_dir)
    % Group-level ANOVA on REACT values across PET templates without averaging voxels
    %
    % Input:
    %   root_dir - the root directory containing subdirectories for each PET template
    
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
        
        num_subjects = length(condition_n_files);
        num_voxels = nnz(stage2_mask);  % Number of voxels in the stage 2 mask
        
        % Preallocate arrays to store data and factors for ANOVA
        voxel_data = [];
        condition_factors = [];
        subject_factors = [];
        voxel_factors = [];
        
        % Load data for each subject and each condition, and reshape
        for s = 1:num_subjects
            % Load and mask subject data for condition "n"
            n_vol = spm_vol(fullfile(pet_dir, condition_n_files(s).name));
            n_data = spm_read_vols(n_vol);
            n_data_masked = n_data(stage2_mask);
            
            % Load and mask subject data for condition "m"
            m_vol = spm_vol(fullfile(pet_dir, condition_m_files(s).name));
            m_data = spm_read_vols(m_vol);
            m_data_masked = m_data(stage2_mask);
            
            % Append data to voxel_data array
            voxel_data = [voxel_data; n_data_masked; m_data_masked];
            
            % Define factors for each entry
            condition_factors = [condition_factors; repmat({'n'}, num_voxels, 1); repmat({'m'}, num_voxels, 1)];
            subject_factors = [subject_factors; repmat(s, num_voxels * 2, 1)];
            voxel_factors = [voxel_factors; (1:num_voxels)'; (1:num_voxels)'];
        end
        
        % Create a table for mixed-effects model
        tbl = table(voxel_data, condition_factors, subject_factors, voxel_factors, ...
            'VariableNames', {'VoxelData', 'Condition', 'Subject', 'Voxel'});
        
        % Fit mixed-effects model treating "Subject" as a random effect
        lme = fitlme(tbl, 'VoxelData ~ Condition + (1|Subject) + (1|Voxel)', ...
            'DummyVarCoding', 'effects');
        
        % Get ANOVA results
        anova_results = anova(lme, 'DFMethod', 'Satterthwaite');
        
        % Display results for this PET template
        fprintf('ANOVA results for PET template: %s\n', pet_template_dirs(d).name);
        disp(anova_results);
        
        % Extract p-value and F-statistic for the Condition effect
        condition_row = anova_results(strcmp(anova_results.Term, 'Condition'), :);
        F_stat = condition_row.FStat;
        p_value = condition_row.pValue;
        
        fprintf('Condition Effect - F-statistic: %.4f, P-value: %.4f\n', F_stat, p_value);
    end
end
