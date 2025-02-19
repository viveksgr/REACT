function REACT_glm_multi_template(save_dir, subj_list, pet_atlas_file, mask_file, subjs)
    % Adapted REACT GLM function with simultaneous processing of multiple templates
    % using SPM for compatibility with Windows.

    % Parse optional inputs
    opt = struct('mnibrain', '', 'norm_method', 'demean', 'save_name', '', 'nuisance', []);
    opt.save_name = subjs;
    norm_meth = lower(opt.norm_method);

    % Load subject list
    subj_list = importdata(subj_list);
    if any(~cellfun(@(x) exist(x, 'file'), subj_list))
        error('Some fMRI files are missing.');
    end

    % Ensure PET atlas is a cell array
    if ~iscell(pet_atlas_file)
        pet_atlas_file = {pet_atlas_file};
    end
    nb_templates = numel(pet_atlas_file);

    % Ensure mask file is a cell array
    if isstruct(mask_file)
        mask_file = {mask_file};
    end

    % Load MNI brain template
    if isempty(opt.mnibrain)
        spm_dir = spm('Dir');
        opt.mnibrain = fullfile(spm_dir, 'canonical', 'avg152T1.nii');
    end
    stdbrain = spm_vol(opt.mnibrain);
    stdbrain_data = spm_read_vols(stdbrain);

    % Subject loop
    for sidx = 1:length(subj_list)
        subj_file = subj_list{sidx};
        if endsWith(subj_file, '.gz')
            subj_file = decompress_nifti(subj_file);
        end
        if isempty(opt.save_name)
            [~, subj_nam, ~] = fileparts(in_fmri);
        else
            subj_nam = opt.save_name{sidx};
        end


        % Load fMRI data
        vol_fmri = spm_vol(subj_file);
        fmri_data = spm_read_vols(vol_fmri);

        % Load masks
        stage2_mask = spm_read_vols(spm_vol(mask_file{1}.stage2_mask)) > 0;

        % stage1_mask = stage2_mask;
        % Load all PET templates and reshape into design matrix
        x_stage1 = [];
        for pet_idx = 1:nb_templates
             
            pet_file = pet_atlas_file{pet_idx};
            if endsWith(pet_file, '.gz')
                pet_file = decompress_nifti(pet_file);
            end
            pet_data = spm_read_vols(spm_vol(pet_file));
            x_stage1 = [x_stage1, pet_data(stage2_mask)];
        end

        % Normalize fMRI data
        rsfmri = reshape(fmri_data, [], size(fmri_data, 4));
        switch norm_meth
            case 'zscore'
                y_stage1 = (rsfmri(stage2_mask(:), :) - mean(rsfmri(stage2_mask(:), :), 1)) ./ std(rsfmri(stage2_mask(:), :), 0, 1);
            case 'demean'
                y_stage1 = rsfmri(stage2_mask(:), :) - mean(rsfmri(stage2_mask(:), :), 1);
            otherwise
                y_stage1 = rsfmri(stage2_mask(:), :);
        end

        % Stage 1 GLM: Regress spatial beta maps against PET templates for each time step
        x_stage1 = [ones(size(x_stage1, 1), 1), x_stage1]; % Add intercept
        beta_time_series = zeros(size(y_stage1, 2), nb_templates); % Templates x time points

        for time_idx = 1:size(y_stage1, 2) % Iterate over time points
            % Run GLM for spatial beta maps at each time point
            beta = x_stage1 \ y_stage1(:, time_idx); % Regress against PET templates
            beta_time_series(time_idx, :) = beta(2:end)'; % Store coefficients (exclude intercept)
        end

        % Stage 2 GLM: Regress receptor time series (beta_time_series) for each voxel time series
        y_stage2 = rsfmri(stage2_mask(:), :)'; % Time series for each voxel in stage 2 mask
        x_stage2 = beta_time_series; % Receptor time series as regressors

        beta_maps_stage2 = zeros(size(y_stage2, 2), nb_templates);

        for voxel_idx = 1:size(y_stage2, 2) % Iterate over voxels
            % Run GLM for each voxel time series
            tmp = regress(y_stage2(:, voxel_idx), [ones(size(x_stage2, 1), 1), x_stage2]); % Add intercept
            beta_maps_stage2(voxel_idx, :) = tmp(2:end)'; % Store coefficients (exclude intercept)
        end

        % Save beta maps for stage 2
        for pet_idx = 1:nb_templates
             pet_file = pet_atlas_file{pet_idx};
            [~, pet_nam, ~] = fileparts(pet_file);
             pet_save_dir = fullfile(save_dir, ['react_mask_', pet_nam(1:end-4)]);
            if ~exist(pet_save_dir, 'dir')
                mkdir(pet_save_dir);
            end



            tmp_vol = stdbrain;
            tmp_vol.fname = fullfile(pet_save_dir, sprintf('%s_react_stage2_map%d.nii', subj_nam, 1));
              
            % tmp_vol.fname = fullfile(save_dir, sprintf('stage2_map_template%d.nii', pet_idx));
            tmp_vol.dt = [spm_type('float32'), spm_platform('bigend')];
            tmp_data = zeros(size(stdbrain_data));
            tmp_data(stage2_mask) = beta_maps_stage2(:, pet_idx);
            spm_write_vol(tmp_vol, tmp_data);
        end


    end
end
