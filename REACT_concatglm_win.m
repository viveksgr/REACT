function REACT_concatglm_win(save_dir, subj_list, pet_atlas_file, mask_file, subjs)
    % Adapted REACT GLM function using SPM for compatibility with Windows.
    % Input:
    %   save_dir: directory for saving output files
    %   subj_list: path to text file listing subject fMRI NIfTI files
    %   pet_atlas_file: path to PET atlas NIfTI file(s)
    %   mask_file: structure containing paths to stage1 and stage2 masks
    %   varargin: additional optional parameters (e.g., nuisance regressors)

    % Parse optional inputs
    opt = struct('mnibrain', '', 'norm_method', 'demean', 'save_name', '', 'nuisance', [], 'param', []);
    % opt = G_SparseArgs(opt, varargin);
    opt.save_name = subjs;
    norm_meth = lower(opt.norm_method);

    % Load subject list
    subj_list = importdata(subj_list);
    tmp = cellfun(@(x) exist(x, 'file') == 2, subj_list);
    if ~all(tmp)
        error('The following fMRI files were not found:\n%s', strjoin(subj_list(~tmp), '\n'));
    end

    % Handle PET atlas files (support for multiple templates)
    if ~iscell(pet_atlas_file)
        all_pet_atlas_file = {pet_atlas_file};
    else
        all_pet_atlas_file = pet_atlas_file;
    end

    % Ensure mask file is a cell array if provided as a struct
    if isstruct(mask_file)
        mask_file = {mask_file};
    end
    if numel(mask_file) ~= numel(all_pet_atlas_file)
        error('The number of PET atlas files must match the number of mask files.');
    end

    % Load MNI brain standard file (e.g., from SPM's templates)
    if isempty(opt.mnibrain)
        % Use SPM’s default templates if not provided
        spm_dir = spm('Dir');
        if ~isempty(spm_dir)
            opt.mnibrain = fullfile(spm_dir, 'canonical', 'avg152T1.nii');
        else
            error('MNI brain file is required.');
        end
    end
    stdbrain = spm_vol(opt.mnibrain);
    stdbrain_data = spm_read_vols(stdbrain);

    % Subject loop
    for sidx = 1:length(subj_list)
        if ~isempty(opt.nuisance)
            nui_mat = opt.nuisance{sidx};
        else
            nui_mat = [];
        end

        in_fmri = subj_list{sidx};
        if isempty(opt.save_name)
            [~, subj_nam, ~] = fileparts(in_fmri);
        else
            subj_nam = opt.save_name{sidx};
        end

        % Load fMRI data
        if endsWith(in_fmri, '.gz')
            in_fmri = decompress_nifti(in_fmri);
        end
        vol_fmri = spm_vol(in_fmri);
        fmri_data = spm_read_vols(vol_fmri);
        [num_voxels, nb_fmri_vol] = size(reshape(fmri_data, [], size(fmri_data, 4)));

        % Load all PET templates and reshape into design matrix
        x_stage1 = [];
        for pet_idx = 1:numel(all_pet_atlas_file)
            pet_file = pet_atlas_file{pet_idx};
            if endsWith(pet_file, '.gz')
                pet_file = decompress_nifti(pet_file);
            end
            pet_data = spm_read_vols(spm_vol(pet_file));
            x_stage1 = [x_stage1, pet_data(stage1_mask)];
        end

        % Normalize fMRI data
        rsfmri = reshape(fmri_data, [], size(fmri_data, 4));
        switch norm_meth
            case 'zscore'
                y_stage1 = (rsfmri(stage1_mask(:), :) - mean(rsfmri(stage1_mask(:), :), 1)) ./ std(rsfmri(stage1_mask(:), :), 0, 1);
            case 'demean'
                y_stage1 = rsfmri(stage1_mask(:), :) - mean(rsfmri(stage1_mask(:), :), 1);
            otherwise
                y_stage1 = rsfmri(stage1_mask(:), :);
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


       %%%%%%%%%%%%%%%%%%%

        % Iterate through each PET atlas file
        beta_mat = zeros(numel(all_pet_atlas_file),nb_fmri_vol);
        for pet_idx = 1:numel(all_pet_atlas_file)
            pet_file = all_pet_atlas_file{pet_idx};
            [~, pet_nam, ~] = fileparts(pet_file);
            pet_save_dir = fullfile(save_dir, ['react_mask_', pet_nam(1:end-4)]);
            if ~exist(pet_save_dir, 'dir')
                mkdir(pet_save_dir);
            end

            % Load masks
            stage1_mask_vol = spm_vol(mask_file{pet_idx}.stage1_mask);
            stage2_mask_vol = spm_vol(mask_file{pet_idx}.stage2_mask);
            stage1_mask = spm_read_vols(stage1_mask_vol) > 0;
            stage2_mask = spm_read_vols(stage2_mask_vol) > 0;

            % Load PET atlas data
            if endsWith(pet_file, '.gz')
                pet_file = decompress_nifti(pet_file);
            end

            pet_vol = spm_vol(pet_file);
            pet_data = spm_read_vols(pet_vol);
            nb_pet_vol = size(pet_data, 4);
            pet_data_reshaped = reshape(pet_data, [], nb_pet_vol);

            % Normalization of fMRI data across time
            rsfmri = reshape(fmri_data, [], nb_fmri_vol);
            switch norm_meth
                case 'zscore'
                    y = (rsfmri - mean(rsfmri, 1)) ./ std(rsfmri, 0, 1);
                case 'demean'
                    y = rsfmri - mean(rsfmri, 1);
                otherwise
                    y = rsfmri;
            end
            y_stage1 = y(stage1_mask(:), :);
            x_stage1 = pet_data_reshaped(stage1_mask(:), :);

            % Stage 1 regression to estimate beta1
            beta1 = nan(nb_fmri_vol, 1);
            for k = 1:nb_fmri_vol
                tmp = regress(y_stage1(:, k), [ones(size(x_stage1, 1), 1), x_stage1]);
                beta1(k) = tmp(2); % Store beta coefficient
            end
            beta_mat(pet_idx,:) = beta1;
            % Save beta1 as a regressor from stage 1
            writematrix(beta1, fullfile(pet_save_dir, [subj_nam, '_react_stage1_beta1.txt']));


            % Stage 2: Normalize and apply beta1 across stage2_mask
            y = y';
            y_stage2 = y(:,stage2_mask(:));
            x_stage2 = (beta1 - mean(beta1)) / std(beta1);

            % Perform regression on each voxel within the stage2 mask
            beta2 = nan(size(y_stage2, 2), nb_pet_vol);
            for k = 1:size(y_stage2, 2)
                tmp = regress(y_stage2(:, k), [ones(size(x_stage2)), x_stage2, nui_mat]);
                beta2(k, :) = tmp(2:size(x_stage2, 2)+1);
            end

            % Save beta2 as maps within the MNI brain template
            for k = 1:nb_pet_vol
                tmp_vol = stdbrain;
                tmp_vol.fname = fullfile(pet_save_dir, sprintf('%s_react_stage2_map%d.nii', subj_nam, k));
                tmp_vol.dt = [spm_type('float32'), spm_platform('bigend')];
                tmp_data = zeros(size(stdbrain_data));
                tmp_data(stage2_mask) = beta2(:, k);
                spm_write_vol(tmp_vol, tmp_data);
            end
        end
    end
end
