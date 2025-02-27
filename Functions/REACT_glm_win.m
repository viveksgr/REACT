function REACT_glm_win(save_dir, subj_list, pet_atlas_file, mask_file, subjs,opt)
    % Adapted REACT GLM function using SPM for compatibility with Windows.
    % Input:
    %   save_dir: directory for saving output files
    %   subj_list: path to text file listing subject fMRI NIfTI files
    %   pet_atlas_file: path to PET atlas NIfTI file(s)
    %   mask_file: structure containing paths to stage1 and stage2 masks
    %   varargin: additional optional parameters (e.g., nuisance regressors)

    % Parse optional inputs
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
    fprintf('\n')
    num_pcs = zeros(length(subj_list),numel(all_pet_atlas_file));
    for sidx = 1  :length(subj_list)
        fprintf('Running subject: %02d/%02d\n',sidx,length(subj_list))
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
        
        % nb_fmri_vol2 = 500; 
        nb_fmri_vol2 = nb_fmri_vol;
        % Iterate through each PET atlas file
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

            % nb_fmri_vol = 500;
            rsfmri = rsfmri(:,1:  nb_fmri_vol2);
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
            switch opt.reg_method
                case 'regression'
                    beta1 = nan(nb_fmri_vol2, 1);
                    for k = 1:nb_fmri_vol2
                        tmp = regress(y_stage1(:, k), [ones(size(x_stage1, 1), 1), x_stage1]);
                        beta1(k) = tmp(2); % Store beta coefficient
                        nb_reg = 1;
                    end

                case 'weighted_sum'
                    % Stage 1 regression to estimate beta1
                     beta1 = y_stage1'*x_stage1;
                     nb_reg = 1;

                case 'pca_make'
                    [c,v,~,~,eps_]=pca( y_stage1','VariableWeights', x_stage1');
                    var_exp = cumsum(eps_);
                    num_pcs(sidx,pet_idx) = sum(var_exp<70);
                    beta1 = v(:,1: num_pcs(sidx,pet_idx));
                    nb_reg =   num_pcs(sidx,pet_idx) ;
            end
            % Save beta1 as a regressor from stage 1
            writematrix(beta1, fullfile(pet_save_dir, [subj_nam, '_react_stage1_beta1.txt']));


            % Stage 2: Normalize and apply beta1 across stage2_mask
            y = y';
            y_stage2 = y(:,stage2_mask(:));
            x_stage2 = (beta1 - mean(beta1)) ./ std(beta1);

            % Perform regression on each voxel within the stage2 mask
            beta2 = nan(size(y_stage2, 2), nb_reg);
            for k = 1:size(y_stage2, 2)
                tmp = regress(y_stage2(:, k), [ones(size(x_stage2,1),1), x_stage2, nui_mat]);
                beta2(k, :) = tmp(2:size(x_stage2, 2)+1);
            end

            % Save beta2 as maps within the MNI brain template
            for k = 1:nb_reg
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
