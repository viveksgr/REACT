function REACT_glm_constant(save_dir, subj_list, subjs,stage2maskpath)
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
subj_list = cellfun(@(x) x(1:end-3),subj_list,'UniformOutput',false);


tmp = cellfun(@(x) exist(x, 'file') == 2, subj_list);
if ~all(tmp)
    error('The following fMRI files were not found:\n%s', strjoin(subj_list(~tmp), '\n'));
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


% Stage 2 mask
stage2_mask_vol = spm_vol(stage2maskpath);
stage2_mask = spm_read_vols(stage2_mask_vol) > 0;

pet_save_dir = fullfile(save_dir, ['react_mask_', 'Normalized_Const']);
if ~exist(pet_save_dir, 'dir')
    mkdir(pet_save_dir);
end
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
    [~, nb_fmri_vol] = size(reshape(fmri_data, [], size(fmri_data, 4)));

    % Normalization of fMRI data across time
    rsfmri = reshape(fmri_data, [], nb_fmri_vol);

    rsfmri = rsfmri(:,1:  nb_fmri_vol);
    switch norm_meth
        case 'zscore'
            y = ((rsfmri - mean(rsfmri, 1)) ./ std(rsfmri, 0, 1))';
        case 'demean'
            y = (rsfmri - mean(rsfmri, 1))';
        otherwise
            y = rsfmri';
    end

    % Stage 2: Normalize and apply beta1 across stage2_mask
    y_stage2 = y(:,stage2_mask(:));
    beta2 = mean(y_stage2)';


    tmp_vol = stdbrain;
    tmp_vol.fname = fullfile(pet_save_dir, sprintf('%s_react_stage2_map1.nii', subj_nam));
    tmp_vol.dt = [spm_type('float32'), spm_platform('bigend')];
    tmp_data = zeros(size(stdbrain_data));
    tmp_data(stage2_mask) = beta2;
    spm_write_vol(tmp_vol, tmp_data);
end
end
