function out_file = REACT_Mask_SPM(mask_savedir,  subj_list, pet_atlas_file, gm_mask_file,thresh)
    % Create masks for multivariate regression using SPM.
    % Input:
    %   save_dir: directory for saving output
    %   subj_list: path to text file listing subject fMRI NIfTI files
    %   pet_atlas_file: path to PET atlas file(s) in NIfTI format
    %   gm_mask_file: path to gray matter mask file in NIfTI format
    %
    % Output:
    %   out_file: structure containing paths to stage1 and stage2 masks

    % Load gray matter mask and check dimensions
    if nargin<5
        thresh = zeros(1,length( pet_atlas_file))+0.01;
    end

    gray_mask = spm_vol(gm_mask_file);
    gm_dim = gray_mask.dim;
    gm_resol = spm_imatrix(gray_mask.mat);
    gm_orient = spm_get_space(gm_mask_file);

    if length(gm_dim) > 3
        error('Gray matter mask must be a 3D volume.');
    end

    % Load PET atlas files and check dimensions
    if ischar(pet_atlas_file) || isstring(pet_atlas_file)
        pet_atlas_file = {pet_atlas_file};
    end

    nb_pet_atlas = numel(pet_atlas_file);
    pet_atlas_vol = cell(nb_pet_atlas, 1);

    for k = 1:nb_pet_atlas
        pet_file = pet_atlas_file{k};
        if endsWith(pet_file, '.gz')
            pet_file = decompress_nifti(pet_file);
        end
        pet_vol = spm_vol(pet_file);
      
        if any(pet_vol.dim(1:3) ~= gm_dim) || ~isequal(spm_imatrix(pet_vol.mat), gm_resol)
            error('%s has a different dimension or resolution from gray matter mask', pet_atlas_file{k});
        end

        pet_atlas_vol{k} = spm_read_vols(pet_vol);
    end

    % Load subject list and initialize fMRI mask
    subj_list = importdata(subj_list);
    fmri_mask = zeros(gm_dim);
    nb_subj = length(subj_list);

    for k = 1:nb_subj
        fprintf('Reading fMRI %d/%d\n', k, nb_subj);
        subj_file = subj_list{k};
        if endsWith(subj_file, '.gz')
            subj_file = decompress_nifti(subj_file);
        end
        subj_vol = spm_vol(subj_file);

        % if any(subj_vol.dim(1:3) ~= gm_dim) || ~isequal(spm_imatrix(subj_vol.mat), gm_resol)
        %     error('%s has a different dimension or resolution from gray matter mask', subj_list{k});
        % end

        % Read fMRI data
        subj_data = spm_read_vols(subj_vol);
        
        % Create standard deviation mask
        if length(subj_vol) > 1  % 4D data
            subj_data_std = std(subj_data, 0, 4);
        else
            subj_data_std = subj_data;
        end
        ind = subj_data_std > 0;
        
        % Update mask where data has non-zero standard deviation
        fmri_mask(ind) = fmri_mask(ind) + 1;
    end

    % Final fMRI mask (across all subjects)
    fmri_mask = fmri_mask == nb_subj;
    stage2_mask = fmri_mask & spm_read_vols(gray_mask) > 0;

    % Prepare save directory
    if isempty(mask_savedir)
        mask_savedir = fullfile(fileparts(pet_atlas_file{1}), 'React_mask');
    elseif ~exist(mask_savedir, 'dir')
        mkdir(mask_savedir);
    end

    % 
    % % Find threshholds manually
    % figure()
    % hold on
    % k_ind = 0;
    % for k = 1:nb_pet_atlas
    %     k_ind = k_ind+1;
    %     subplot(4,5,k_ind)
    %     vec = pet_atlas_vol{k};
    %     vec = vec(vec>0);
    %     histogram(vec,"EdgeAlpha",0)
    %     title(neurolabels{k})
    % end
    % 

    % Save output masks
    out_file = cell(nb_pet_atlas, 1);
    for k = 1:nb_pet_atlas
        % Intersection of all PET masks
        
        % -- - - - - 

        % Dubious stuff here:
        % pet_mask1 = sum(pet_atlas_vol{k} > 0, 4) == size(pet_atlas_vol{k}, 4);
        pet_mask = pet_atlas_vol{k} > thresh(k);

        % - - - - -

        stage1_mask = stage2_mask & pet_mask;

        % Define file paths for saving
        [~, pet_nam, ~] = fileparts(pet_atlas_file{k}(1:end-4));
        cur_savedir = fullfile(mask_savedir, ['react_mask_', pet_nam]);
        mkdir(cur_savedir)

        % Save stage1 mask
        save_file = fullfile(cur_savedir, 'mask_stage1.nii');
        save_nifti(stage1_mask, gray_mask, save_file);
        out_file{k}.stage1_mask = save_file;

        % Save stage2 mask
        save_file = fullfile(cur_savedir, 'mask_stage2.nii');
        save_nifti(stage2_mask, gray_mask, save_file);
        out_file{k}.stage2_mask = save_file;
    end
end

