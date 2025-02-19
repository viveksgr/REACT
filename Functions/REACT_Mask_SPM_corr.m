function pet_atlas_vol = REACT_Mask_SPM_corr(pet_atlas_file, gm_mask_file)
% Create masks for multivariate regression using SPM.
% Input:
%   save_dir: directory for saving output
%   subj_list: path to text file listing subject fMRI NIfTI files
%   pet_atlas_file: path to PET atlas file(s) in NIfTI format
%   gm_mask_file: path to gray matter mask file in NIfTI format
%
% Output:
%   out_file: structure containing paths to stage1 and stage2 masks

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

pet_file = pet_atlas_file{1};
if endsWith(pet_file, '.gz')
    pet_file = decompress_nifti(pet_file);
end
pet_vol = spm_read_vols(spm_vol(pet_file));
pet_atlas_vol = zeros(numel( pet_vol), nb_pet_atlas);

for k = 1:nb_pet_atlas
    pet_file = pet_atlas_file{k};
    if endsWith(pet_file, '.gz')
        pet_file = decompress_nifti(pet_file);
    end
    pet_vol = spm_vol(pet_file);

    if any(pet_vol.dim(1:3) ~= gm_dim) || ~isequal(spm_imatrix(pet_vol.mat), gm_resol)
        error('%s has a different dimension or resolution from gray matter mask', pet_atlas_file{k});
    end

    temp_file = spm_read_vols(pet_vol);
    pet_atlas_vol(:,k) =    temp_file (:);
end


end

