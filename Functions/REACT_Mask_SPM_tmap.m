function pet_atlas_vol = REACT_Mask_SPM_tmap(pet_atlas_file)

dirs = dir(fullfile(pet_atlas_file,'react_mask_Normalized*'));
nb_pet_atlas = numel(dirs);

pet_file = dirs(1).name;
pet_vol = spm_read_vols(spm_vol(fullfile(pet_atlas_file,pet_file,'tmap_pos.nii')));
pet_atlas_vol = zeros(numel( pet_vol), nb_pet_atlas);

for k = 1:nb_pet_atlas
    pet_file = dirs(k).name;

    pet_vol = spm_read_vols(spm_vol(fullfile(pet_atlas_file,pet_file,'tmap_pos.nii')));
    pet_atlas_vol(:,k) =  pet_vol(:);
end


end

