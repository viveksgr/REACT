function save_statistical_maps_neutral(ref_vol, save_dir, map, fname)
    % Save continuous maps as NIfTI files
    %
    % Inputs:
    %   ref_vol: Reference SPM volume structure (e.g., from spm_vol)
    %   save_dir: Directory to save the output NIfTI file
    %   map: Continuous variable to write to the NIfTI file
    %   fname: Output filename

    % Create a new volume structure to avoid modifying the original
    new_vol = ref_vol;

    % Update data type to float32 for continuous data
    new_vol.dt = [spm_type('float32'), spm_platform('bigend')];

    % Prepare output volume
    new_vol.fname = fullfile(save_dir, fname);
    t_stat_data = zeros(size(spm_read_vols(ref_vol)));

    % Assign map values to valid voxels
    valid_voxels = ref_vol.private.dat(:) > 0;
    t_stat_data(valid_voxels) = map;

    % Write the new volume
    spm_write_vol(new_vol, t_stat_data);
end
