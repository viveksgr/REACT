function save_statistical_maps(ref_vol, save_dir, map, fname)
    % Save t-statistic and p-value maps as NIfTI files
    t_stat_vol = ref_vol;
    t_stat_vol.fname = fullfile(save_dir, fname);
    t_stat_data = zeros(size(spm_read_vols(ref_vol)));
    t_stat_data(ref_vol.private.dat(:) > 0) = map;
    spm_write_vol(t_stat_vol, t_stat_data);
end
