function save_nifti(mask_data, ref_vol, filename)
    % Save a binary mask to a NIfTI file using SPM format.
    ref_vol.fname = filename;
    ref_vol.dt = [spm_type('uint8') spm_platform('bigend')];
    spm_write_vol(ref_vol, mask_data);
end
