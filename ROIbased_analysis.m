function [Tmat] = ROIbased_analysis(save_dir,aal_file,aal_labels,opt)
    % Perform group-level analysis on REACT beta files for neurotransmitters.
    %
    % Input:
    %   root_dir - root directory containing subdirectories for each neurotransmitter.

    % opt.num = 90;
    % opt.fname = 'rtmap_pos.nii'
    
    % List all neurotransmitter directories
    neurotransmitter_dirs = dir(fullfile(save_dir, 'react_mask_Normalized*'));
    neurotransmitter_dirs = neurotransmitter_dirs([neurotransmitter_dirs.isdir]);

    % Load AAL data
    aal_masks = spm_read_vols(spm_vol(aal_file));
    aal_labels_list = importdata(aal_labels);
    
    % Initialize results
    Tmat = zeros(length( neurotransmitter_dirs ),opt.num/2);
    
    % Loop through each neurotransmitter
    for n_idx = 1:length(neurotransmitter_dirs)
        neurotransmitter_dir = fullfile(save_dir, neurotransmitter_dirs(n_idx).name);
        fprintf('Processing neurotransmitter: %s\n', neurotransmitter_dirs(n_idx).name);
        
        file_tmap = spm_read_vols(spm_vol(fullfile( neurotransmitter_dir ,opt.fname )));
        Tmat(n_idx,:) = compute_combined_roi_t_scores( file_tmap ,   aal_masks);

        % % Apply TFCE
        % % TFCE parameters (these may vary depending on the library):
        % H = 1; % Height exponent
        % E = 0.1; % Extent exponent
        % dh = 0.1; % Step size
        % connectivity = 26; % 3D connectivity for clusters
        % 
        % % % Apply TFCE
        % file_tmap(isnan(file_tmap))=0;
        % tfce_scores = matlab_tfce_transform(file_tmap , H, E, dh, connectivity);
        % ref_vol = spm_vol(fullfile( neurotransmitter_dir ,opt.fname ));
        % 
        % % save_statistical_maps_neutral(ref_vol, save_dir,  tfce_scores, sprintf('tfce_%s', neurotransmitter_dirs(n_idx).name))
        % 
        % tfce_vol = ref_vol; % Copy metadata
        % tfce_vol.fname = fullfile(save_dir,neurotransmitter_dirs(n_idx).name,'tfce.nii');
        % spm_write_vol(tfce_vol, tfce_scores);
    end
end


