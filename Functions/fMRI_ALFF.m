function tr = fMRI_ALFF( fmri_file, save_dir, psd_param, f_ois, total_foi, varargin)
% fmri_file, full-path-to-4d-nifiti
% 
% fractional amplitude of low-frequency oscillations
% 
% file_name = 'filtered_func_data_dspk_WmCSFWb24Motion_globalnorm_s3_rlt2std_nonlinear';
% 
% psd_param = struct( 'pw_nfft', pw_nfft, ...
%     'pw_win', hamming( pw_nfft), ...
%     'pw_overlap', 0.9, ...
%     'ft_nfft', 2048,...
%     'ft_win', hamming( ft_nfft));
% 
% GZ

fsldir = getenv( 'FSLDIR');
if ~isempty( fsldir)
    fsldir = fullfile( fsldir, 'data', 'standard');
    default_std = fullfile( fsldir, 'MNI152_T1_2mm_brain.nii.gz');
    default_stdmask = fullfile( fsldir, 'MNI152_T1_2mm_brain_mask.nii.gz');
else
    default_std = '';
    default_stdmask = '';
end


opt = struct( 'standardbrain', default_std, ...
    'brainmask', default_stdmask,...
    'surfix', '');

opt = G_SparseArgs( opt, varargin);

surfix = sprintf( '%s.nii.gz', opt.surfix);


stdbrain_file = opt.standardbrain;
std_brain = MRIread( stdbrain_file);
std_brain.vol = zeros( size( std_brain.vol));

result =  ( fmri_file, opt.brainmask, 'param', psd_param);
ts_ind = result.vox_ind;
tr = result.tr;

if ~exist( save_dir, 'dir')
    mkdir( save_dir);
end

for fidx = 1 : numel( f_ois)
    f_oi = f_ois{ fidx};

    for j = 1 : numel( total_foi)
        cur_total_foi = total_foi{ j};
        if isempty( cur_total_foi)
            cur_total_foi = [0, Inf];
        end

        prefix = sprintf( 'Sum%sHz_FOI%sHz', ...
            strjoin( string( cur_total_foi), '-'), ...
            strjoin( string( f_oi), '-'));

        if ~isempty( result.pw_pxx)
            cur_prefix = sprintf( 'PW%s', prefix);
            total_foiloc = G_RangeLoc( result.pw_freq, cur_total_foi);
            foiloc = G_RangeLoc( result.pw_freq, f_oi);

            ALFF = sum( result.pw_pxx( foiloc, :), 1);
            total_pow = sum( result.pw_pxx( total_foiloc, :), 1);
            fALFF = sum( result.pw_pxx( foiloc, :), 1) ./ total_pow;

            alff_z = zscore( ALFF);
            falff_z = zscore( fALFF);
            falff_z( ~isfinite( falff_z)) = 0;
            fALFF( ~isfinite( fALFF)) = 0;

            if j == 1
                std_brain.vol( ts_ind) = ALFF;
                sf = fullfile( save_dir, sprintf( '%s_%s_%s', cur_prefix, 'ALFF', surfix));
                MRIwrite( std_brain, sf);
                system( sprintf( 'fslcpgeom %s %s', stdbrain_file, sf));

                std_brain.vol( ts_ind) = alff_z;
                sf = fullfile( save_dir, sprintf( '%s_%s_%s', cur_prefix, 'ALFFZ', surfix));
                MRIwrite( std_brain, sf);
                system( sprintf( 'fslcpgeom %s %s', stdbrain_file, sf));
            end

            std_brain.vol( ts_ind) = fALFF;
            sf = fullfile( save_dir, sprintf( '%s_%s_%s', cur_prefix, 'fALFF', surfix));
            MRIwrite( std_brain, sf);
            system( sprintf( 'fslcpgeom %s %s', stdbrain_file, sf));

            std_brain.vol( ts_ind) = falff_z;
            sf = fullfile( save_dir, sprintf( '%s_%s_%s', cur_prefix, 'fALFFZ', surfix));
            MRIwrite( std_brain, sf);
            system( sprintf( 'fslcpgeom %s %s', stdbrain_file, sf));
        end


        if ~isempty( result.ft_pxx)
            cur_prefix = sprintf( 'FFT%s', prefix);

            total_foiloc = G_RangeLoc( result.ft_freq, cur_total_foi);
            foiloc = G_RangeLoc( result.ft_freq, f_oi);

            ALFF = sum( result.ft_pxx( foiloc, :), 1);
            total_pow = sum( result.ft_pxx( total_foiloc, :), 1);
            fALFF = ALFF ./ total_pow;

            alff_z = zscore( ALFF);
            falff_z = zscore( fALFF);
            falff_z( ~isfinite( falff_z)) = 0;
            fALFF( ~isfinite( fALFF)) = 0;

            if j == 1
                std_brain.vol( ts_ind) = ALFF;
                sf = fullfile( save_dir, sprintf( '%s_%s_%s', cur_prefix, 'ALFF', surfix));
                MRIwrite( std_brain, sf);
                system( sprintf( 'fslcpgeom %s %s', stdbrain_file, sf));

                std_brain.vol( ts_ind) = alff_z;
                sf = fullfile( save_dir, sprintf( '%s_%s_%s', cur_prefix, 'ALFFZ', surfix));
                MRIwrite( std_brain, sf);
                system( sprintf( 'fslcpgeom %s %s', stdbrain_file, sf));
            end

            std_brain.vol( ts_ind) = fALFF;
            sf = fullfile( save_dir, sprintf( '%s_%s_%s', cur_prefix, 'fALFF', surfix));
            MRIwrite( std_brain, sf);
            system( sprintf( 'fslcpgeom %s %s', stdbrain_file, sf));

            std_brain.vol( ts_ind) = falff_z;
            sf = fullfile( save_dir, sprintf( '%s_%s_%s', cur_prefix, 'fALFFZ', surfix));
            MRIwrite( std_brain, sf);
            system( sprintf( 'fslcpgeom %s %s', stdbrain_file, sf));
        end
    end
end


end % function


