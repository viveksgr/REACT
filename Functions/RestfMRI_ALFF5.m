function RestfMRI_ALFF5( idx2run)
% idx2run = 1:17;

save_dir = '/Volumes/iEEG/RestOlfAnaly/Results/fALF_PET';
wkpath = '/Volumes/iEEG/RestOlfAnaly/SubjectData';
subj_list = '/Volumes/iEEG/RestOlfAnaly/nasal_oral_subject_list.txt';

folders = {'Nasal10Min', 'Mouth10Min'};   

% linear transform
file_name = 'filtered_func_data_dspk_WmCSFWb24Motion_globalnorm_s3_rlt2std';
% nonlinear transform
file_name = 'filtered_func_data_dspk_WmCSFWb24Motion_globalnorm_s3_rlt2std_nonlinear';
file_name = 'filtered_func_data_dspk_WmCSFWb24Motion_globalnorm_rlt2std_nonlinear_s3';

% pwelch parameters
pw_nfft = 256;
ft_nfft = 2048;
psd_param = struct( 'pw_nfft', pw_nfft, ...
    'pw_win', hamming( pw_nfft), ...
    'pw_overlap', 0.9, ...
    'ft_nfft', ft_nfft, ...
    'ft_win', hamming( ft_nfft));

f_ois = {[0.008, 0.1], [0.01, 0.1]};
f_ois = {[0.008, 0.07], [0.07, 0.2]};
total_foi = {[], [0, 0.25], [0, 0.2]};

subjs = importdata( subj_list);
save_dir = fullfile( save_dir, file_name);
for subj_idx = 1 : numel( idx2run)
    subj = subjs{ idx2run( subj_idx)};
    for folder_idx = 1 : numel( folders)
        fn = folders{ folder_idx};
        file = fullfile( wkpath, subj, fn, [fn, '.feat'], [file_name, '.nii.gz']);
        if ~isfile( file)
            continue;
        end
        
        surfix = [fn, '_', subj];
        % % Comment the 3 lines below to re-run
        % falff_file = fullfile( save_dir, ['PWSum0-InfHz_FOI0.008-0.1Hz_fALFF_', surfix, '.nii.gz']);
        % if isfile( falff_file)
        %     continue;
        % end

        tr = fMRI_ALFF( file, save_dir, psd_param, f_ois, total_foi, 'surfix', surfix);
    end
end

save( fullfile( save_dir, 'psd_parameters'), 'psd_param', 'tr', 'f_ois', 'total_foi');

end % function


