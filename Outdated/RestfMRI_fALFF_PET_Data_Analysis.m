%% RestfMRI_fALFF_PET_Data_Analysis

%% PET atlas pre-processing
% % % PET atlas can be downloaded from https://github.com/netneurolab/hansen_receptors
% % 
% Python for REACT anaylysis: https://github.com/ottaviadipasquale/react-fmri
% % 
% 
% # 5HT1a 5HT1b 5HT2a 5HT4 5HT6 5HTT A4B2 CB1 D1 D2 DAT FDOPA GABAa H3 M1 MU NAT NMDA VAChT mGluR5
% # Serotonin: 5HT1a 5HT1b 5HT2a 5HT4 5HT6 5HTT
% # Dopamine: D1 D2 DAT
% # Norepinephrine: NET? or NAT
% # Histamine: H3
% # Acetylcholine: A4B2 M1 VAChT
% # Cannabinoid: CB1
% # Opioid: MOR? MU
% # Glutamate: NMDA mGluR5
% # GABA: GABA A/BZ
% # FDOPA MU MOR? NAT NET?



% 1. Downsample raw images to 2 mm resolution
% 2. Merge all images of the same receptor into a 4d volume
% 3. Normalize each 3d image and average across all images for each receptor
% See REACT_Combine_PET.sh for steps 1&2


% % clear; clc;
% %
% % % images were resampled to 2 mm and combined (into 4d) for each receptor
% % pet_dir = '~/Desktop/test/REACT/PET_Atlas/CombineByReceptor';
% %
% % d = fileparts( pet_dir);
% % save_dir = G_Fullfile( d, 'CombineByReceptor_Normalized');
% % files = G_UFind( pet_dir, 'name', 'Combined*.nii.gz');
% % ext = '.nii.gz';
% % for k = 1 : length( files)
% %     [~, nam] = G_Fileparts( files{ k}, ext);
% %     nam = extractAfter( nam, '_');
% %     nam = sprintf( 'Normalized_%s_atlas%s', nam, ext);
% %     out_file = fullfile( save_dir, nam);
% %     if isfile( out_file)
% %         error( '%s already exists.', out_file);
% %     end
% %
% %     REACT_Normalize( files{ k}, 'savefile', out_file, 'average', 'yes', 'intersect', 'yes');
% % end


%% fMRI data preprocessing
% Use FSL's FEAT_GUI for motion correction
% RestfMRI_ALFF_Preprocessing.sh

% Calculate fALFF
% RestfMRI_ALFF5( 1:17)

% Combine fALFF maps across subjects and compare between conditions.
% RestfMRI_ALFF_Matlab_Stats.sh


%% Compare oral to nasal including interaction in GLM
% hippocampal memory related receptors
% # Acetylcholine: A4B2 M1 VAChT
% # Glutamate: NMDA mGluR5

clear; clc; close all

wkpath = ['E:\RestOlfAnaly\Results\fALF_PET\' ...
    'filtered_func_data_dspk_WmCSFWb24Motion_globalnorm_rlt2std_nonlinear_s3/stats/Mouth_Nasal'];

% file_names = {'merged_PWSum0-0.25Hz_FOI0.008-0.1Hz_fALFFZ_Mouth_Nasal'
%     'merged_PWSum0-InfHz_FOI0.008-0.1Hz_fALFFZ_Mouth_Nasal'};

file_names = {'merged_PWSum0-0.25Hz_FOI0.01-0.1Hz_fALFFZ_Mouth_Nasal'
    'merged_PWSum0-InfHz_FOI0.01-0.1Hz_fALFFZ_Mouth_Nasal'};

gray_mask = '~/Desktop/Data/MNIROIs/avg152T1_gray_thresh100.nii.gz';
pet_atlas_dir = '~/Desktop/test/REACT/PET_Atlas/CombineByReceptor_Normalized';
load( fullfile( pet_atlas_dir, 'MRIread_CombineByReceptor_Normalized'));


fig_dir = fullfile( wkpath, 'Figures');

pet_ind = [1, 10, 12];
pet_ind = [7, 20];
% pet_ind = [9, 14, 15];
targ_pet_ind = {[1, 10, 12], [7 20]};
gray_mask_vol = MRIread( gray_mask);


for file_idx = 1 : numel( file_names)
    pval = [];
    tval = [];
    result_txt = {};

    for ind_idx = 1 : numel( targ_pet_ind)
        pet_ind = targ_pet_ind{ ind_idx};

        pet_receptor_names = pet_atlas_name( pet_ind);

        file_name = file_names{ file_idx};
        fmri_file = fullfile( wkpath, file_name);
        [~, ts, gray_mask_ind] = G_ExtractTS( fmri_file, gray_mask_vol);

        mr = MRIread( fullfile( wkpath, [file_name, '.nii.gz']));
        nbsubjs = floor( size( mr.vol, 4)/2);
        bval = nan( 2, nbsubjs, numel( pet_vol));
        % ace_bval = nan( 2, nbsubjs, 3);
        ace_bval = [];
        for k = 1 : nbsubjs
            % spatial regression
            pet_receptor_names = pet_atlas_name( pet_ind);
            pet_mat = arrayfun( @(x) pet_vol{ x}.vol( gray_mask_ind), pet_ind, 'un', 0);
            pet_mat = cell2mat( pet_mat);

            nasal_falff = ts( :, k+nbsubjs);
            oral_falff = ts( :, k);


            [X, X_nam] = GetInteraction( pet_mat, pet_receptor_names);
            b = regress( nasal_falff, X);
            ace_bval( 1, k, :) = b( 2:end);

            b = regress( oral_falff, X);
            ace_bval( 2, k, :) = b( 2:end);
            X_nam = X_nam( 2:end);

                
            % % fitlm modelling
            % b_tbl = array2table( [pet_mat, nasal_falff], 'VariableNames', [pet_receptor_names(:)', {'ALFF'}]);
            % % mdl = fitlm( b_tbl, 'ALFF~NMDA+mGluR5+mGluR5*NMDA')
            % mdl = fitlm( b_tbl, 'interactions')
            % ace_bval( 1, k, :) = mdl.Coefficients.Estimate( 2:end);
            % % % mdl = fitglm( b_tbl, 'ALFF ~ NMDA + mGluR5 + NMDA*mGluR5')
            % 
            % b_tbl = array2table( [pet_mat, oral_falff], 'VariableNames', [pet_receptor_names(:)', {'ALFF'}]);
            % % ace_bval( 2, k, :) = mdl.Coefficients.Estimate( 2:end);
            % % mdl = fitlm( b_tbl, 'ALFF~NMDA+mGluR5+mGluR5*NMDA')
            % mdl = fitlm( b_tbl, 'interactions')
            % ace_bval( 2, k, :) = mdl.Coefficients.Estimate( 2:end);
            % 
            % X_nam = mdl.CoefficientNames( 2:end);
        end


        result_txt = [result_txt; {'  '}; {file_name}];
        figure
        for k = 1 : size( ace_bval, 3)
            subplot( 2, 4, k);
            [tmp_p, tmp_t, s] = Plot_Ttest( ace_bval( :, :, k)');
            txt = sprintf( '%s: %s', X_nam{ k}, s);
            title( txt);
            set( gca, 'xticklabel', {'Nasal', 'Oral'});
            ylabel( 'Beta values');

            result_txt = [result_txt; {txt}];
            pval = [pval; tmp_p];
            tval = [tval; tmp_t.tstat];
        end
        
        exportgraphics( fullfile( fig_dir, sprintf( '%s_%s.pdf', file_names{ file_idx}, strjoin( pet_receptor_names, '_')));

        writetable( cell2table( result_txt), ...
            fullfile( fig_dir, sprintf( '%s_%s.txt', file_names{ file_idx}, strjoin( pet_receptor_names, '_'))), ...
            'WriteVariableNames', false);
    end

    fdr( pval)

end




%%
% % Python toolbox: https://github.com/ottaviadipasquale/react-fmri

% clear; clc
% save_dir = '/Volumes/GZ_Backup/RestfMRI_REACT_24Motion';
% % each row is full-path-to-preprocessed fmri.nii.gz
% subj_list = fullfile( save_dir, 'subject_list.txt');
% subjs = importdata( subj_list);% 
% subjs = cellfun( @(x) regexp( x, 'ro_subj[0-9]*\/[a-zA-Z0-9]*\/', 'match'), subjs);
% subjs = cellfun( @(x) strrep( x, '/', '_'), subjs, 'un', 0);
% 
% % gm_file = '~/Desktop/Toolbox/react-fmri-main/data/gm_mask.nii.gz';
% gm_file = '~/Desktop/Data/MNIROIs/avg152T1_gray_thresh100.nii.gz';
% pet_atlas = G_UFind( '~/Desktop/test/REACT/PET_Atlas/CombineByReceptor_Normalized', 'name', 'Norm*.nii.gz');
% mask_savedir = fullfile( save_dir, 'React_mask');
% mask_file = REACT_Mask( mask_savedir, subj_list, pet_atlas, gm_file);
% 
% save( fullfile( save_dir, 'mask_file'), 'mask_file', 'subjs');
% REACT_glm( save_dir, subj_list, pet_atlas, mask_file, 'save_name', subjs);

