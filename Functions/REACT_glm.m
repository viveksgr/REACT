function REACT_glm( save_dir, subject_list, pet_atlas_file, stage_mask, varargin)
% clear; clc;
% %b1 = importdata( '~/Desktop/test/rfMRI_subj001_Nasal_CombinedPet/REACT_5HT1b/ro_subj001_react_stage1.txt');
% %t1 = MRIread( '~/Desktop/test/rfMRI_subj001_Nasal_CombinedPet/REACT_5HT1b/ro_subj001_react_stage2_map1.nii.gz');
%
% stdbrain_file = '/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz';
% subj_list_file = '~/Desktop/test/rfMRI_subj001_Nasal/subject_list.txt';
% pet_atlas_file = '~/Desktop/test/rfMRI_subj001_Nasal/Pet_Image/Combined/Normazlied/5HT1b_combined_2mm_atlas.nii.gz';
%
% stage1_mask_file = '~/Desktop/test/rfMRI_subj001_Nasal_CombinedPet/Test/mask_stage1.nii.gz';
% stage2_mask_file = '~/Desktop/test/rfMRI_subj001_Nasal_CombinedPet/Test/mask_stage2.nii.gz';
%
% save_dir = '~/Desktop/test/rfMRI_subj001_Nasal_CombinedPet/Test';
%
% See also REACT_Normalize REACT_Mask
%
% GZ
%
% save_dir, parent saving directory
% 'subject_list', '',... % full-path-to subject list, each row is one file
% 'pet_atlas', full-path-to-pet-atlas.nii.gz or cell array
% stage_mask, data structure with the following fields:
%       stage1_mask, full-path-to-stage1_mask.nii.gz 
%       stage2_mask, full-path-to-stage2_mask.nii.gz 
%   or cell array of data structure, must be same length to pet_atlas.
% 


opt = struct( 'mnibrain', '', ...
    'norm_method', 'demean', ... 'demean' | 'zscore'
    'save_name', '',...
    'nuisance', [],... % nuisance variable such as motion
    'param', []) ;

opt = G_SparseArgs( opt, varargin);


% zscore or demean fmri time series across time.
norm_meth = lower( opt.norm_method);
if ~ismember( norm_meth, {'demean', 'zscore', 'none'})
    error( 'unknown normalization method.');
end

subj_list = importdata( subject_list);
tmp = cellfun( @(x) isfile( x) || isfile( [x, '.nii']) || isfile( [x, '.nii.gz']), subj_list);
if ~all( tmp)
    error( 'The following frmri files were not found.\n%s', strjoin( subj_list( ~tmp), '\n'));
end


if ~iscell( pet_atlas_file)
    all_pet_atlas_file = {pet_atlas_file};
else
    all_pet_atlas_file = pet_atlas_file;
end

if isstruct( stage_mask)
    stage_mask = {stage_mask};
end

if numel( stage_mask) ~= numel( all_pet_atlas_file)
    error( 'The length of pet atlas must be equal to the number of stage masks.');
end


pet_atlas = MRIread( all_pet_atlas_file{1});
vol_res = pet_atlas.volres( 1:3);
vol_res = vol_res(:);
if isempty( opt.mnibrain)
    fsldir = getenv( 'FSLDIR');
    if isempty( fsldir)
        error( 'mnibrain is required.');
    end

    if all( vol_res == [1;1;1])
        mni_brain = 'MNI152_T1_1mm_brain.nii.gz';

    elseif all( vol_res == [2;2;2])
        mni_brain = 'MNI152_T1_2mm_brain.nii.gz';

    elseif all( vol_res == [0.5;0.5;0.5])
        mni_brain = 'MNI152_T1_0.5mm.nii.gz';

    else
        error( 'No MNI brain was found.');
    end

    opt.mnibrain = fullfile( fsldir, 'data', 'standard', mni_brain);
end

stdbrain = MRIread( opt.mnibrain);

for sidx = 1 : length( subj_list)
    if ~isempty( opt.nuisance)
        nui_mat = opt.nuisance{ sidx};
    else
        nui_mat = [];
    end

    in_fmri = subj_list{ sidx};
    if isempty( opt.save_name)
        [~, subj_nam, ext] = fileparts( in_fmri);
        subj_nam = extractBefore( [subj_nam, ext], '.');
    else
        subj_nam = opt.save_name{ sidx};
    end

    vol_fmri = MRIread( in_fmri);
    nb_fmri_vol = size( vol_fmri.vol, 4);
    rsfmri = reshape( vol_fmri.vol, [], nb_fmri_vol);

    for pet_idx = 1 : numel( all_pet_atlas_file)
        pet_atlas_file = all_pet_atlas_file{ pet_idx};

        if ~isfile( pet_atlas_file)
            error( 'PET atlas file was not found.');
        end
        [~, pet_atlas_nam] = fileparts( pet_atlas_file);
        pet_atlas_nam = extractBefore( pet_atlas_nam, '.');
        pet_save_dir = G_Fullfile( save_dir, ['REACT_', pet_atlas_nam]);

        if ~isfile( stage_mask{ pet_idx}.stage1_mask) || ~isfile( stage_mask{ pet_idx}.stage2_mask)
            error( 'Stage mask file was not found.');
        end

        stage1_mask = MRIread( stage_mask{ pet_idx}.stage1_mask);
        stage1_mask = stage1_mask.vol > 0;

        stage2_mask = MRIread( stage_mask{ pet_idx}.stage2_mask);
        stage2_mask = stage2_mask.vol > 0;


        pet_atlas = MRIread( pet_atlas_file);
        nb_pet_vol = size( pet_atlas.vol, 4);
        pet = reshape( pet_atlas.vol, [], nb_pet_vol);

        % demean, before masking?
        y = (rsfmri - mean( rsfmri, 1));
        x = (pet - mean( pet, 1));

        y = y( stage1_mask(:), :);
        x = x( stage1_mask(:), :);

        beta1 = nan( nb_fmri_vol, 1);
        for k = 1 : nb_fmri_vol
            tmp = regress( y( :, k), [ones( size( x)), x]);
            beta1( k, 1) = tmp( 2);
            % beta1( k, 1) = single( tmp( 2));
        end

        % save beta1 as regressor from stage 1
        writematrix( beta1, fullfile( pet_save_dir, [subj_nam, '_react_stage1_beta1.txt']));
        % isequal( single( beta1), single( b1))

        y = rsfmri';
        if strcmpi( norm_meth, 'zscore')
            y = (y - mean( y, 1)) ./ std( y, 0, 1);
        elseif strcmpi( norm_meth, 'demean')
            y = y - mean( y, 1);
        else
            % do nothing
        end

        y = y( :, stage2_mask(:));

        % z-score normalize beta1
        x = (beta1 - mean( beta1)) / std( beta1);

        nb_vox = size( y, 2);
        beta2 = nan( nb_vox, nb_pet_vol);
        for k = 1 : nb_vox
            tmp = regress( y( :, k), [ones( size( x)), x, nui_mat]);
            beta2( k, :) = tmp( 2:size( x, 2)+1);
        end

        for k = 1 : nb_pet_vol
            tmp = stdbrain;
            tmp.vol = zeros( size( tmp.vol));
            tmp.vol( stage2_mask) = beta2( :, k);
            save_nam = sprintf( '%s_react_stage2_map%d.nii.gz', subj_nam, k);
            MRIwrite( tmp, fullfile( pet_save_dir, save_nam), 'double');
        end
    end % 3/4d pet images

end % subject loop

