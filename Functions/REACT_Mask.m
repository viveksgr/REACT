function out_file = REACT_Mask( save_dir, subj_list, pet_atlas_file, gm_mask_file)
% Create masks for multivariate regression.
% 
% Input
%   save_dir = '~/Desktop/test/REACT';
%   subject list. Each row is a full-path-to-fMRI.nii.gz
%       subj_list = '~/Desktop/test/subject_list.txt';
% 
%   Normalized pet atlas, 4d or multiple 4d volumes
%   pet_atlas_file = '~/Desktop/test/Normazlied/5HT1b_combined_2mm_atlas.nii.gz';
%   Intersect across all 4d volumes for each file
%  
%   Gray matter mask
%   gm_mask_file = '~/Desktop/Toolbox/react-fmri-main/data/gm_mask.nii.gz';
% 
% Output
%   out_file, data structure with the following fields
%       stage1_mask
%       stage2_mask
% 
% See also REACT_Normalize REACT_glm
% 
% GZ


gray_mask = MRIread( gm_mask_file);
[gm_dim, gm_resol, gm_orient] = MRInfo( gray_mask);
if length( gm_dim) > 3
    error( 'Gray matter mask must be a 3d volume.');
end

% PET atlas
if ischar( pet_atlas_file) || isstring( pet_atlas_file)
    pet_atlas_file = {pet_atlas_file};
end

nb_pet_atlas = numel( pet_atlas_file);
pet_atlas_vol = cell( nb_pet_atlas, 1);
for k = 1 : nb_pet_atlas
    mri = MRIread( pet_atlas_file{ k});
    [mri_dim, mri_resol, mri_orient] = MRInfo( mri);
    if ~strcmpi( gm_orient, mri_orient)
        error( '%s has a different orientation from gray matter mask.', pet_atlas_file{ k});
    end

    if ~isequal( gm_dim, mri_dim( 1:3))
        error( '%s has a different dimension from gray matter mask', pet_atlas_file{ k});
    end

    if ~isequal( gm_resol, mri_resol( 1:3))
        error( '%s has a different resolution from gray matter mask', pet_atlas_file{ k});
    end

    pet_atlas_vol{ k} = mri.vol;
end

subj_list = importdata( subj_list);

% fMRI mask, all subjects have std > 0
fmri_mask = zeros( gm_dim);
nb_subj = length( subj_list);
for k = 1 : nb_subj
    fprintf( 'Reading fmri %d/%d\n', k, nb_subj);
    
    subj_file = subj_list{ k};    
    [p, nam, ext] = fileparts( subj_file);
    nam = [nam, ext];
    ext = extractAfter( nam, '.');
    nam = extractBefore( nam, '.');
    tstd_file = fullfile( p, [nam, '_Tstd.', ext]);
    if isfile( tstd_file)
        mri = MRIread( tstd_file);
    else
        mri = MRIread( subj_file);
    end

    [mri_dim, mri_resol, mri_orient] = MRInfo( mri);
    if ~strcmpi( gm_orient, mri_orient)
        error( '%s has a different orientation from gray matter mask.', subj_file);
    end

    if ~isequal( gm_dim, mri_dim( 1:3))
        error( '%s has a different dimension from gray matter mask', subj_file);
    end

    if ~isequal( gm_resol, mri_resol( 1:3))
        error( '%s has a different resolution from gray matter mask', subj_file);
    end

    if isfile( tstd_file)
        ind = mri.vol > 0;
    else
        ind = std( mri.vol, [], 4) > 0;
    end
    
    fmri_mask( ind) = fmri_mask( ind) + 1;
end

fmri_mask = fmri_mask == nb_subj;
stage2_mask = fmri_mask & (gray_mask.vol > 0);

if isempty( save_dir)
    p = fileparts( pet_atlas_file{1});
    save_dir = fullfile( p, 'React_mask');
else
    if ~exist( save_dir, 'dir')
        mkdir( save_dir);
    end
end

out_file = cell( nb_pet_atlas, 1);
for k = 1 : nb_pet_atlas
    % intersection of all pet masks
    pet_mask = sum( pet_atlas_vol{ k} > 0, 4) == size( pet_atlas_vol{ k}, 4);
    stage1_mask = stage2_mask & pet_mask;

    [~, pet_nam] = G_Fileparts( pet_atlas_file{ k}, '.nii.gz');
    cur_savedir = G_Fullfile( save_dir, ['react_mask_', pet_nam]);

    % write stage masks to file
    save_file = fullfile( cur_savedir, 'mask_stage1.nii.gz');
    m = gray_mask;
    m.vol = stage1_mask;
    MRIwrite( m, save_file, 'double');
    out_file{ k}.stage1_mask = save_file;

    save_file = fullfile( cur_savedir, 'mask_stage2.nii.gz');
    m = gray_mask;
    m.vol = stage2_mask;
    MRIwrite( m, save_file, 'double');
    out_file{ k}.stage2_mask = save_file;
end

end % function


function [d, r, orient] = MRInfo( x)
M = x.vox2ras0( 1:3, 1:3);
sz = sqrt( sum( M.^2, 1));
M = M ./ sz;
[~, loc] = max( abs( M), [], 1);

txt = {{'L', 'R'}, {'P', 'A'}, {'I', 'S'}};
M = sign( M);
orient = arrayfun( @(a) txt{ a}{ M( loc(a), a) == [-1, 1]}, 1:3, 'un', 0);
orient = strjoin( orient, '');

d = permute( x.volsize, [2, 3, 1, 4]);
d = d(:)';
r = permute( x.volres, [2, 3, 1, 4]);
end % function
