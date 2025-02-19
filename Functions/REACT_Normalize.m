function out_file = REACT_Normalize( x, varargin)
% Matlab implementation of REACT normalization of 
% https://github.com/ottaviadipasquale/react-fmri
% Note. react-fmri uses a single-precision.
% 
% Normalize pet image to 0-1 for all positive values.
% 
% Usage
%   out_file = REACT_Normalize( x);
%   out_file = REACT_Normalize( x, KEY, VALUE);
% 
% Input
%   x, MRIread data structure or full-path-to-nifti-file (including extension)
% 
%   Key-value
%       'intersect', 'yes' | 'no' (default); intersecting all volumes if x is a 4d
%       'average', 'yes' | 'no' (default); average across all images after normalization
%       'savefile', ''; full-path-to-file-to-save.nii.gz
% 
% See also REACT_Mask REACT_glm
% 
% GZ

opt = struct( 'intersect', 'no', ...
    'average', 'no', ...
    'savefile', '', ...
    'param', []);

opt = G_SparseArgs( opt, varargin);

if ~isfile( x)
    error( 'The following file does not exist.\n    %s', x);
end

% saving file
out_file = opt.savefile;
if isempty( out_file)
    if endsWith( x, '.nii')
        ext = '.nii';
    elseif endsWith( x, '.nii.gz')
        ext = '.nii.gz';
    else
        ext = '';
    end
    
    [p, nam, nam_ext] = filepart( x);
    nam = [nam, nam_ext];
    nam = regexprep( nam, [ext, '$'], '');
    out_file = fullfile( p, [nam, '_atlas.nii.gz']);
end

mri = MRIread( x);

% % react-fmri used a single precision
% mri.vol = single( mri.vol);

% Normalize positive non-zeros to [0, 1]
inters_img = lower( opt.intersect);
if ~ismember( inters_img, {'yes', 'no'})
    error( 'intersect must be yes or no');
end

do_avg = lower( opt.average);
if ~ismember( do_avg, {'yes', 'no'})
    error( 'average must be yes or no');
end

vol_dim = size( mri.vol);
if strcmpi( inters_img, 'yes')
    mask = all( mri.vol > 0, 4);
else
    mask = ones( vol_dim( 1:3));
end

nb_vol = size( mri.vol, 4);
for k = 1 : nb_vol
    m = mri.vol( :, :, :, k);
    zero_mask = m <= 0;
    zero_mask = zero_mask | ~mask;

    m( zero_mask) = 0;
    ind = find( ~zero_mask);
    lo = min( m( ind));
    hi = max( m( ind));
    m( ind) = (m( ind) - lo) / (hi - lo);
    mri.vol( :, :, :, k) = m;
end

if strcmpi( do_avg, 'yes')
    mri.vol = mean( mri.vol, 4);
end

% % avoid overwriting
% out_file = G_FileExist( out_file, 'ext', ext);
MRIwrite( mri, out_file, 'double');

% pet_atlas = MRIread( pet_atlas_file);
% isequal( mri.vol, python_atlas.vol)

end % function
