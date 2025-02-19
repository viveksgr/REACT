function [var_mat, var_nam] = GetInteraction( mat, nam)
% x, N x M matrix
% y, 1 x M cell array



% pet_ind = [7 20];
% pet_name = pet_atlas_name( pet_ind);
% 
% pet_mat = arrayfun( @(x) pet_vol{ x}.vol( gray_mask_ind), pet_ind, 'un', 0);
% pet_mat = cell2mat( pet_mat);
% 
% 
% nasal_falff = ts( :, k+nbsubjs);
% oral_falff = ts( :, k);

% 
% mat = pet_mat;
% nam = pet_name;

[N, M] = size( mat);

var_mat = [];
var_nam = [];
for k = 2 : M
    ind = nchoosek( 1:M, k);
    nb_combs = size( ind, 1);

    for cidx = 1 : nb_combs 
        cn = nam{ ind( cidx, 1)};
        v = mat( :, ind( cidx, 1));
        for j = 2 : k
            v = v .* mat( :, ind( cidx, j));
            cn = sprintf( '%s*%s', cn, nam{ ind( cidx, j)});
        end
        var_nam = [var_nam, {cn}];

        var_mat = [var_mat, v];
    end

end

var_mat = [ones( N, 1), mat, var_mat];
var_nam = cat( 2, {'Constant'}, nam(:)', var_nam);



% X = [ones( numel( gray_mask_ind), 1), ...
%     pet_mat, pet_mat(:, 1) .* pet_mat(:, 2)];
% X_name = cat( 2, pet_name', [pet_name{1}, '_', pet_name{2}]);





% end

