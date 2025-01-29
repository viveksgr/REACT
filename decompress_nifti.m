function decompressed_file = decompress_nifti(gz_file)
    % Decompress a .nii.gz file to a .nii file using MATLAB
    % Input:
    %   gz_file - path to the .nii.gz file
    % Output:
    %   decompressed_file - path to the decompressed .nii file

    % Check if the file has .gz extension
    if endsWith(gz_file, '.gz')
        % Define the output file by removing the .gz extension
        decompressed_file = gz_file(1:end-3); 

        % Read the compressed file, uncompress, and save as .nii
        gunzip(gz_file, fileparts(decompressed_file)); 
    else
        error('File must have a .gz extension');
    end
end
