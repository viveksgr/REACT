function coreg_images(root_dir,T1_file,fname)

 matlabbatch = [];
% List all neurotransmitter directories
neurotransmitter_dirs = dir(fullfile(root_dir, 'react_mask_Normalized*'));
neurotransmitter_dirs = neurotransmitter_dirs([neurotransmitter_dirs.isdir]);

file_list = cell(length(neurotransmitter_dirs),1);
for ff = 1:length(file_list)
    file_list{ff} = sprintf('%s,1',fullfile(neurotransmitter_dirs(ff).folder,neurotransmitter_dirs(ff).name,fname));
end

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {sprintf('%s,1',T1_file)};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = file_list(1);
matlabbatch{1}.spm.spatial.coreg.estwrite.other = file_list(2:end);
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

spm_jobman('run', matlabbatch);