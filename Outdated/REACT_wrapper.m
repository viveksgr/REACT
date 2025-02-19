%%% % Python toolbox: https://github.com/ottaviadipasquale/react-fmri

% save_dir = 'E:\RestOlfAnaly\Results\RestfMRI_REACT_24Motion';
% Actual subject data: E:\RestOlfAnaly\SubjectData\ro_subj002\Nasal10Min\Nasal10Min.feat

save_dir = 'C:\Work\REACT\Data\Prelim_res';
subj_list = fullfile( save_dir, 'subject_list.txt');
% subjs = importdata( subj_list);
subjs = {'S2_n','S2_m','S3_n','S3_m'};

gm_mask_file = 'C:\Work\REACT\Data\avg152T1_gray_thresh100.nii';

fdir = dir( fullfile('C:\Work\REACT\Data\PET_Atlas\CombineByReceptor_Normalized', 'Norm*.nii.gz'));
files = fullfile({fdir.folder}, {fdir.name});
files = string(files);
pet_atlas_file = cellstr(files);

mask_savedir = fullfile( save_dir, 'React_mask');
mask_file = REACT_Mask_SPM(mask_savedir, subj_list, pet_atlas_file, gm_mask_file);

REACT_glm_win( save_dir, subj_list, pet_atlas_file, mask_file, subjs);

