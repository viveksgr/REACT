% Define subject info
% subject_id = 'sub-001';  % Change for each subject
% data_dir = 'C:\fMRI_data';  % Root data directory
% output_dir = fullfile(data_dir, subject_id, 'first_level'); % Output directory
% func_dir1 = fullfile(data_dir, subject_id, 'func', 'condition1'); % Path to Condition 1 scans
% func_dir2 = fullfile(data_dir, subject_id, 'func', 'condition2'); % Path to Condition 2 scans
% realignment_file1 = fullfile(func_dir1, 'rp.txt'); % Motion regressors for session 1
% realignment_file2 = fullfile(func_dir2, 'rp.txt'); % Motion regressors for session 2

ss = 2;
TR = 0.555; % Repetition time (modify if needed)


% Get functional images
str_nasal = sprintf('F:\\RestOlfAnaly\\SubjectData\\ro_subj0%02d\\Nasal10Min\\Nasal10Min.feat',ss);
str_mouth = sprintf('F:\\RestOlfAnaly\\SubjectData\\ro_subj0%02d\\Mouth10Min\\Mouth10Min.feat',ss);
func_images1_lab = fullfile(str_nasal,'filtered_func_data_dspk_WmCSFWb24Motion_globalnorm_rlt_filtBP008012std_nonlinear_s3.nii'); 
func_images2_lab = fullfile(str_mouth,'filtered_func_data_dspk_WmCSFWb24Motion_globalnorm_rlt_filtBP008012std_nonlinear_s3.nii');
output_dir = 'C:\Work\REACT\Results\output_dir';
output_dir_sub = fullfile(output_dir,sprintf('sub%02d',ss));
mkdir(output_dir_sub)

% Expand 4D NIfTI to a list of 3D volumes
mkdir(fullfile(str_nasal,'3d_vols'))
func_images1_lab2 = spm_file_split_modified(func_images1_lab,fullfile(str_nasal,'3d_vols'));
func_images1 = cellfun(@(x) strcat(x,',1'),func_images1_lab2,'UniformOutput',false);

mkdir(fullfile(str_mouth,'3d_vols'))
func_images2_lab2 = spm_file_split_modified(func_images2_lab,fullfile(str_mouth,'3d_vols'));
func_images2 = cellfun(@(x) strcat(x,',1'),func_images2_lab2,'UniformOutput',false);


%% Step 1: Model Specification
matlabbatch = {};
matlabbatch{1}.spm.stats.fmri_spec.dir = {output_dir_sub};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

% Session 1 (Condition 1)
matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = func_images1;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''}; % No task-related regressors
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {''}; % Motion regressors
matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128; % High-pass filter

% Session 2 (Condition 2)
matlabbatch{1}.spm.stats.fmri_spec.sess(2).scans = func_images2;
matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi = {''}; % No task-related regressors
matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi_reg = {''}; % Motion regressors
matlabbatch{1}.spm.stats.fmri_spec.sess(2).hpf = 128; % High-pass filter

% Specify model settings
matlabbatch{1}.spm.stats.fmri_spec.fact = struct([]); % No factorial design
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; % Use canonical HRF
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mask = {''}; % No explicit mask
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'none';
spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);
disp('First-level BOLD contrast analysis completed.');


%% Step 2: Model Estimation
matlabbatch = {};
matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(output_dir_sub, 'SPM.mat')};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

% Step 3: Define Contrast (Condition 1 - Condition 2)
matlabbatch{2}.spm.stats.con.spmmat = {fullfile(output_dir_sub, 'SPM.mat')};

% Condition1 > Condition2
matlabbatch{2}.spm.stats.con.consess{1}.tcon.name = 'Condition1 > Condition2';
matlabbatch{2}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{2}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

% Condition2 > Condition1
matlabbatch{2}.spm.stats.con.consess{2}.tcon.name = 'Condition2 > Condition1';
matlabbatch{2}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{2}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

matlabbatch{2}.spm.stats.con.delete = 1; % Delete previous contrasts

% Run the batch
spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);
disp('First-level BOLD contrast analysis completed.');
