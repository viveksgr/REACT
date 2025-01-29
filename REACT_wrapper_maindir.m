%%% % Python toolbox: https://github.com/ottaviadipasquale/react-fmri

% save_dir = 'E:\RestOlfAnaly\Results\RestfMRI_REACT_24Motion';
% Actual subject data: E:\RestOlfAnaly\SubjectData\ro_subj002\Nasal10Min\Nasal10Min.feat

save_dir = 'C:\Work\REACT\Data\group_res - adj thr';
atlas_dir = 'C:\Work\REACT\Data\PET_Atlas';
subj_list = fullfile( save_dir, 'subject_list_full.txt');

subjs = importdata( subj_list);
subjs = REACT_transform_strings(subjs);

gm_mask_file = 'C:\Work\REACT\Data\avg152T1_gray_thresh100.nii';

fdir = dir( fullfile('C:\Work\REACT\Data\PET_Atlas\CombineByReceptor_Normalized', 'Norm*.nii.gz'));
files = fullfile({fdir.folder}, {fdir.name});
files = string(files);
pet_atlas_file = cellstr(files);

mask_savedir = fullfile( save_dir, 'React_mask');
thresh = [0.05 0.05 0.1 0.05 0.1 0.005 0.02 0.05 0.05 0.05...
    0.1 0.2 0.1 0.05 0.1 0.05 0.05 0.05 0.05 0.2];

% Compute masks
% mask_file = REACT_Mask_SPM(mask_savedir, subj_list, pet_atlas_file, gm_mask_file,thresh);

% Compute GLM
% REACT_glm_win( save_dir, subj_list, pet_atlas_file, mask_file, subjs);

% Compute composite GLM - not functional
% REACT_glm_multi_template( save_dir, subj_list, pet_atlas_file, mask_file, subjs);

% Whole brain t-maps
[results, groupstats, corrvox] = neurotransmitter_analysis( save_dir);

% Corregistration
% T1_file = fullfile(atlas_dir,'mniT1.nii');
% coreg_images(save_dir,T1_file,'tmap_neg.nii')

% ROI based maps
opt.num = 90;
aal_file = fullfile(atlas_dir,'raal.nii');
aal_labels = fullfile(atlas_dir,'aal.txt');
opt.fname = 'rtmap_pos.nii';
Tmat = ROIbased_analysis(save_dir,aal_file,aal_labels,opt);

%% Plotting
neurolabels = cellfun(@(x) x(23:end-10), groupstats.map_labels,'UniformOutput', false);

% p-value adjustment
cluster_labels = kmeans(Tmat,4, 'Distance', 'sqeuclidean', 'Replicates', 10);
[labels, sorted_indices] = sort(cluster_labels);
sortedmat = Tmat(sorted_indices,:);

figure()
imagesc(sortedmat )
xticks(1:1:45)
aal_labels_list = importdata(aal_labels);
xticklabels(aal_labels_list.textdata(1:2:end, 2))
ax = gca;
ax.TickLabelInterpreter = 'none';
xtickangle(90)
yticks(1:1:20)
yticklabels(neurolabels(sorted_indices))

% p_val = arrayfun(@(x) 2*(1 - tcdf(abs(x), 17)),sortedmat); 
figure()
sig_mat = sortedmat;
p_val = sig_mat;
sig_mat(abs(p_val)<1.975)=0;
imagesc(sig_mat)
xticks(1:1:45)
aal_labels_list = importdata(aal_labels);
xticklabels(aal_labels_list.textdata(1:2:end, 2))
ax = gca;
ax.TickLabelInterpreter = 'none';
xtickangle(90)
yticks(1:1:20)
yticklabels(neurolabels(sorted_indices))

%% figure()
% hold on
% bar(groupstats.map_score*100)
% xticks(1:length(groupstats.map_score))
% xtickangle(90)
% xticklabels(neurolabels)
% ylabel('% significant voxels')
% 
% % figure
% methods = 'normedsignificance';
% switch methods
%     case 'kmeans'
%         cluster_labels = kmeans(corrvox.p',4, 'Distance', 'correlation', 'Replicates', 10);
%         [~, sorted_indices] = sort(cluster_labels);
%         ad = corrvox.p(:,sorted_indices);
%         sorted_adj = corrcoef(ad);
% 
%     case 'spectral'
%         distmat = corrcoef(corrvox.p);
%         [spec, sorted_indices] = spectral_reorder(1-distmat);
%         sorted_adj = 1-spec;
% 
%     case 'significance'
%         [sorted_adj, sorted_indices] = cluster_templates_by_signed_significance(corrvox.p, 1.65, 3);
% 
%     case 'normedsignificance'
%         [sorted_adj, sorted_indices] = cluster_templates_by_normalized_significance(corrvox.p, 1.65, 3);
% end
% 
% xticks(1:length(groupstats.map_score))
% xtickangle(90)
% xticklabels(neurolabels(sorted_indices))
% yticks(1:length(groupstats.map_score))
% % ytickangle(90)
% yticklabels(neurolabels(sorted_indices))
% xlabel('Receptors')
% ylabel('Receptors')



