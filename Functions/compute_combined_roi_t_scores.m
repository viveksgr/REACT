function median_t_scores = compute_combined_roi_t_scores( file_tmap,  aal_masks,method)
    % Compute median t-scores for combined ROIs based on the AAL atlas.
    %
    % Input:
    %   file_tmap: Path to the 3D NIfTI file with t-scores
    %   file_aal_mask: Path to the 3D NIfTI AAL mask file
    %
    % Output:
    %   median_t_scores: Vector of median t-scores for 45 combined ROIs

    % Initialize output vector
    num_combined_rois = 45;
    median_t_scores = zeros(num_combined_rois, 1);

    % Iterate through combined ROIs
    for roi_idx = 1:num_combined_rois
        % Calculate ROI indices
        left_roi = 2 * roi_idx - 1; % Left hemisphere ROI
        right_roi = 2 * roi_idx;   % Right hemisphere ROI

        % Create combined ROI mask
        combined_mask = (aal_masks == left_roi) | (aal_masks == right_roi);

        % Extract t-scores for the combined ROI
        roi_t_scores = file_tmap(combined_mask);

        switch method
            case 'median_t'
                % Compute median t-score for the combined ROI
                median_t_scores(roi_idx) = median(roi_t_scores, 'omitnan'); % Omit NaNs

            case 'significance'
                median_t_scores(roi_idx) = sum(abs(roi_t_scores)>tinv(0.975,17))/length(roi_t_scores); % Omit NaNs

            case 'significance_unn'
                 median_t_scores(roi_idx) = sum(abs(roi_t_scores)>tinv(0.975,17)); % Omit NaNs

            case 'weighted_mean'
                % % Compute weighted mean t-score for the combined ROI
                weights = abs(roi_t_scores); % Use absolute t-scores as weights
                weighted_mean = sum(roi_t_scores .* weights, 'omitnan') / sum(weights, 'omitnan');
                median_t_scores(roi_idx) = weighted_mean;

            otherwise
                median_t_scores(roi_idx) = mean(roi_t_scores);

        end
    end
end
