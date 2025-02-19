function filtered_files = filter_image_files(image_files,string_search,numpc)
    % Extract unique subject prefixes and count how many maps each has
    subject_map_counts = containers.Map(); 
    subject_files = containers.Map(); 

    if strcmp(string_search,'nasal')
        % Regular expression to extract subject prefix and map index
        pattern = '^(subject_\d+)_nasal_react_stage2_map(\d+).nii$';
    elseif strcmp(string_search,'mouth')
        pattern = '^(subject_\d+)_mouth_react_stage2_map(\d+).nii$';
    end

    for i = 1:length(image_files)
        file_name = image_files(i).name;
        tokens = regexp(file_name, pattern, 'tokens');
        if ~isempty(tokens)
            subj_prefix = tokens{1}{1};
            map_index = str2double(tokens{1}{2});

            % Store the map indices for each subject
            if ~isKey(subject_files, subj_prefix)
                subject_files(subj_prefix) = [];
            end
            subject_files(subj_prefix) = [subject_files(subj_prefix), map_index];
        end
    end

    % Find the minimum number of maps across all subjects
    min_map_count = inf;
    subjects = keys(subject_files);
    for i = 1:length(subjects)
        maps = subject_files(subjects{i});
        min_map_count = min(min_map_count, length(maps));
    end

    if nargin>2
        min_map_count = numpc;
    end

    % Create a new filtered list with equal number of maps per subject
    filtered_files = {};
    for i = 1:length(subjects)
        subj_prefix = subjects{i};
        sorted_maps = sort(subject_files(subj_prefix)); % Sort indices
        selected_maps = sorted_maps(1:min_map_count);  % Keep only min count

        for j = 1:length(selected_maps)
            map_idx = selected_maps(j);
            file_name = sprintf('%s_nasal_react_stage2_map%d.nii', subj_prefix, map_idx);
            filtered_files{end+1} = file_name; % Append to list
        end
    end
end
