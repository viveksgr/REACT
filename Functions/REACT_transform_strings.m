function formatted_strings = REACT_transform_strings(file_paths)
    % Transforms file paths into formatted strings like "subject_02_nasal"
    %
    % Input:
    %   file_paths - cell array of file paths (strings)
    %
    % Output:
    %   formatted_strings - cell array of formatted strings
    
    % Initialize output
    formatted_strings = cell(size(file_paths));
    
    for i = 1:length(file_paths)
        % Extract subject ID
        subject_match = regexp(file_paths{i}, 'ro_subj(\d+)', 'tokens');
        subject_id = sprintf('%02d', str2double(subject_match{1}{1})); % Format to 2 digits
        
        % Extract condition (nasal or mouth)
        condition_match = regexp(file_paths{i}, '(Nasal|Mouth)', 'match', 'once');
        condition = lower(condition_match); % Convert to lowercase
        
        % Combine into the desired format
        formatted_strings{i} = sprintf('subject_%s_%s', subject_id, condition);
    end
end
