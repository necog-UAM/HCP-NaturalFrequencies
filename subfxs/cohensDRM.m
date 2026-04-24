function cohen_d_rm = cohensDRM(subj_condition1, subj_condition2)
    % COHENS_DRM Calculates Cohen's d for a repeated-measures design
    %
    % Inputs:
    %   subj_condition1 - Data for Condition 1 (e.g., motor task)
    %   subj_condition2 - Data for Condition 2 (e.g., fixation)
    %
    % Output:
    %   cohen_d_rm - Cohen's d for the repeated-measures design

    % Ensure the input matrices have the same dimensions
    if ~isequal(size(subj_condition1), size(subj_condition2))
        error('Input matrices must have the same dimensions.');
    end

    % Calculate the difference scores
    differences = subj_condition1 - subj_condition2;

    % Calculate the mean of the differences
    mean_diff = mean(differences, 1);

    % Calculate the standard deviation of the differences
    std_diff = std(differences, [], 1);

    % Compute Cohen's d for repeated measures
    cohen_d_rm = mean_diff ./ std_diff;
end