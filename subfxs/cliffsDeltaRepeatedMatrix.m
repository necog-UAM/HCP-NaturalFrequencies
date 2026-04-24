function deltas = cliffsDeltaRepeatedMatrix(subj_condition1, subj_condition2)
    % Ensure the input matrices have the same dimensions
    if ~isequal(size(subj_condition1), size(subj_condition2))
        error('Input matrices must have the same dimensions.');
    end

    % Get the number of subjects and frequencies
    [numSubjects, numFrequencies] = size(subj_condition1);

    % Initialize a vector to store Cliff's Delta for each frequency
    deltas = zeros(1, numFrequencies);

    % Loop over each frequency
    for freq = 1:numFrequencies
        % Extract the data for the current frequency
        data1 = subj_condition1(:, freq); % Condition 1 for this frequency
        data2 = subj_condition2(:, freq); % Condition 2 for this frequency

        % Compute Cliff's Delta for the repeated-measures design
        count = 0;
        total = numSubjects;

        for i = 1:numSubjects
            if data1(i) > data2(i)
                count = count + 1; % Count cases where Condition 1 > Condition 2
            elseif data1(i) < data2(i)
                count = count - 1; % Count cases where Condition 1 < Condition 2
            end
        end

        % Compute Cliff's Delta for this frequency
        deltas(freq) = count / total;
    end
end