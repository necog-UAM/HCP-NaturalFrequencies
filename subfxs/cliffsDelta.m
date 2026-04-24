function delta = cliffsDeltaRepeated(data1, data2)
    % Ensure data1 and data2 have the same number of elements (one per subject)
    if length(data1) ~= length(data2)
        error('Data vectors must have the same length.');
    end

    % Compute all pairwise comparisons
    n = length(data1);
    count = 0;
    total = n;

    for i = 1:n
        if data1(i) > data2(i)
            count = count + 1; % Count cases where Condition 1 > Condition 2
        elseif data1(i) < data2(i)
            count = count - 1; % Count cases where Condition 1 < Condition 2
        end
    end

    % Compute Cliff's Delta
    delta = count / total;
end