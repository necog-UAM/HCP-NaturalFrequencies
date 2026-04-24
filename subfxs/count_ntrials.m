load('ntrials_Wrkmem.mat')
% Compute mean per column
col_mean = mean(ntrials, 1);

% Compute standard deviation per column
col_std = std(ntrials, 0, 1);

% Display results
disp('Column means:')
disp(col_mean)

disp('Column standard deviations:')
disp(col_std)

% Compute column-wise mean and standard deviation
col_mean = mean(ntrials, 1);
col_std  = std(ntrials, 0, 1);

% Get ranges
mean_range = [min(col_mean) max(col_mean)];
std_range  = [min(col_std)  max(col_std)];

% Display results
fprintf('Range of column means: %.4f – %.4f\n', mean_range(1), mean_range(2));
fprintf('Range of column standard deviations: %.4f – %.4f\n', std_range(1), std_range(2));
