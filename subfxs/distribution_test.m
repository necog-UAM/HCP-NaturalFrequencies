%%
% Obtener todos los field names
fields = fieldnames(fnat_struct);

% Comparar Motort_Fixation contra todas las condiciones que contienen Motort_
motort_fields = fields(contains(fields, 'Motort_'));
reference_condition_Motort = 'Motort_Fixation';

for i = 1:length(motort_fields)
    current_condition = motort_fields{i};
    
    if strcmp(current_condition, reference_condition_Motort)
        continue;
    end
    
    cohen_d_per_voxel = NaN(1, 1925);
    T_sum_per_voxel = NaN(1, 1925);
    
    for vox = 1:1925
        %[T_sum, cohen_d_M] = HCP9_Naturalfreq_cohen_d_foi(fnat_struct, reference_condition_Motort, current_condition, vox, foi, neighbor_matrix, substsk);
        [T_sum, cohen_d_M] = HCP9_Naturalfreq_cohen_d(fnat_struct, reference_condition_Motort, current_condition, vox, foi1, foi1x, foi2, neighbor_matrix, substsk);
        cohen_d_per_voxel(vox) = cohen_d_M;
        T_sum_per_voxel(vox) = T_sum;
    end
    
    comparison_name = [reference_condition_Motort '_vs_' current_condition];
    cohen_d_struct.(comparison_name).T_sum = T_sum_per_voxel;
    cohen_d_struct.(comparison_name).cohen_D = cohen_d_per_voxel;
end
%%
data1 = subj_condition1(1,:);
data2 = subj_condition2(1,:);
% Crear un nuevo gráfico
figure;

% Pintar el histograma de la primera distribución
histogram(data1, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;

% Pintar el histograma de la segunda distribución
histogram(data2, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Añadir leyenda y etiquetas
legend('Distribución 1', 'Distribución 2');
xlabel('Valor');
ylabel('Densidad de probabilidad');
title('Comparación de dos distribuciones');

% Mostrar el gráfico
hold off;
%%
% Load your data
data1 = subj_condition1(1,:);
data2 = subj_condition2(1,:);

% Create a new figure
figure;

% Plot the histogram for the first distribution
histogram(data1, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;

% Plot the histogram for the second distribution
histogram(data2, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Calculate the means of the two distributions
mean_data1 = mean(data1);
mean_data2 = mean(data2);
median_data1 = prctile(data1, 50);
median_data2 = prctile(data2, 50);

% Plot the mean of the first distribution as a vertical line or marker
plot([mean_data1, mean_data1], [0, max(histcounts(data1, 'Normalization', 'pdf'))], 'b--', 'LineWidth', 1.5); % Red dashed line
% Plot the mean of the second distribution as a vertical line or marker
plot([mean_data2, mean_data2], [0, max(histcounts(data2, 'Normalization', 'pdf'))], 'r--', 'LineWidth', 1.5); % Green dashed line

% Plot the mean of the first distribution as a vertical line or marker
plot([median_data1, median_data1], [0, max(histcounts(data1, 'Normalization', 'pdf'))], 'b-', 'LineWidth', 1.5); % Red dashed line
% Plot the mean of the second distribution as a vertical line or marker
plot([median_data2, median_data2], [0, max(histcounts(data2, 'Normalization', 'pdf'))], 'r-', 'LineWidth', 1.5); % Green dashed line

% Add legend and labels
legend('Distribución 1', 'Distribución 2', 'Media Distribución 1', 'Media Distribución 2', 'Mediana 1', 'Mediana 2');
xlabel('Valor');
ylabel('Densidad de probabilidad');
title('Comparación de dos distribuciones con medias marcadas');

% Show the plot
hold off;
%%
% Load your data
data1 = subj_condition1(1,:);
data2 = subj_condition2(1,:);

% Create a new figure
figure;

% Plot the kernel density estimate for the first distribution
[f1, xi1] = ksdensity(data1);
plot(xi1, f1, 'b', 'LineWidth', 1.5);

hold on;

% Plot the kernel density estimate for the second distribution
[f2, xi2] = ksdensity(data2);
plot(xi2, f2, 'r', 'LineWidth', 1.5);

% Calculate the means and medians of the two distributions
mean_data1 = mean(data1);
mean_data2 = mean(data2);
median_data1 = prctile(data1, 50);
median_data2 = prctile(data2, 50);

% Plot the mean of the first distribution as a vertical line
plot([mean_data1, mean_data1], [0, max(f1)], 'b--', 'LineWidth', 1.5); % Blue dashed line

% Plot the mean of the second distribution as a vertical line
plot([mean_data2, mean_data2], [0, max(f2)], 'r--', 'LineWidth', 1.5); % Red dashed line

% Plot the median of the first distribution as a vertical line
plot([median_data1, median_data1], [0, max(f1)], 'b-', 'LineWidth', 1.5); % Blue solid line

% Plot the median of the second distribution as a vertical line
plot([median_data2, median_data2], [0, max(f2)], 'r-', 'LineWidth', 1.5); % Red solid line

% Add legend and labels
legend('Distribución 1', 'Distribución 2', 'Media Distribución 1', 'Media Distribución 2', 'Mediana 1', 'Mediana 2');
xlabel('Valor');
ylabel('Densidad de probabilidad');
title('Comparación de dos distribuciones usando KDE');

% Show the plot
hold off;

%%
% normality
diff = data1-data2;
[~, p_value] = lillietest(diff);
%%
% Load your data
data1 = subj_condition1(1,:);
data2 = subj_condition2(1,:);

% Create a new figure
figure;

% Plot the histogram for the first distribution (normalized to PDF)
histogram(data1, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Hist 1');
hold on;
histogram(data2, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Hist 2');

% Estimate parameters for Gaussian fit
mu1 = mean(data1); sigma1 = std(data1);
mu2 = mean(data2); sigma2 = std(data2);

% Generate x values for Gaussian curves
x_gauss1 = linspace(min(data1), max(data1), 100);
x_gauss2 = linspace(min(data2), max(data2), 100);

% Compute Gaussian PDFs
y_gauss1 = normpdf(x_gauss1, mu1, sigma1);
y_gauss2 = normpdf(x_gauss2, mu2, sigma2);

% Plot Gaussian fits
plot(x_gauss1, y_gauss1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Gaussian Fit 1');
plot(x_gauss2, y_gauss2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Gaussian Fit 2');

% Kernel density estimates
[f1, xi1] = ksdensity(data1);
[f2, xi2] = ksdensity(data2);

% Plot kernel density estimates
plot(xi1, f1, 'b--', 'LineWidth', 1.5, 'DisplayName', 'KDE 1');
plot(xi2, f2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'KDE 2');

% Log-normal fits
params_lognorm1 = lognfit(data1);
params_lognorm2 = lognfit(data2);

x_lognorm1 = linspace(min(data1), max(data1), 100);
x_lognorm2 = linspace(min(data2), max(data2), 100);

y_lognorm1 = lognpdf(x_lognorm1, params_lognorm1(1), params_lognorm1(2));
y_lognorm2 = lognpdf(x_lognorm2, params_lognorm2(1), params_lognorm2(2));

% Plot log-normal fits
plot(x_lognorm1, y_lognorm1, 'c-', 'LineWidth', 1.5, 'DisplayName', 'Log-Normal Fit 1');
plot(x_lognorm2, y_lognorm2, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Log-Normal Fit 2');

% Gamma fits
params_gamma1 = gamfit(data1);
params_gamma2 = gamfit(data2);

x_gamma1 = linspace(min(data1), max(data1), 100);
x_gamma2 = linspace(min(data2), max(data2), 100);

y_gamma1 = gampdf(x_gamma1, params_gamma1(1), params_gamma1(2));
y_gamma2 = gampdf(x_gamma2, params_gamma2(1), params_gamma2(2));

% Plot gamma fits
plot(x_gamma1, y_gamma1, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Gamma Fit 1');
plot(x_gamma2, y_gamma2, 'y-', 'LineWidth', 1.5, 'DisplayName', 'Gamma Fit 2');

% Beta fits (if data is between 0 and 1)
if all(data1 >= 0 & data1 <= 1) && all(data2 >= 0 & data2 <= 1)
    params_beta1 = betafit(data1);
    params_beta2 = betafit(data2);

    x_beta1 = linspace(0, 1, 100);
    x_beta2 = linspace(0, 1, 100);

    y_beta1 = betapdf(x_beta1, params_beta1(1), params_beta1(2));
    y_beta2 = betapdf(x_beta2, params_beta2(1), params_beta2(2));

    % Plot beta fits
    plot(x_beta1, y_beta1, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Beta Fit 1');
    plot(x_beta2, y_beta2, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Beta Fit 2');
end

% Calculate means and medians
mean_data1 = mean(data1); mean_data2 = mean(data2);
median_data1 = prctile(data1, 50); median_data2 = prctile(data2, 50);

% Plot means and medians
plot([mean_data1, mean_data1], [0, max(f1)], 'b:', 'LineWidth', 1.5, 'DisplayName', 'Mean 1');
plot([mean_data2, mean_data2], [0, max(f2)], 'r:', 'LineWidth', 1.5, 'DisplayName', 'Mean 2');
plot([median_data1, median_data1], [0, max(f1)], 'b-.', 'LineWidth', 1.5, 'DisplayName', 'Median 1');
plot([median_data2, median_data2], [0, max(f2)], 'r-.', 'LineWidth', 1.5, 'DisplayName', 'Median 2');

% Add legend and labels
legend('Location', 'best', 'FontSize', 8);
xlabel('Valor');
ylabel('Densidad de probabilidad');
title('Comparación de distribuciones con ajustes múltiples');

% Show the plot
hold off;


% 5. What If the Data Is Not Normally Distributed?
% If the normality assumption is violated, you have a few options:
% 
% Non-parametric Alternative : Use the Wilcoxon signed-rank test , which does not assume normality.
% Transformation : Apply a transformation (e.g., log, square root) to the data to make it more normal-like.
% Bootstrap Methods : Use resampling techniques like bootstrapping to estimate the significance of the test without assuming normality.
%%
% Load your data
data1 = subj_condition1(1,:);
data2 = subj_condition2(1,:);

% Create a new figure
figure;

% Plot the histogram for the first distribution (normalized to PDF)
histogram(data1, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Hist 1');
hold on;
histogram(data2, 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Hist 2');

% Estimate parameters for Gaussian fit
mu1 = mean(data1); sigma1 = std(data1);
mu2 = mean(data2); sigma2 = std(data2);

% Generate x values for Gaussian curves
x_gauss1 = linspace(min(data1), max(data1), 100);
x_gauss2 = linspace(min(data2), max(data2), 100);

% Compute Gaussian PDFs
y_gauss1 = normpdf(x_gauss1, mu1, sigma1);
y_gauss2 = normpdf(x_gauss2, mu2, sigma2);

% Plot Gaussian fits
plot(x_gauss1, y_gauss1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Gaussian Fit 1');
plot(x_gauss2, y_gauss2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Gaussian Fit 2');

% Kernel density estimates
[f1, xi1] = ksdensity(data1);
[f2, xi2] = ksdensity(data2);

% Plot kernel density estimates
plot(xi1, f1, 'b--', 'LineWidth', 1.5, 'DisplayName', 'KDE 1');
plot(xi2, f2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'KDE 2');

% % Log-normal fits
% params_lognorm1 = lognfit(data1);
% params_lognorm2 = lognfit(data2);
% 
% x_lognorm1 = linspace(min(data1), max(data1), 100);
% x_lognorm2 = linspace(min(data2), max(data2), 100);
% 
% y_lognorm1 = lognpdf(x_lognorm1, params_lognorm1(1), params_lognorm1(2));
% y_lognorm2 = lognpdf(x_lognorm2, params_lognorm2(1), params_lognorm2(2));
% 
% % Plot log-normal fits
% plot(x_lognorm1, y_lognorm1, 'c-', 'LineWidth', 1.5, 'DisplayName', 'Log-Normal Fit 1');
% plot(x_lognorm2, y_lognorm2, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Log-Normal Fit 2');

% Gamma fits
params_gamma1 = gamfit(data1);
params_gamma2 = gamfit(data2);

x_gamma1 = linspace(min(data1), max(data1), 100);
x_gamma2 = linspace(min(data2), max(data2), 100);

y_gamma1 = gampdf(x_gamma1, params_gamma1(1), params_gamma1(2));
y_gamma2 = gampdf(x_gamma2, params_gamma2(1), params_gamma2(2));

% Plot gamma fits
plot(x_gamma1, y_gamma1, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Gamma Fit 1');
plot(x_gamma2, y_gamma2, 'y-', 'LineWidth', 1.5, 'DisplayName', 'Gamma Fit 2');

% Beta fits (if data is between 0 and 1)
if all(data1 >= 0 & data1 <= 1) && all(data2 >= 0 & data2 <= 1)
    params_beta1 = betafit(data1);
    params_beta2 = betafit(data2);

    x_beta1 = linspace(0, 1, 100);
    x_beta2 = linspace(0, 1, 100);

    y_beta1 = betapdf(x_beta1, params_beta1(1), params_beta1(2));
    y_beta2 = betapdf(x_beta2, params_beta2(1), params_beta2(2));

    % Plot beta fits
    plot(x_beta1, y_beta1, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Beta Fit 1');
    plot(x_beta2, y_beta2, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Beta Fit 2');
end

% Calculate means and medians
mean_data1 = mean(data1); mean_data2 = mean(data2);
median_data1 = prctile(data1, 50); median_data2 = prctile(data2, 50);

% Plot means and medians
plot([mean_data1, mean_data1], [0, max(f1)], 'b:', 'LineWidth', 1.5, 'DisplayName', 'Mean 1');
plot([mean_data2, mean_data2], [0, max(f2)], 'r:', 'LineWidth', 1.5, 'DisplayName', 'Mean 2');
plot([median_data1, median_data1], [0, max(f1)], 'b-.', 'LineWidth', 1.5, 'DisplayName', 'Median 1');
plot([median_data2, median_data2], [0, max(f2)], 'r-.', 'LineWidth', 1.5, 'DisplayName', 'Median 2');

% Add legend and labels
legend('Location', 'best', 'FontSize', 8);
xlabel('Valor');
ylabel('Densidad de probabilidad');
title('Comparación de distribuciones con ajustes múltiples');

% Show the plot
hold off;
%%
% Calculate Cliff's Delta for each frequency
deltas = cliffsDeltaRepeatedMatrix(subj_condition1, subj_condition2);
% Plot Cliff's Delta values
figure;
plot(foi2, deltas, 'LineWidth', 1.5);
xlabel('Frequency Index');
ylabel("Cliff's Delta");
title('Cliff''s Delta Across Frequencies');
grid on;
%--------------------------------------------------------------------------
%% Gráfico tamaño del efecto (Cohen's d)
% Definir colores amigables para personas con daltonismo
color1 = [0, 0.45, 0.70]; % Azul
color2 = [0.85, 0.33, 0.10]; % Rojo anaranjado
color3 = [0.93, 0.69, 0.13]; % Amarillo
color4 = [0.49, 0.18, 0.56]; % Púrpura
color5 = [0.47, 0.67, 0.19]; % Verde
color6 = [0.30, 0.75, 0.93]; % Cian
color7 = [0.64, 0.08, 0.18]; % Marrón
figure;
hold on;

% Obtener los límites del eje x y y para ajustar los parches de color
x_limits = [min(foi2), max(foi2)];
y_limits = [min(cohen_d), max(cohen_d)]; % Ajustar automáticamente al rango actual del eje y

% Colorear la región superior (color2) y la región inferior (color5)
red_patch = fill([x_limits, fliplr(x_limits)], [0, 0, y_limits(2), y_limits(2)], ...
    color1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
green_patch = fill([x_limits, fliplr(x_limits)], [y_limits(1), y_limits(1), 0, 0], ...
    color2, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
line_plot = plot(foi2, cohen_d, 'Color', color4, 'LineWidth', 2);
title('Tamaño del Efecto (Cohen''s d)');
xlabel('Frecuencia (Hz)');
ylabel('Tamaño del Efecto');
% Añadir la leyenda
legend([red_patch, green_patch, line_plot], ...
    {'Favors Motort Fixation', 'Favors Motort Right hand', 'Cohen''s d'}, ...
    'Location', 'Best');

grid on;