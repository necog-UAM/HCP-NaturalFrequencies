function [cohen_d_M, delta_M, T_sum, p_values] = HCP9_Naturalfreq_Plot_Conditions_max_final(fnat_struct, condition1, condition2, voxel, foi1, foi1x, foi2, neighbor_matrix, substsk, plot_flag, measure_type)
    % HCP9_Naturalfreq_Plot_Conditions genera gráficos comparativos para dos condiciones,
    % calculando el promedio del voxel y sus vecinos por sujeto y luego promediando entre sujetos.
    %
    % Inputs:
    %   - fnat_struct: la estructura que contiene los datos de frecuencia natural
    %   - condition1: string para la primera condición (ej. 'Motort_Left_hand')
    %   - condition2: string para la segunda condición (ej. 'Motort_Right_hand')
    %   - voxel: número del voxel a analizar
    %   - foi1x: rango de frecuencia para histograma
    %   - foi2: rango de frecuencia interpolada para el gráfico
    %   - neighbor_matrix: matriz de vecinos para los voxeles
    %   - measure_type: tipo de medidas ('rep_measures' o 'indep_measures')

    % Verificar que el tipo de medida sea válido
    if ~ismember(measure_type, {'rep_measures', 'indep_measures'})
        error('measure_type debe ser ''rep_measures'' o ''indep_measures''.');
    end

    % Si no se pasan argumentos opcionales, asignar valores predeterminados
    if nargin < 10
        plot_flag = true; % Por defecto se activan los plots
    end

    % Obtener el número de sujetos
    Nsub1 = size(fnat_struct.(condition1), 1); 
    Nsub2 = size(fnat_struct.(condition2), 1); 

    % Pre-alocación de matrices para almacenar los promedios por sujeto
    subj_condition1 = zeros(Nsub1, length(foi2));
    subj_condition2 = zeros(Nsub2, length(foi2));

    % Calcular el promedio para cada sujeto y condición usando vecinos
    for subj = 1:Nsub1
        task_condition_voxel1 = HCP8_VoxelNeighbors_Histogram_flag(fnat_struct, condition1, subj, voxel, foi1, foi1x, foi2, neighbor_matrix, substsk, false);
        subj_condition1(subj, :) = mean(task_condition_voxel1, 1); % Promedio por sujeto
    end

    for subj = 1:Nsub2
        task_condition_voxel2 = HCP8_VoxelNeighbors_Histogram_flag(fnat_struct, condition2, subj, voxel, foi1, foi1x, foi2, neighbor_matrix, substsk, false);
        subj_condition2(subj, :) = mean(task_condition_voxel2, 1); % Promedio por sujeto
    end

    % Filtrar solo sujetos comunes (solo para medidas repetidas)
    if strcmp(measure_type, 'rep_measures')
        % Identificar sujetos comunes en ambas condiciones
        tasks = {'Restin', 'StoryM', 'Wrkmem', 'Motort'};
        task_idx1 = find(contains(tasks, strtok(condition1, '_'))); % Encontrar índice de la tarea 1
        task_idx2 = find(contains(tasks, strtok(condition2, '_'))); % Encontrar índice de la tarea 2

        if isempty(task_idx1) || isempty(task_idx2)
            error('No se encontraron tareas correspondientes a las condiciones.');
        end

        common_subjects = intersect(substsk{task_idx1}, substsk{task_idx2});
        if isempty(common_subjects)
            warning('No se encontraron sujetos comunes entre las condiciones.');
            cohen_d_M = [];
            delta_M = [];
            T_sum = [];
            p_values = [];
            return;
        end

        % Filtrar datos solo para sujetos comunes
        subj_condition1_common = zeros(length(common_subjects), size(subj_condition1, 2));
        subj_condition2_common = zeros(length(common_subjects), size(subj_condition2, 2));

        for i = 1:length(common_subjects)
            subj_idx1 = find(ismember(substsk{task_idx1}, common_subjects(i)));
            subj_idx2 = find(ismember(substsk{task_idx2}, common_subjects(i)));

            if ~isempty(subj_idx1) && ~isempty(subj_idx2)
                subj_condition1_common(i, :) = subj_condition1(subj_idx1, :);
                subj_condition2_common(i, :) = subj_condition2(subj_idx2, :);
            else
                warning('Sujeto común no encontrado en ambos conjuntos de datos.');
            end
        end

        % Reemplazar las matrices originales con las filtradas
        subj_condition1 = subj_condition1_common;
        subj_condition2 = subj_condition2_common;

    elseif strcmp(measure_type, 'indep_measures')
        % Para medidas independientes, no es necesario filtrar sujetos
        if any(ismember(substsk{find(contains(tasks, strtok(condition1, '_')))}, ...
               substsk{find(contains(tasks, strtok(condition2, '_')))}))
            warning('Se detectó solapamiento de sujetos entre condiciones. Considera usar ''rep_measures''.');
        end
    end

    % Calcular la media y la desviación estándar entre sujetos para cada condición
    mean_condition1 = mean(subj_condition1, 1);
    std_condition1 = std(subj_condition1, 0, 1);

    mean_condition2 = mean(subj_condition2, 1);
    std_condition2 = std(subj_condition2, 0, 1);

    % Calcular el Error de la Media Estándar (SEM) y los Intervalos de Confianza (IC)
    sem_condition1 = std_condition1 / sqrt(size(subj_condition1, 1));
    sem_condition2 = std_condition2 / sqrt(size(subj_condition2, 1));

    cil1 = mean_condition1 - 1.96 * sem_condition1;
    cih1 = mean_condition1 + 1.96 * sem_condition1;

    cil2 = mean_condition2 - 1.96 * sem_condition2;
    cih2 = mean_condition2 + 1.96 * sem_condition2;

    % Realizar prueba t según el tipo de medida
    if strcmp(measure_type, 'rep_measures')
        [h, p_values, ci, stats] = ttest(subj_condition1, subj_condition2); % Paired t-test
        % Cálculo del tamaño del efecto (Cohen's d)
        cohen_d = cohensDRM(subj_condition1, subj_condition2);
        delta = cliffsDeltaMatrix(subj_condition1, subj_condition2, 'rep_measures'); % Cliff's Delta para medidas repetidas
    elseif strcmp(measure_type, 'indep_measures')
        [h, p_values, ci, stats] = ttest2(subj_condition1, subj_condition2); % Independent t-test
        % Cohen's d para medidas independientes
        sd_pooled = sqrt((var(subj_condition1, 0, 1) + var(subj_condition2, 0, 1)) / 2); % Varianza agrupada
        cohen_d = (mean(subj_condition1, 1) - mean(subj_condition2, 1)) ./ sd_pooled; % Effect size
        delta = cliffsMatrix(subj_condition1, subj_condition2, 'indep_measures'); % Cliff's Delta para medidas independientes

    end

    cohen_d(isnan(cohen_d)) = 0; % Reemplaza NaN por 0
    cohen_d_M = nanmean(abs(cohen_d)); % Media del valor absoluto de Cohen's d
    T_sum = nansum(abs(stats.tstat)); % Suma absoluta de estadísticas t
    delta_M = nanmean(abs(delta)); % Media del valor absoluto de Cliff's Delta

    % Código nuevo para calcular porcentaje de solapamiento
    overlap_percentage = 100 * (1 - normcdf(abs(cohen_d) / 2));

    % Ploteo de resultados
    if plot_flag
        % Definir colores amigables para personas con daltonismo
        color1 = [0, 0.45, 0.70]; % Azul
        color2 = [0.85, 0.33, 0.10]; % Rojo anaranjado
        color3 = [0.93, 0.69, 0.13]; % Amarillo
        color4 = [0.49, 0.18, 0.56]; % Púrpura
        color5 = [0.47, 0.67, 0.19]; % Verde
        color6 = [0.30, 0.75, 0.93]; % Cian
        color7 = [0.64, 0.08, 0.18]; % Marrón

        % Graficar resultados
        figure;
        hold on;

        % Área sombreada para la condición 1 con transparencia
        fill1 = fill([foi2, fliplr(foi2)], [cil1, fliplr(cih1)], color1, 'EdgeColor', 'none', 'FaceAlpha', 0.3);

        % Líneas de intervalo de confianza
        %plot(foi2, cil1, '--', 'Color', color1, 'LineWidth', 1.5);
        %plot(foi2, cih1, '--', 'Color', color1, 'LineWidth', 1.5);

        % Línea de media
        mean1 = plot(foi2, mean_condition1, 'Color', color1, 'LineWidth', 3);

        fill2 = fill([foi2, fliplr(foi2)], [cil2, fliplr(cih2)], color2, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        %plot(foi2, cil2, '--', 'Color', color2, 'LineWidth', 1.5);
        %plot(foi2, cih2, '--', 'Color', color2, 'LineWidth', 1.5);
        mean2 = plot(foi2, mean_condition2, 'Color', color2, 'LineWidth', 3);

        % Personalización del gráfico
        %title(['Comparación entre ', strrep(condition1, '_', ' '), ' y ', strrep(condition2, '_', ' '), ' para el voxel ', num2str(voxel)]);
        xlabel('Frequency (Hz)', 'FontSize', 14); % Aumentar tamaño de fuente en eje X
        % ylabel eliminado
        ax = gca;
        ax.YTick = 0:0.1:1; % Mostrar ticks del eje Y cada 0.1
        ax.XTick = 0:10:80;
        ax.XLim = [0 80];
        ax.YLim = [0 0.8];
        ax.XGrid = 'off'; % Turn off vertical grid lines
        ax.YGrid = 'off';% Turn off horizontal grid lines
        ax.XAxis.FontSize = 20; % Tamaño de número de eje X
        ax.YAxis.FontSize = 20; % Tamaño de número de eje Y si quieres ajustarlo también
        ax.YAxis.TickLabels = ''; % Ocultar los valores numéricos del eje Y, pero mantener los tics

        % Leyenda solo para las líneas de media
        tasks = {'Restin', 'StoryM', 'Wrkmem', 'Motort'};  % Lista de posibles tasks

        % Función anónima mejorada para limpiar nombres
        cleanCond = @(cond, t) ...
            strrep(...
            strtrim(strrep(cond, ...
            t{find(cellfun(@(x) strncmp(cond, [x '_'], length(x)+1), t), 1)}, '')), ...
            '_', ' ');

        % Limpiamos las etiquetas
        if(strcmp(condition1, 'Restin_Restin'))
            label1 = 'Resting State';
        else
            label1 = cleanCond(condition1, tasks);  % Resultado: 'Fixation'
        end
        label2 = cleanCond(condition2, tasks);  % Resultado: 'Encoding Task'
        legend([mean1, mean2], ...
            {label1, label2}, ...
            'Location', 'northeast', 'Interpreter', 'none', 'FontSize', 20);

        hold off;


        % Gráfico de p-valores
%         figure;
%         plot(foi2, p_values, 'Color', color4, 'LineWidth', 1.5);
%         title('P-valores del T-test');
%         xlabel('Frecuencia (Hz)');
%         ylabel('P-value');
%         
%         % Añadir líneas horizontales para umbrales de significancia
%         hold on;
%         yline(0.05, '--', 'Color', color3, 'LineWidth', 1.5); % Línea amarilla discontinua en 0.05
%         yline(0.01, '--', 'Color', color7, 'LineWidth', 1.5); % Línea marrón discontinua en 0.01
%         hold off;

        % Gráfico tamaño del efecto (Cohen's d y Cliff's Delta)
        figure;
        hold on;
        
        % Obtener los límites del eje x y y para ajustar los parches de color
        x_limits = [min(foi2), max(foi2)];
        y1_limits = [min(cohen_d), max(cohen_d)]; % Límites para Cohen's d
        y2_limits = [min(delta), max(delta)];    % Límites para Cliff's Delta
        
        % Colorear la región superior (color1) y la región inferior (color2) para Cohen's d
        red_patch = fill([x_limits, fliplr(x_limits)], [0, 0, y1_limits(2), y1_limits(2)], ...
            color1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        green_patch = fill([x_limits, fliplr(x_limits)], [y1_limits(1), y1_limits(1), 0, 0], ...
            color2, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        % Plot Cohen's d
        line_cohen = plot(foi2, cohen_d, 'Color', color4, 'LineWidth', 2);
        
        % Plot Cliff's Delta (e.g., in a different color or style)
        line_delta = plot(foi2, delta, 'Color', color3, 'LineStyle', '--', 'LineWidth', 2);
        
        % Añadir título y etiquetas
        title('Tamaño del Efecto (Cohen''s d y Cliff''s Delta)');
        xlabel('Frecuencia (Hz)');
        ylabel('Tamaño del Efecto');
        
        % Ajustar el rango del eje y para incluir ambos efectos
        y_combined_limits = [min([y1_limits(1), y2_limits(1)]), max([y1_limits(2), y2_limits(2)])];
        ylim(y_combined_limits);
        
        % Añadir la leyenda
        legend([red_patch, green_patch, line_cohen, line_delta], ...
            {['Favors ' label1], ['Favors ' label2], 'Cohen''s d', 'Cliff''s Delta'}, ...
            'Location', 'Best');

        grid on;
    end
end

% source.pos(find(source.inside == 1),:)