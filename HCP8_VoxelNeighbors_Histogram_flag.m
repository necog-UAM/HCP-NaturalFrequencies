function task_condition_voxel = HCP8_VoxelNeighbors_Histogram_flag(fnat_struct, condition_label, subj, voxel_idx, foi1, foi1x, foi2, neighbor_matrix, substsk, plot_flag)
%     if nargin < 10
%         plot_flag = true; % Si no se pasa, por defecto se activan los plots
%     end

    % Obtener los vecinos del voxel seleccionado
    voxel_neighbors = find(neighbor_matrix(voxel_idx, :) == 1);

    % Añadir el voxel central a la lista de vecinos
    %voxel_neighbors = [voxel_neighbors]; % Cambiar y eliminar voxel_idx %% voxel_idx, 

    % Inicializar matriz para almacenar los resultados de cada voxel
    task_condition_voxel = zeros(length(voxel_neighbors), length(foi2));

    % Calcular el histograma acumulado para todos los voxeles (central y vecinos)

    total_counts = zeros(1, length(foi2));
    for i = 1:length(voxel_neighbors)
        current_voxel = voxel_neighbors(i);

        % Obtener las cuentas del histograma para el voxel actual
        [counts, ~] = histcounts(fnat_struct.(condition_label)(subj, current_voxel), foi1, 'Normalization', 'probability');
        counts = interp1(foi1x,counts,foi2,'pchip');   % Interpolar para suavizar
       
        % Sumar al total de counts acumulado
        total_counts = total_counts + counts;

        % Guardar los resultados de cada voxel
        task_condition_voxel(i, :) = counts;
    end

    % Encontrar el índice de la tarea usando la etiqueta condition_label
    tasks = {'Restin', 'StoryM', 'Wrkmem', 'Motort'}; % Nombres de las tareas

    % Dividir condition_label para obtener solo la tarea (parte antes del '_')
    task_name = split(condition_label, '_');
    task_name = task_name{1}; % La tarea es la primera parte antes del '_'

    % Encontrar el índice de la tarea en la lista 'tasks'
    task_idx = find(strcmp(tasks, task_name));

    % Extraer el número de sujeto desde substsk usando task_idx
    subject_number = substsk{task_idx}(subj);

    % Cálculo de los valores (sin plot si plot_flag es false)
    if plot_flag
        % Primer gráfico: Histograma suavizado
        % Primer gráfico: Histograma suavizado
        figure;
        hold on;

        plot(foi2, total_counts/length(voxel_neighbors), ...
            'k', 'LineWidth', 2);
        ax = gca;

        % Límites
        ax.XLim = [0 80];
        ax.YLim = [0 0.8];

        % Ticks
        ax.XTick = 0:10:80;
        ax.YTick = 0:0.1:0.8;

        % Solo abajo e izquierda
        ax.Box = 'off';
        ax.TickDir = 'out';         % ticks hacia fuera
        ax.XAxisLocation = 'bottom';
        ax.YAxisLocation = 'left';

        % Quitar ticks superiores y derechos
        ax.XRuler.TickLabelGapOffset = 0;
        ax.YRuler.TickLabelGapOffset = 0;

        % Grosor líneas
        ax.LineWidth = 1.5;

        % Fuente
        ax.FontSize = 20;

        % Sin grid
        grid off
        % Segundo gráfico: Promedio y área sombreada de la desviación típica
        figure;
        hold on;
        mean_voxel = mean(task_condition_voxel, 1);
        std_voxel = std(task_condition_voxel, 0, 1);

        % Graficar el promedio
        plot(foi2, mean_voxel, 'b', 'LineWidth', 2); % Línea azul para el promedio

        % Área sombreada de la desviación estándar
        fill([foi2, fliplr(foi2)], [mean_voxel + std_voxel, fliplr(mean_voxel - std_voxel)], ...
            'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Área azul transparente

        title(['Promedio y desviación estándar para voxel ', num2str(voxel_idx), 'y vecinos del sujeto ', subject_number, ' en ', strrep(condition_label, '_', ' ')]);
        xlabel('Frecuencia');
        ylabel('Cuentas');
        hold off;
    end
end


