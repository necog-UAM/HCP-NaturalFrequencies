function [normality_p, T_sum, cohen_d_M, delta_M] = HCP9_Naturalfreq_cohen_d_max(fnat_struct, condition1, condition2, voxel, foi1, foi1x, foi2, neighbor_matrix, substsk, measure_type)
    % HCP9_Naturalfreq_cohen_d genera gráficos comparativos para dos condiciones,
    % calculando el promedio del voxel y sus vecinos por sujeto y luego promediando entre sujetos.
    %
    % Inputs:
    %   - fnat_struct: la estructura que contiene los datos de frecuencia natural
    %   - condition1: string para la primera condición (ej. 'Motort_Left_hand')
    %   - condition2: string para la segunda condición (ej. 'Motort_Right_hand')
    %   - voxel: número del voxel a analizar
    %   - foi1: rango de frecuencia para histograma
    %   - foi2: rango de frecuencia interpolada para el gráfico
    %   - neighbor_matrix: matriz de vecinos para los voxeles
    %   - substsk: lista de sujetos disponibles (cell array 1x4)
    %   - measure_type: tipo de medidas ('rep_measures' o 'indep_measures')

    if ~ismember(measure_type, {'rep_measures', 'indep_measures'})
        error('measure_type debe ser ''rep_measures'' o ''indep_measures''.');
    end

    % Obtener el número de sujetos en cada condición
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

    % Selección de sujetos comunes (solo para medidas repetidas)
    if strcmp(measure_type, 'rep_measures')
        % Identificar sujetos comunes en ambas condiciones
        % Asumiendo que substsk es un cell array donde cada cell corresponde a una tarea
        tasks = {'Restin', 'StoryM', 'Wrkmem', 'Motort'};
        task_idx1 = find(contains(tasks, strtok(condition1, '_'))); % Encontrar índice de la tarea 1
        task_idx2 = find(contains(tasks, strtok(condition2, '_'))); % Encontrar índice de la tarea 2

        if isempty(task_idx1) || isempty(task_idx2)
            error('No se encontraron tareas correspondientes a las condiciones.');
        end

        % Intersección de sujetos comunes
        common_subjects = intersect(substsk{task_idx1}, substsk{task_idx2});
        if isempty(common_subjects)
            warning('No se encontraron sujetos comunes entre las condiciones.');
            normality_p = [];
            T_sum = [];
            cohen_d_M = [];
            delta_M = [];
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
        % No es necesario filtrar sujetos para medidas independientes
        % Asegúrate de que no haya superposición entre sujetos en ambas condiciones
        task_idx1 = find(contains(tasks, strtok(condition1, '_'))); % Encontrar índice de la tarea 1
        task_idx2 = find(contains(tasks, strtok(condition2, '_'))); % Encontrar índice de la tarea 2

        if any(ismember(substsk{task_idx1}, substsk{task_idx2}))
            warning('Se detectó solapamiento de sujetos entre condiciones. Considera usar ''rep_measures''.');
        end
    end

    % Normality check
    %normality_p = checkNormality(subj_condition1, subj_condition2);

    % Cálculo T Student y tamaño del efecto (Cohen's d)
    if strcmp(measure_type, 'rep_measures')
        normality_p = checkNormality(subj_condition1, subj_condition2, 'rep_measures');
        % Para medidas repetidas
        [h, p, ci, stats] = ttest(subj_condition1, subj_condition2); % Paired t-test
        T_sum = nansum(abs(stats.tstat)); % Suma absoluta de estadísticas t

        % Cohen's d para medidas repetidas
        cohen_d_rm = cohensDRM(subj_condition1, subj_condition2);
        cohen_d_M = nanmax(abs(cohen_d_rm)); % Máximo del valor absoluto de Cohen's d

        % Cliff's delta para medidas repetidas
        delta = cliffsDeltaMatrix(subj_condition1, subj_condition2, 'rep_measures');
        delta_M = nanmax(abs(delta)); % Máximo del valor absoluto de Cliff's delta

    elseif strcmp(measure_type, 'indep_measures')
        normality_p = checkNormality(subj_condition1, subj_condition2, 'indep_measures');
        % Para medidas independientes
        [h, p, ci, stats] = ttest2(subj_condition1, subj_condition2); % Independent t-test
        T_sum = nansum(abs(stats.tstat)); % Suma absoluta de estadísticas t

        % Cohen's d para medidas independientes
        sd_pooled = sqrt((var(subj_condition1, 0, 1) + var(subj_condition2, 0, 1)) / 2); % Varianza agrupada
        cohen_d = (mean(subj_condition1, 1) - mean(subj_condition2, 1)) ./ sd_pooled; % Effect size
        cohen_d_M = nanmax(abs(cohen_d)); % Máximo del valor absoluto de Cohen's d

        % Cliff's delta para medidas independientes
        delta = cliffsDeltaMatrix(subj_condition1, subj_condition2, 'indep_measures');
        delta_M = nanmax(abs(delta)); % Máximo del valor absoluto de Cliff's delta
    end
end