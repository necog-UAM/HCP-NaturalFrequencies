function [max_t_values, max_cohen_d_values, max_delta_M_values] = HCP9_Naturalfreq_permuted_all_vox_max(fnat_struct, condition1, condition2, foi1, foi1x, foi2, neighbor_matrix, substsk, num_permutations, measure_type)
    % Obtener el número de sujetos y voxeles
    Nsub1 = size(fnat_struct.(condition1), 1); 
    Nsub2 = size(fnat_struct.(condition2), 1); 
    num_voxels = 1925; % Número total de voxeles

    % Validar el tipo de medida
    if ~ismember(measure_type, {'rep_measures', 'indep_measures'})
        error('measure_type debe ser ''rep_measures'' o ''indep_measures''.');
    end

    % Inicializar matrices para almacenar los valores máximos de T, Cohen's d y Cliff's delta
    max_t_values = NaN(num_permutations, 1);
    max_cohen_d_values = NaN(num_permutations, 1);
    max_delta_M_values = NaN(num_permutations, 1);

    % Pre-alocación de matrices para almacenar los promedios por sujeto
    subj_condition1 = zeros(max(Nsub1, Nsub2), length(foi2), num_voxels);
    subj_condition2 = zeros(max(Nsub1, Nsub2), length(foi2), num_voxels);

    % Calcular los datos de cada sujeto y condición fuera del bucle de permutación
    for vox = 1:num_voxels
        for subj = 1:max(Nsub1, Nsub2)
            if subj <= Nsub1 && subj <= Nsub2
                task_condition_voxel1 = HCP8_VoxelNeighbors_Histogram_flag(fnat_struct, condition1, subj, vox, foi1, foi1x, foi2, neighbor_matrix, substsk, false);
                task_condition_voxel2 = HCP8_VoxelNeighbors_Histogram_flag(fnat_struct, condition2, subj, vox, foi1, foi1x, foi2, neighbor_matrix, substsk, false);
                
                % Almacenar los promedios por sujeto
                subj_condition1(subj, :, vox) = mean(task_condition_voxel1, 1);
                subj_condition2(subj, :, vox) = mean(task_condition_voxel2, 1);
            else
                warning('Sujeto no encontrado en ambas condiciones.');
            end
        end
    end

    % Selección de sujetos comunes (solo para medidas repetidas)
    if strcmp(measure_type, 'rep_measures')
        % Identificar sujetos comunes en ambas condiciones
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
            return;
        end

        % Filtrar datos solo para sujetos comunes
        Nsub_common = length(common_subjects);
        subj_condition1_common = zeros(Nsub_common, length(foi2), num_voxels);
        subj_condition2_common = zeros(Nsub_common, length(foi2), num_voxels);

        for i = 1:Nsub_common
            subj_idx1 = find(ismember(substsk{task_idx1}, common_subjects(i)));
            subj_idx2 = find(ismember(substsk{task_idx2}, common_subjects(i)));

            if ~isempty(subj_idx1) && ~isempty(subj_idx2)
                subj_condition1_common(i, :, :) = subj_condition1(subj_idx1, :, :);
                subj_condition2_common(i, :, :) = subj_condition2(subj_idx2, :, :);
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
        if any(ismember(substsk{task_idx1}, substsk{task_idx2}))
            warning('Se detectó solapamiento de sujetos entre condiciones. Considera usar ''rep_measures''.');
        end
    end

    % Realizar permutaciones
    for perm = 1:num_permutations
        % Generar una matriz de permutación aleatoria de 1 y -1
        perm_labels = randi([0, 1], size(subj_condition1, 1), 1) * 2 - 1;

        % Inicializar matrices para almacenar los datos permutados
        perm_subj_condition1 = zeros(size(subj_condition1));
        perm_subj_condition2 = zeros(size(subj_condition2));

        % Asignar los datos del sujeto a una condición u otra según las etiquetas de permutación
        for subj = 1:size(subj_condition1, 1)
            if perm_labels(subj) == 1
                perm_subj_condition1(subj, :, :) = subj_condition1(subj, :, :);
                perm_subj_condition2(subj, :, :) = subj_condition2(subj, :, :);
            else
                perm_subj_condition1(subj, :, :) = subj_condition2(subj, :, :);
                perm_subj_condition2(subj, :, :) = subj_condition1(subj, :, :);
            end
        end

        % Inicializar matrices para almacenar los valores de T, Cohen's d y Cliff's delta para todos los voxeles
        t_values = NaN(num_voxels, 1);
        cohen_d_values = NaN(num_voxels, 1);
        delta_M_values = NaN(num_voxels, 1);

        % Calcular T Student y tamaño del efecto (Cohen's d y Cliff's delta) para cada voxel
        for vox = 1:num_voxels
            data1 = squeeze(perm_subj_condition1(:, :, vox));
            data2 = squeeze(perm_subj_condition2(:, :, vox));

            % Paired t-test para medidas repetidas
            if strcmp(measure_type, 'rep_measures')
                [h, p, ci, stats] = ttest(data1, data2); % Paired t-test
                cohen_d_rm = cohensDRM(data1, data2); % Cohen's d para medidas repetidas
                delta = cliffsDeltaMatrix(data1, data2, 'rep_measures'); % Cliff's delta para medidas repetidas

            % Independent t-test para medidas independientes
            elseif strcmp(measure_type, 'indep_measures')
                [h, p, ci, stats] = ttest2(data1, data2); % Independent t-test
                sd_pooled = sqrt((var(data1, 0, 1) + var(data2, 0, 1)) / 2); % Varianza agrupada
                cohen_d = (mean(data1, 1) - mean(data2, 1)) ./ sd_pooled; % Effect size para medidas independientes
                delta = cliffsDeltaMatrix(data1, data2, 'indep_measures'); % Cliff's delta para medidas independientes
            end

            % Guardar el valor máximo absoluto de T, Cohen's d y Cliff's delta para este voxel
            t_values(vox) = nanmax(abs(stats.tstat));
            cohen_d_values(vox) = nanmax(abs(cohen_d_rm)); % Para medidas repetidas
            delta_M_values(vox) = nanmax(abs(delta));
        end

        % Guardar el valor máximo absoluto de T, Cohen's d y Cliff's delta considerando todos los voxeles
        max_t_values(perm) = max(t_values);
        max_cohen_d_values(perm) = max(cohen_d_values);
        max_delta_M_values(perm) = max(delta_M_values);
    end
end