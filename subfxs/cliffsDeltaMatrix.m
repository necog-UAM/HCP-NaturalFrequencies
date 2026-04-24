function deltas = cliffsDeltaMatrix(subj_condition1, subj_condition2, measure_type)
    % Verificar que el tipo de medida sea válido
    if ~ismember(measure_type, {'rep_measures', 'indep_measures'})
        error('measure_type debe ser ''rep_measures'' o ''indep_measures''.');
    end

    % Get the number of subjects and frequencies
    [numSubjects, numFrequencies] = size(subj_condition1);

    % Inicializar un vector para almacenar Cliff's Delta para cada frecuencia
    deltas = zeros(1, numFrequencies);

    % Calcular Cliff's Delta según el tipo de medida
    for freq = 1:numFrequencies
        % Extraer los datos para la frecuencia actual
        data1 = subj_condition1(:, freq); % Condición 1 para esta frecuencia
        data2 = subj_condition2(:, freq); % Condición 2 para esta frecuencia

        if strcmp(measure_type, 'rep_measures')
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
        elseif strcmp(measure_type, 'indep_measures')
            % Obtener el número de observaciones en cada grupo
            n1 = length(data1);
            n2 = length(data2);

            % Conteo de diferencias
            count = 0;
            total = n1 * n2;

            % Comparar todos los pares posibles entre los dos grupos
            for i = 1:n1
                for j = 1:n2
                    if data1(i) > data2(j)
                        count = count + 1; % Contar casos donde Condición 1 > Condición 2
                    elseif data1(i) < data2(j)
                        count = count - 1; % Contar casos donde Condición 1 < Condición 2
                    end
                end
            end

            % Calcular Cliff's Delta
            deltas(freq) = count / total;
        end
    end
end