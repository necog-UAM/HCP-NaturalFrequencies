% function normality_p = checkNormality(subj_condition1, subj_condition2)
%     % Ensure the input matrices have the same dimensions
%     if ~isequal(size(subj_condition1), size(subj_condition2))
%         error('Input matrices must have the same dimensions.');
%     end
% 
%     % Calculate difference scores
%     differences = subj_condition1 - subj_condition2;
% 
%     % Initialize results matrix
%     [numSubjects, numFrequencies] = size(differences);
%     normalityResults = zeros(1, numFrequencies);
% 
%     % Loop over each frequency
%     for freq = 1:numFrequencies
%         % Extract the difference scores for the current frequency
%         diffFreq = differences(:, freq);
% 
%         % Perform Shapiro-Wilk test
%         [~, p] = lillietest(diffFreq); %shapiro(diffFreq);
% 
%         % Store result (1 if normal, 0 if not)
%         normality_p(freq) = p;
%     end
% end

function normality_p = checkNormality(subj_condition1, subj_condition2, measure_type)
    % checkNormality verifica el supuesto de normalidad para análisis estadísticos
    % Inputs:
    %   - subj_condition1: matriz [sujetos x frecuencias] con datos condición 1
    %   - subj_condition2: matriz [sujetos x frecuencias] con datos condición 2
    %   - measure_type: 'rep_measures' o 'indep_measures'
    
    if nargin < 3
        warning('measure_type no especificado. Usando ''rep_measures'' por defecto.');
        measure_type = 'rep_measures';
    end
    
    if ~ismember(measure_type, {'rep_measures', 'indep_measures'})
        error('measure_type debe ser ''rep_measures'' o ''indep_measures''.');
    end
    
    [nsub1, nfreq] = size(subj_condition1);
    [nsub2, ~]     = size(subj_condition2);
    
    normality_p = zeros(1, nfreq);  % Valor-p por frecuencia
    
    if strcmp(measure_type, 'rep_measures')
        % Verificar dimensiones iguales
        if ~isequal(size(subj_condition1), size(subj_condition2))
            error('Para medidas repetidas, matrices deben tener mismas dimensiones.');
        end
        
        differences = subj_condition1 - subj_condition2;
        
        for freq = 1:nfreq
            diffFreq = differences(:, freq);
            [~, p] = lillietest(diffFreq);
            normality_p(freq) = p;
        end
        
    elseif strcmp(measure_type, 'indep_measures')
        % Caso de medidas independientes: probar normalidad en ambos grupos por separado
        normality_p_cond1 = zeros(1, nfreq);
        normality_p_cond2 = zeros(1, nfreq);
        
        for freq = 1:nfreq
            data1 = subj_condition1(:, freq);
            data2 = subj_condition2(:, freq);
            
            % Test de Lilliefors (normalidad) por grupo
            [~, p1] = lillietest(data1);
            [~, p2] = lillietest(data2);
            
            % Combinar resultados (ej: usando mínimo de p-values, u otra estrategia)
            normality_p_cond1(freq) = p1;
            normality_p_cond2(freq) = p2;
            
            % Aquí decidimos qué criterio usar para "global" normalidad
            % Ejemplo: promedio o mín de los p-valores
            normality_p(freq) = min([p1, p2]); % Ser más conservador
        end
    end
end