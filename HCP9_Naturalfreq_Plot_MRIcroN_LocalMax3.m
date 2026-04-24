function [pks, locs, locmx, voxel_numbers] = HCP9_Naturalfreq_Plot_MRIcroN_LocalMax3(cohen_d_struct, name, source, source_forward, stat_name, cutoff_point, min_distance)
    % Validar que stat_name sea uno de los valores permitidos
    valid_stat_names = {'T_sum', 'cohen_D', 'cohen_D_min', 'cohen_D_max','delta_M', 'delta_min', 'delta_max'};
    if ~ismember(stat_name, valid_stat_names)
        error('stat_name must be either ''T_sum'' or ''cohen_D_min''or''cohen_D_max''or''delta_min''or''delta_max''.');
    end
    
    min_names = {'cohen_D_min', 'delta_min'};
    max_names = {'cohen_D_max', 'delta_max', 'cohen_D', 'delta_M'};
    
    % Extraer la condición de referencia del nombre
    reference_condition = strsplit(name, '_vs_');
    reference_condition = reference_condition{1};
    
    % Encontrar los vóxeles dentro de la fuente
    voxel_inside = find(source.inside == 1);
    
    % Find local maxima in the interpolated data
    for i = 1:3
        pos(:, i) = source.pos(:, i) - min(source.pos(:, i)) + 1;
        dim(i) = max(pos(:, i));
    end

    voxel_inside = find(source.inside == 1);
    Nvox = length(voxel_inside);
    vtempl = zeros(dim);
    dtempl = zeros(Nvox, prod(dim));

    for v = 1:Nvox
        temp = zeros(dim);
        p = pos(voxel_inside(v), :);
        vtempl(p(1), p(2), p(3)) = v;
        temp(p(1), p(2), p(3)) = 1;
        dtempl(v, :) = temp(:);
    end
    
    % Crear una copia de source y asignar valores de Cohen's d a avg.pow
    source2 = source;
    cohen_d2 = cohen_d_struct.(name).(stat_name);

    % Filtrar según el punto de corte
    if ismember(stat_name, max_names)
        cohen_d2(cohen_d2 < cutoff_point) = 0;
        dinterp = cohen_d2 * dtempl;
        dinterp2 = reshape(dinterp, [dim(1), dim(2), dim(3)]);
        % Encontrar los máximos locales y sus valores
        [pks, locs, locmx] = findlocalmax(dinterp2, 3);
    elseif ismember(stat_name, min_names)
        cohen_d2(cohen_d2 > cutoff_point) = 0;
        dinterp = cohen_d2 * dtempl;
        dinterp2 = reshape(dinterp, [dim(1), dim(2), dim(3)]);
        % Encontrar los máximos locales y sus valores
        [pks, locs, locmx] = findlocalmin(dinterp2, 3);
    end

    % Crear una nueva matriz para almacenar los valores de los extremos locales
    max_values = zeros(size(dinterp2));
    for idx = 1:length(pks)
        max_values(locs(idx,1), locs(idx,2), locs(idx,3)) = pks(idx);
    end

    % Ordenar los máximos locales por su valor de pico (pks) en orden descendente
    if ~isempty(pks)
        [sorted_pks, sort_indices] = sort(pks, 'descend'); % Ordenar pks en orden descendente
        locs = locs(sort_indices, :); % Reordenar locs según los índices de ordenamiento
        pks = sorted_pks; % Actualizar pks con los valores ordenados
    end

    % Filtrar los máximos locales que están demasiado cerca
    locs_filtered = [];
    pks_filtered = [];
    if ~isempty(locs) % Verificar que locs no esté vacío
        for i = 1:length(locs)
            keep = true;
            % Solo intentar filtrar si locs_filtered no está vacío
            if ~isempty(locs_filtered)
                for j = 1:size(locs_filtered, 1) % Use size instead of length for clarity
                    dist = sqrt(sum((locs(i, :) - locs_filtered(j, :)).^2)); % Calcular distancia entre máximos locales
                    if dist < min_distance
                        keep = false;
                        break;
                    end
                end
            end
            % Si este máximo local debe ser mantenido, añadirlo a locs_filtered y pks_filtered
            if keep
                locs_filtered = [locs_filtered; locs(i, :)]; % Mantener el máximo local
                pks_filtered = [pks_filtered; pks(i)]; % Mantener el valor correspondiente
            end
        end
    end

    % Actualizar los valores de locs y pks con los máximos locales filtrados
    locs = locs_filtered;
    pks = pks_filtered;

    % Map 3D indices (locs) back to the original 1D voxel indices
    voxel_numbers = [];
    for v = 1:Nvox
        p = pos(voxel_inside(v), :); % Get 3D position of voxel v
        try
            if any(ismember(locs, p, 'rows')) % Check if this voxel's 3D position matches any local extremum
            voxel_numbers = [voxel_numbers, v]; % Store the original voxel index
            end
        end
    end

    
    % Guardar los números de voxel como salida
    fprintf('Voxel numbers corresponding to local extrema: %s\n', mat2str(voxel_numbers));

    source2 = source;
    source2.avg.pow = max_values; % Guardar los valores de los máximos locales
    source2.avg.mom = cell(length(source2.avg.noise), 1);
    source2.time = 1;

    cfg = [];
    cfg.parameter = 'avg.pow';
    cfg.downsample = 2;
    cfg.interpmethod = 'linear';
    source_interp = ft_sourceinterpolate(cfg, source2, source_forward.mri);
    
    % Definir sufijo del archivo
    if contains(name, '_pow')
        suffix = '_pow';
    else
        suffix = '';
    end
    
    % Guardar como NIfTI
    cfg = [];
    cfg.filetype = 'nifti';
    cfg.parameter = 'pow';
    cfg.filename = [name, '_', stat_name, '_localmax', suffix];
    ft_sourcewrite(cfg, source_interp);
    
    fprintf('NIfTI file saved as: %s_%s_localmax%s.nii\n', name, stat_name, suffix);
end
