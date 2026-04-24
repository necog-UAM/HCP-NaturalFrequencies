function [pks, locs, locmx] = HCP9_Naturalfreq_Plot_MRIcroN_LocalMax(cohen_d_struct, name, source, source_forward, stat_name, cutoff_point)
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

%     dinterp = cohen_d2 * dtempl;
%     dinterp2 = reshape(dinterp, [dim(1), dim(2), dim(3)]);
% 
%     % Encontrar los máximos locales y sus valores
%     [pks, locs, locmx] = findlocalmax(dinterp2, 3);

    % Crear una nueva matriz para almacenar los valores de los máximos locales
    max_values = zeros(size(dinterp2));
    for idx = 1:length(pks)
        max_values(locs(idx,1), locs(idx,2), locs(idx,3)) = pks(idx);
    end

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
