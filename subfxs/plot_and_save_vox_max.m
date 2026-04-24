function plot_and_save_vox_max(comparison_name, struct_data, stat_name, source, source_forward, cutoff_point_D, voxel_index, condition_dir)

    stats = {'cohen_D'};
    for stat = stats
        stat_name = stat{1};
        if isfield(struct_data.(comparison_name), stat_name)

            voxel_inside = find(source.inside == 1);

            % Buscar la posición del voxel en los índices internos
            idx_inside = voxel_inside(voxel_index);
            if isempty(idx_inside)
                fprintf('El voxel %d no está dentro de source.inside\n', voxel_index);
                return
            end

            % Crear un nuevo source con todos los valores en 0, excepto el voxel deseado en 1
            source2 = source;
            source2.avg.pow = zeros(length(source.inside),1);
            source2.avg.pow(voxel_inside(voxel_index)) = 1;
            source2.avg.mom = cell(length(source.inside), 1);

            % Interpolación (para proyectar a superficie)
            cfg = [];
            cfg.parameter = 'pow';
            cfg.downsample = 2;
            cfg.interpmethod = 'nearest';  % Importante para mantener el valor discreto
            source_interp = ft_sourceinterpolate(cfg, source2, source_forward.mri);

            % Visualización
            figure('WindowState', 'maximized', 'Color', [1 1 1]);

            cfg = [];
            cfg.figure = 'gca';
            cfg.method = 'surface';
            cfg.funparameter = 'pow';
            cfg.maskparameter = cfg.funparameter;
            cfg.funcolorlim = [0 1];  % Escala fija: 0 a 1
            cfg.funcolormap = 'hot_omega_mod';
            cfg.projmethod = 'nearest';
            cfg.opacitymap = 'rampup';
            cfg.camlight = 'no';
            cfg.colorbar = 'yes';

            % Superficie izquierda
            cfg.surffile = 'surface_pial_left.mat';
            cfg.surfinflated = 'surface_inflated_left_caret_white.mat';
            subplot(2, 2, 1), ft_sourceplot(cfg, source_interp), view([-90 0]), camlight('left');
            subplot(2, 2, 3), ft_sourceplot(cfg, source_interp), view([90 0]), camlight('left');

            % Superficie derecha
            cfg.surffile = 'surface_pial_right.mat';
            cfg.surfinflated = 'surface_inflated_right_caret_white.mat';
            subplot(2, 2, 2), ft_sourceplot(cfg, source_interp), view([90 0]), camlight('right');
            subplot(2, 2, 4), ft_sourceplot(cfg, source_interp), view([-90 0]), camlight('right');

            % Guardar imagen
            filename = sprintf('%s_%s_voxel_%d.tiff', comparison_name, stat_name, voxel_index);
            cd(condition_dir);
            %saveas(gcf, filename);
            print(gcf, filename, '-dtiff' ,'-r300');
            close(gcf)

        else
            fprintf('Campo %s no encontrado en %s\n', stat_name, comparison_name);
        end
    end
end
