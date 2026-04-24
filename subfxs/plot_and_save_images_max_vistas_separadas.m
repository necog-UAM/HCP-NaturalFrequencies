function plot_and_save_images_max_vistas_separadas(comparison_name, struct_data, stat_name, source, source_forward, cutoff_point_D, condition_dir)
    stats = {'cohen_D'};
    for stat = stats
        stat_name = stat{1};
        if isfield(struct_data.(comparison_name), stat_name)

            voxel_inside = find(source.inside == 1);
            source2 = source;
            pow_values = struct_data.(comparison_name).(stat_name);

            % Poner como NaN los valores que no superan el punto de corte
            if strcmp(stat_name, 'cohen_D')
                pow_values(pow_values < cutoff_point_D) = NaN;
            end

            % --- 1. Eliminar NaNs para calcular el máx real (opcional) ---
            valid_values = pow_values(~isnan(pow_values));

            if isempty(valid_values)
                warning('No hay valores válidos para visualizar en %s', comparison_name);
                continue;
            end

            % --- 2. Calcular el máximo y redondearlo a 1 decimal ---
            max_val_raw = max(valid_values);
            max_val = round(max_val_raw, 2, 'significant');  % Redondea a 1 decimal significativo
            
            source2.avg.pow(voxel_inside) = pow_values;
            source2.avg.mom = cell(length(source2.inside), 1);

            % Interpolación
            cfg = [];
            cfg.parameter = 'pow';
            cfg.downsample = 2;
            cfg.interpmethod = 'linear';
            source_interp = ft_sourceinterpolate(cfg, source2, source_forward.mri);


            % Parámetros base de visualización
            cfg = [];
            cfg.method = 'surface';
            cfg.funparameter = 'pow';
            cfg.maskparameter = cfg.funparameter;
            % --- 3. Configurar FieldTrip para usar ese rango ---
            cfg.funcolorlim = [cutoff_point_D - 0.2 1.7];            % Desde 0 hasta el máximo redondeado
            cfg.colorbarlimits = cfg.funcolorlim;     % Hacer que el colorbar use esos límites
            %cfg.funcolormap = 'jet_omega_mod';
            cfg.funcolormap   = 'hot_omega_mod';
            cfg.projmethod = 'nearest';
            cfg.opacity = 0.8;
            cfg.camlight = 'no';
            %cfg.funcolorlim = [0 2.05];      % Aquí defines el rango de colores
            %cfg.colorbarlimits = cfg.funcolorlim;  % Asegura que el colorbar use esos límites
            cfg.colorbar = 'yes';            % Muestra el colorbar
            


            % ---- LEFT MEDIAL ----
            cfg.surffile = 'surface_pial_left.mat';
            cfg.surfinflated = 'surface_inflated_left_caret_white.mat';

            f1 = figure('Visible', 'off', 'Color', [1 1 1]);
            ft_sourceplot(cfg, source_interp), view([-90 0]), camlight('left');
            % --- Ajustar tamaño del colorbar ---
            cb = colorbar;  % Captura el colorbar recién creado
            if ~isempty(cb)
                % Definir los ticks que queremos mostrar
                tick_values = [0, cutoff_point_D, 1.7];

                % Asignar ticks y etiquetas con 1 decimal
                cb.Ticks = tick_values;
                cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x), tick_values, 'UniformOutput', false);
                cb.FontSize = 20;
            end
            filename = sprintf('%s_%s_left_lateral.tiff', comparison_name,  stat_name);
            %exportgraphics(f1, fullfile(condition_dir, filename), 'BackgroundColor', 'none', 'ContentType', 'image');
            print(f1, fullfile(condition_dir, filename), '-dtiff' ,'-r300');
            close(f1)

            % ---- LEFT LATERAL ----
            f2 = figure('Visible', 'off', 'Color', [1 1 1]);
            ft_sourceplot(cfg, source_interp), view([90 0]), camlight('left');
            % --- Ajustar tamaño del colorbar ---
            cb = colorbar;  % Captura el colorbar recién creado
            if ~isempty(cb)
                % Definir los ticks que queremos mostrar
                tick_values = [0, cutoff_point_D, max_val];

                % Asignar ticks y etiquetas con 1 decimal
                cb.Ticks = tick_values;
                cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x), tick_values, 'UniformOutput', false);
                cb.FontSize = 20;
            end
            filename = sprintf('%s_%s_left_medial.tiff', comparison_name,  stat_name);
            %exportgraphics(f2, fullfile(condition_dir, filename), 'BackgroundColor', 'none', 'ContentType', 'image');
            print(f2, fullfile(condition_dir, filename), '-dtiff' ,'-r300');
            close(f2)

            % ---- RIGHT LATERAL ----
            cfg.surffile = 'surface_pial_right.mat';
            cfg.surfinflated = 'surface_inflated_right_caret_white.mat';

            f3 = figure('Visible', 'off', 'Color', [1 1 1]);
            ft_sourceplot(cfg, source_interp), view([90 0]), camlight('right');
            % --- Ajustar tamaño del colorbar ---
            cb = colorbar;  % Captura el colorbar recién creado
            if ~isempty(cb)
                % Definir los ticks que queremos mostrar
                tick_values = [0, cutoff_point_D, max_val];

                % Asignar ticks y etiquetas con 1 decimal
                cb.Ticks = tick_values;
                cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x), tick_values, 'UniformOutput', false);
                cb.FontSize = 20;
            end
            filename = sprintf('%s_%s_right_lateral.tiff', comparison_name,  stat_name);
            %exportgraphics(f3, fullfile(condition_dir, filename), 'BackgroundColor', 'none', 'ContentType', 'image');
            print(f3, fullfile(condition_dir, filename), '-dtiff' ,'-r300');
            close(f3)

            % ---- RIGHT MEDIAL ----
            f4 = figure('Visible', 'off', 'Color', [1 1 1]);
            ft_sourceplot(cfg, source_interp), view([-90 0]), camlight('right');
            % --- Ajustar tamaño del colorbar ---
            cb = colorbar;  % Captura el colorbar recién creado
            if ~isempty(cb)
                % Definir los ticks que queremos mostrar
                tick_values = [0, cutoff_point_D, max_val];

                % Asignar ticks y etiquetas con 1 decimal
                cb.Ticks = tick_values;
                cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x), tick_values, 'UniformOutput', false);
                cb.FontSize = 20;
            end
            filename = sprintf('%s_%s_right_medial.tiff', comparison_name,  stat_name);
            %exportgraphics(f4, fullfile(condition_dir, filename), 'BackgroundColor', 'none', 'ContentType', 'image');
            print(f4, fullfile(condition_dir, filename), '-dtiff' ,'-r300');
            close(f4)

        else
            fprintf('Campo %s no encontrado en %s\n', stat_name, comparison_name);
        end
    end
end

% function plot_and_save_images_max_vistas_separadas(comparison_name, struct_data, stat_name, source, source_forward, cutoff_point_D, condition_dir)
%     stats = {'cohen_D'};
%     for stat = stats
%         stat_name = stat{1};
%         if isfield(struct_data.(comparison_name), stat_name)
% 
%             voxel_inside = find(source.inside == 1);
%             source2 = source;
%             pow_values = struct_data.(comparison_name).(stat_name);
% 
%             % --- 1. Eliminar NaNs para calcular el máx real ---
%             valid_values = pow_values(~isnan(pow_values));
%             if isempty(valid_values)
%                 warning('No hay valores válidos para visualizar en %s', comparison_name);
%                 continue;
%             end
% 
%             % --- 2. Calcular el máximo y redondearlo a 1 decimal ---
%             max_val_raw = max(valid_values);
%             max_val = round(max_val_raw, 2, 'significant');  % Redondea a 1 decimal significativo
% 
%             % --- 3. Reemplazar valores < cutoff por un valor fuera del rango visible (gris) ---
%             pow_values_clean = pow_values;
%             pow_values_clean(pow_values_clean < cutoff_point_D) = NaN;
%             %pow_values_clean(isnan(pow_values_clean)) = cutoff_point_D - 1e-6;
% 
%             source2.avg.pow(voxel_inside) = pow_values_clean;
%             source2.avg.mom = cell(length(source2.inside), 1);
% 
%             % Interpolación
%             cfg = [];
%             cfg.parameter = 'pow';
%             cfg.downsample = 2;
%             cfg.interpmethod = 'linear';
%             source_interp = ft_sourceinterpolate(cfg, source2, source_forward.mri);
% 
%             % --- 4. Crear colormap personalizado con 'hot_omega_mod' ---
%             num_colors = 64;
% 
%             % Cargar el colormap base y reescalarlo a num_colors si es necesario
%             base_colormap = hot_omega_mod(num_colors);  % Tu colormap personalizado
% 
%             gray_color = [0.7 0.7 0.7];                 % Gris claro para valores no significativos
%             custom_colormap = [repmat(gray_color, 1, 1); base_colormap];
% 
%             % --- 5. Configurar FieldTrip ---
%             cfg = [];
%             cfg.method = 'surface';
%             cfg.funparameter = 'pow';
%             cfg.maskparameter = cfg.funparameter;
%             cfg.funcolormap = custom_colormap;             % Pasamos directamente la matriz Nx3
%             cfg.funcolorlim = [cutoff_point_D - 0.1 max_val];    % Desde cutoff hasta el máximo redondeado
%             cfg.colorbarlimits = cfg.funcolorlim;          % Hacer que el colorbar use este rango
%             cfg.projmethod = 'nearest';
%             cfg.opacity = 0.8;
%             cfg.camlight = 'no';
%             cfg.colorbar = 'yes';
% 
%             % ---- LEFT MEDIAL ----
%             cfg.surffile = 'surface_pial_left.mat';
%             cfg.surfinflated = 'surface_inflated_left_caret_white.mat';
% 
%             f1 = figure('Visible', 'off', 'Color', [1 1 1]);
%             ft_sourceplot(cfg, source_interp); view([-90 0]); camlight('left');
%             configure_colorbar(f1, cutoff_point_D, max_val);
%             filename = sprintf('%s_%s_left_lateral.tiff', comparison_name, stat_name);
%             print(f1, fullfile(condition_dir, filename), '-dtiff', '-r300');
%             close(f1);
% 
%             % ---- LEFT LATERAL ----
%             f2 = figure('Visible', 'off', 'Color', [1 1 1]);
%             ft_sourceplot(cfg, source_interp); view([90 0]); camlight('left');
%             configure_colorbar(f2, cutoff_point_D, max_val);
%             filename = sprintf('%s_%s_left_medial.tiff', comparison_name, stat_name);
%             print(f2, fullfile(condition_dir, filename), '-dtiff', '-r300');
%             close(f2);
% 
%             % ---- RIGHT LATERAL ----
%             cfg.surffile = 'surface_pial_right.mat';
%             cfg.surfinflated = 'surface_inflated_right_caret_white.mat';
% 
%             f3 = figure('Visible', 'off', 'Color', [1 1 1]);
%             ft_sourceplot(cfg, source_interp); view([90 0]); camlight('right');
%             configure_colorbar(f3, cutoff_point_D, max_val);
%             filename = sprintf('%s_%s_right_lateral.tiff', comparison_name, stat_name);
%             print(f3, fullfile(condition_dir, filename), '-dtiff', '-r300');
%             close(f3);
% 
%             % ---- RIGHT MEDIAL ----
%             f4 = figure('Visible', 'off', 'Color', [1 1 1]);
%             ft_sourceplot(cfg, source_interp); view([-90 0]); camlight('right');
%             configure_colorbar(f4, cutoff_point_D, max_val);
%             filename = sprintf('%s_%s_right_medial.tiff', comparison_name, stat_name);
%             print(f4, fullfile(condition_dir, filename), '-dtiff', '-r300');
%             close(f4);
% 
%         else
%             fprintf('Campo %s no encontrado en %s\n', stat_name, comparison_name);
%         end
%     end
% end