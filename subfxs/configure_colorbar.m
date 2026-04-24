function configure_colorbar(fig, cutoff_point_D, max_val)
    % Cambia el tamaño del colorbar y sus etiquetas
    figure(fig);
    cb = colorbar;

    if ~isempty(cb)
        % Definir los ticks que queremos mostrar
        tick_values = [0, cutoff_point_D, max_val];

        % Asignar ticks y etiquetas con 1 decimal
        cb.Ticks = tick_values;
        cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x), tick_values, 'UniformOutput', false);
        cb.FontSize = 20;
        %cb.Label.String = 'Cohen''s d';
        %cb.Label.FontSize = 22;

        % Añadir texto "N.S." en la parte gris
        x_pos = cb.Position(1);                        % Posición x del colorbar
        y_pos = cb.Position(2) + 0.5 * cb.Position(4);  % Centrado verticalmente en la parte gris

%         text(x_pos, y_pos, 'N.S.', ...
%             'Units', 'normalized', ...
%             'HorizontalAlignment', 'center', ...
%             'VerticalAlignment', 'middle', ...
%             'FontSize', 12, ...
%             'Color', 'k');  % Color negro para mayor legibilidad
    end
end