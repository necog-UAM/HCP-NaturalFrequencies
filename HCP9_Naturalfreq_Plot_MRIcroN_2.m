function HCP9_Naturalfreq_Plot_MRIcroN_2(cohen_d_struct, name, source, source_forward, stat_name, cutoff_point)
    % Validar que stat_name sea 'T_sum' o 'cohen_D'
    valid_stat_names = {'T_sum', 'cohen_D', 'cohen_D_min', 'cohen_D_max','delta_M', 'delta_min', 'delta_max'};
    if ~ismember(stat_name, valid_stat_names)
        error('stat_name must be either ''T_sum'' or ''cohen_D_min''or''cohen_D_max''or''delta_min''or''delta_max''.');
    end
    
    min_names = {'cohen_D_min', 'delta_min'};
    max_names = {'cohen_D_max', 'delta_max', 'cohen_D', 'delta_M'};
    % Extract the conditions being compared from the name
    reference_condition = strsplit(name, '_vs_');
    reference_condition = reference_condition{1}; % Take the first part

    % Find the voxels inside the source
    voxel_inside = find(source.inside == 1);

    % Create a copy of source and assign Cohen's d values to avg.pow
    source2 = source;
    %source2.avg.pow(voxel_inside) = cohen_d_struct.(name).(stat_name);  % Use name as field in cohen_d_struct
    cohen_d2 = zeros(size(cohen_d_struct.(name).(stat_name)));
    if ismember(stat_name, max_names)
        cohen_d2(cohen_d_struct.(name).(stat_name) >= cutoff_point) = cohen_d_struct.(name).(stat_name)(cohen_d_struct.(name).(stat_name) >= cutoff_point);
    elseif ismember(stat_name, min_names)
        cohen_d2(cohen_d_struct.(name).(stat_name) <= cutoff_point) = cohen_d_struct.(name).(stat_name)(cohen_d_struct.(name).(stat_name) <= cutoff_point);
    end
    source2.avg.pow(voxel_inside) = cohen_d2;
    % Reset the avg.mom and set time to 1
    source2.avg.mom = cell(length(source2.inside), 1);

    % Interpolation configuration
    cfg = [];
    cfg.parameter = 'pow';
    cfg.downsample = 2;
    cfg.interpmethod = 'nearest';
    source_interp = ft_sourceinterpolate(cfg, source2, source_forward.mri);

    % Determine the filename suffix based on the presence of '_pow' in the field name
    if contains(name, '_pow')
        suffix = '_pow';
    else
        suffix = '';
    end

    % Write the output as NIfTI file using the provided name
    cfg = [];
    cfg.filetype = 'nifti';
    cfg.parameter = 'pow';
    cfg.filename = [name, '_', stat_name, suffix];  % Use name to save the file
    ft_sourcewrite(cfg, source_interp);

    fprintf('NIfTI file saved as: %s_%s%s_sig.nii\n', name, stat_name, suffix);
end