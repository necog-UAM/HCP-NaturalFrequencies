function HCP6_Naturalfreq_groupmap (cfg, substsk, task, condit, p, rootDir) %, sub)
%% This script is part 3 to obtain the brain's natural frequencies at rest from EEG recordings
% It is needed singlesub_fnat of the subjects to obtain the group's map

datapath = p.datapath; 
kmeanspath = p.kmeanspath; 
outpath = p.outpath; 

t = cfg.numtask;
Nk   = cfg.Nk;
Nsub = cfg.Nsub;
Nvox = cfg.Nvox;

task    = cfg.task;
numtask = cfg.numtask;
condit  = cfg.condit;
condit_label = cfg.condit_label;
condit_code  = cfg.condit_code;
source = cfg.source.inverse;
source_forward = cfg.source.forward;

cfg2 = cfg;
voxel_inside = find(source.inside==1);
%% Folders
% Obtener la lista de carpetas (sujetos) en el directorio principal
subjects = dir(rootDir);
subjects = subjects([subjects.isdir]); % Filtrar solo las carpetas
subjects = subjects(~ismember({subjects.name}, {'.', '..'})); % Eliminar las carpetas '.' y '..'

% Inicializar una tabla para almacenar los resultados
folders = table('Size', [length(subjects), 2], 'VariableTypes', {'string', 'cell'}, 'VariableNames', {'subject', 'tasks'});

% Recorrer cada carpeta de sujeto
for i = 1:length(subjects)
    subjectName = subjects(i).name;
    subjectPath = fullfile(rootDir, subjectName);
    
    % Obtener la lista de subcarpetas dentro de la carpeta del sujeto
    subfolders = dir(subjectPath);
    subfolders = subfolders([subfolders.isdir]); % Filtrar solo las carpetas
    subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'})); % Eliminar las carpetas '.' y '..'
    
    % Filtrar las subcarpetas que contienen '_' en el nombre
    taskNames = {subfolders.name};
    taskNames = taskNames(contains(taskNames, '_'));
    
    % Añadir la información al dataframe
    folders.subject(i) = string(subjectName);
    folders.tasks{i} = taskNames;
end


%% Load singlesubject fnat.

s=1;
fnatall=[];
for folidx=1:numel(folders.subject)
    % Incluir bucle para task y bucle para condit? --> Fuera de la función
    route = fullfile(datapath, folders.subject{folidx}, ['_' task]); % cambiar 1 por indice de task (t)
    if exist(route)
        cd(route)
        %load singlesub_fnat %singlesub_fnat
        fil = sprintf('load singlesub_fnat_%s_%d_new',task,condit);
        eval(fil)
        fnatall(s,:) = fnat.fnatsig;
        s=s+1;
    end
    % fnatpks{s}=fnat.pks;
    % fnatw{s}=fnat.w;
end

%% Define frequencies for histograms.
f1 = 0;
f2 = 4.65;%3.65;
d  = 0.2; % Why 0.2? -> for histograms
foi1 = exp(f1:d:f2);
foi2 = interp(foi1(1:end-1),10);
% Otra alternativa:
% foix = 0.35:0.05:3.75;
% foix = exp(foix);
% foiint = interp(foi(1:end-1),10);
%% groupmaps
rectp = [];
for i = 1:length(foi2)
    c0 = foi2(i);
    rectp(i,:) = rectangularPulse(c0-c0/4 ,c0+c0/4 ,foi2);      % to isolate individual peaks later
end

pp = parpool(4);                  % parallel pool
pp.IdleTimeout = inf;

natf.A       = NaN(Nvox,1);
natf.mu      = NaN(Nvox,1);
natf.sigma   = NaN(Nvox,1);
natf.rsq     = NaN(Nvox,1);

y = zeros(Nsub,Nvox);
selsubs = 1:Nsub;
% selsubs =find(age>=60);          % to select a subgroup of subjects (only elders for example)

for vx = 1:Nvox
    vx;
    [counts,~] = histcounts(fnatall(selsubs,vx),foi1);
    counts  = interp(counts,10);
    centers = foi2;

    [pks,locs] = findpeaks(counts);
    [~,ii] = sort(pks,'descend');
    if length(ii) >= 4
        c0s = centers(locs(ii(1:4)));         % fit up to 4 peaks 
    else
        c0s = centers(locs);
    end

    % gaussian fit of all candidate peak frequencies
    [A,mu,sigma,rsq] = gausfitc0(c0s,counts,centers,rectp,foi2); %error sometimes, changed gausfitc0
    [~,i2] = sort(A,'descend');

    for ii = 1:length(i2)
        i = i2(ii);
        natf.A(vx,i)       = A(i);
        natf.mu(vx,i)      = mu(i);
        natf.sigma(vx,i)   = sigma(i);
        natf.rsq(vx,i)     = rsq(i);
    
        gausf = A(i) * exp(-(centers-mu(i)).^2 /(2*sigma(i).^2));
        plot(centers,counts,'k'), hold on
        plot(centers,squeeze(gausf),'r')
        set(gca,'XLim',[0 35],'YLim',[0 80])
        hold on
        % pause    

        for s=1:Nsub
            x=fnatall(s,vx);
            y(s,vx) = y(s,vx) + A(i) * exp(-(x-mu(i)).^2 /(2*sigma(i).^2));
        end    
    end
    hold off

    [~,idx] = max(natf.A(vx,:));
    fnatgroup(vx) = natf.mu(vx,idx);         % for each voxel, choose the freq value shared by most of the subjects

end

y2 = y./Nsub;  %%          % divide by Nsub to get values from 0 (no participants show this result) to 1 (all participants show exactly the same result)
y2 = 1-y2;              % invert the data: now 0.95 to 0.99 (only 0.05 of participants show this result), 0.99 to 1 (only 0.01 of participants)


% Plot of group results

source2 = source;
source2.avg.pow(voxel_inside) = log(fnatgroup);
source2.avg.mom = cell(length(source2.inside),1);

cfg=[];
cfg.parameter  = 'pow';
cfg.downsample = 2;
cfg.interpmethod = 'nearest';
source_interp = ft_sourceinterpolate (cfg, source2, source_forward.mri);

figure('WindowState','maximized','Color',[1 1 1]);
% figure

cfg               = [];
cfg.figure        = 'gca';
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.7 3.4];  % same as capilla et al. (2022)    
cfg.funcolormap   = 'jet_omega_mod';
cfg.projmethod    = 'nearest';
cfg.opacity       = 0.8;
cfg.camlight      = 'no';
cfg.colorbar      = 'no';
cfg.surffile     = 'surface_pial_left.mat';
cfg.colorbar = 'yes';            % Muestra el colorbar
cfg.surfinflated  = 'surface_inflated_left_caret_white.mat';
subplot(2,2,1), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('left')
subplot(2,2,3), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('left')

cfg.surffile     = 'surface_pial_right.mat';
cfg.surfinflated  = 'surface_inflated_right_caret_white.mat';
subplot(2,2,2), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('right')
subplot(2,2,4), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('right')

if strmatch(condit_label{t}, 'Restin')
    if ~exist([outpath '\' task '\'], 'dir')
        mkdir([outpath '\' task '\']);
    end
    cd([outpath '\' task '\'])
else
    if ~exist([outpath '\' task '\'  condit_label{numtask}{condit}], 'dir')
        mkdir([outpath '\' task '\'  condit_label{numtask}{condit}]);
    end
    cd([outpath '\' task '\'  condit_label{numtask}{condit}])
end
current_dir = pwd
if ~exist([current_dir '\groupmap'], 'dir')
    mkdir('groupmap');
end

cd('groupmap');
print('-dtiff','-r300',[sprintf('group_fnat_%s_%d.tiff',task,condit)]);
save fnatgroup_new fnatgroup
delete(gcp('nocreate'))
end