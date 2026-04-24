%% BATCH_HCP.m
%-------------------------------------------------------------------------%
% AUTOMATED PIPELINE FOR HCP MEG DATA
% Juan J. Herrera-Morueco, Enrique Stern, Lydia Arana & Almudena Capilla 
% ------------------------------------------------------------------------%
% This script is structured in 2 parts (Part 1: preprocessing and natural
% frequency estimation, Part 2: Statistical analysis & hypothesis testing).
% Each part is organized in sections to facilitate execution.

%% HCP (Part 1)
%-------------------------------------------------------------------------%
% Part 1.
%-------------------------------------------------------------------------%
% This script (1) implements a complete batch-processing workflow for Human
% Connectome Project (HCP) MEG data. It sequentially performs:
% (1) subject validation and path configuration,
% (2) source reconstruction and frequency-domain decomposition across four
% HCP tasks (Rest, Story/Math, Working Memory, Motor),
% (3) K-means clustering of power spectra to identify reproducible spectral 
% profiles,
% (4) computation of individual and group-level natural frequencies maps.
% Finally, it extracts, structure-like objects, and saves the resulting 
% natural frequency dataset into labeled matrices for downstream 
% statistical analysis and visualization (part 2, see below).
%-------------------------------------------------------------------------%
%% Paths and files
% Add required toolboxes and custom script directories to MATLAB's search path
addpath('Z:\Toolbox\fieldtrip-20230118\') %fieldtrip              
addpath(genpath('Z:\HCP\Juancho\')) % Work folder
addpath(genpath('Z:\HCP\Juancho\github\subfxs')) % for function: Omega_neighbors_ly 
% % rmpath('C:\toolbox\fieldtrip-20210912\external\lagextraction\')      % if incompatibility with savefig

% Define main directories for raw data, processed results, clustering outputs, and figures
p.rawpath    = 'Z:\HCP\raw_preproc\';           % preprocessed data
p.datapath   = 'Z:\HCP\data\';                  % data results
p.kmeanspath = 'Z:\HCP\Juancho\Kmeans';     % kmeans results: Kmeans (per_condition)
p.outpath    = 'Z:\HCP\Juancho\Kmeans\Figures'; % Figure directory
p.figpath    = 'Z:\HCP\Juancho\Kmeans\Figures';
p.FNAT = 'Z:\HCP\Juancho\Kmeans\FNAT_max';

% Scan the data directory to extract subject IDs (skip '.' and '..')
cd (p.datapath)
temp = ls;
subs = deblank(temp(3:end,:)); % List of valid subject folders, all but  .  and  ..

% Define HCP tasks, their corresponding session numbers, and number of conditions per task
tasks    = {'Restin', 'StoryM', 'Wrkmem', 'Motort'};
sessions = {[3 4 5],[8 9],[6 7],[10 11]};
Ncondit  = [1, 2, 5, 5]; % Number of conditions: Rest, StoryM, Wrkmem, Motort

Nsub = size(subs,1); % Total number of subjects found
Ntask = length(tasks); % Total number of tasks

% Validate which subjects actually have data for each task
% Creates a nested cell array: substsk{task_idx}{valid_subject_idx} = subject_ID
substsk = {};
for t = 1:Ntask
    task = tasks{t};
    ct = 1;
    for s = 1: Nsub
        sub = deblank(subs(s,:));
        if exist([p.rawpath sub '\MEG\' tasks{t}])==7
            substsk{t}{ct} = sub;
            ct = ct+1;
        end
    end
end

% Define human-readable labels for each condition within each task
% Restin
condit_label{1}{1} = 'Restin';

% StoryM
condit_label{2}{1} = 'Story';
condit_label{2}{2} = 'Math';

% Wrkmem
condit_label{3}{1} = '0_back'; %'0-back'
condit_label{3}{2} = '2_back'; %'2-back'
condit_label{3}{3} = 'Face';
condit_label{3}{4} = 'Tools';
condit_label{3}{5} = 'WM_Fixation';

% Motort
condit_label{4}{1} = 'Left_hand';
condit_label{4}{2} = 'Left_foot';
condit_label{4}{3} = 'Right_hand';
condit_label{4}{4} = 'Right_foot';
condit_label{4}{5} = 'Fixation';

% Define numerical condition codes (used internally by custom HCP functions or loops)
condit_code={}; % Restin
condit_code{2} = [1 2];     % StoryM             % 1: Story, 2: Math
condit_code{3} = [1 2 1 2 0];     % Wrkmem     % 1: 0-back, 2: 2-back, 1: Faces, 2: Tools, 0: Fixation % N-back load and item content are mixed in the task.
condit_code{4} = [1 2 4 5 6];     % Motort     % 1-Left Hand,  2 - Left Foot, 4 - Right Hand. 5 - Right Foot, 6 - Fixation
%% Load template forward/inverse models
% Load precomputed source models from a reference subject (100307)
% These serve as spatial templates for all subjects' source-space analyses
cd('Z:\HCP\data\100307\Restin3')
load source_forward % only for template
load source_inverse % only for template
%% Freqanalysis
% Iterate through all validated subjects and tasks to perform source reconstruction and frequency analysis.
% For each subject-task pair, verify data availability, reconstruct source activity,
% and then compute frequency-domain source data while separating it by experimental conditions.
for s = 1:Nsub      
    sub = deblank(subs(s,:));
    for t = 1:Ntask
        task = tasks{t};
        cd([p.rawpath sub '\MEG'])
        if exist(tasks{t})==7
            for ses = sessions{t}(1):sessions{t}(end)
                HCP1_Source_reconstruction (sub, task, ses, p)
            end
             HCP2_Freqsource_100f (sub, task, p, 0)          % separar los datos según condiciones
        end
    end
end

%% K-means powerspectrum
% Perform the first stage of K-means clustering on power spectra for each task-condition combination.
% Configures analysis parameters and passes them to the clustering function to identify reproducible
% spectral profiles across subjects.
% Prepares data from all subjects
% and sessions of task to subsequetly perform Kmeans clustering in HCP4
for t = [1 2 3 4]
    for condit = 1:Ncondit(t)
        task = tasks{t};
        cfg.numtask = t;
        cfg.condit  = condit;
        cfg.task    = tasks{cfg.numtask};
        cfg.condit_label = condit_label;
        cfg.condit_code = condit_code{t};
        cfg.source.forward = source_forward;
        cfg.source.inverse = source;
        cfg.test    = i;
        HCP3_Kmeans_powsp_100f (substsk{t}, task, p, cfg)
    end
end
%% K-means powerspectrum clustering
% Finalize the K-means clustering process. Changes to the results directory, 
% assigns data points to clusters, and saves the organized clustering 
% outputs for each task and condition.
% Performs Kmeans clustering for each task and condition with data from all
% subjects pooled together.
for t = [1 2 3 4]
    for condit = 1:Ncondit(t)
        task = tasks{t};
        cfg.numtask = t;
        cfg.condit  = condit;
        cfg.task    = tasks{cfg.numtask};
        cfg.condit_label = condit_label;
        cfg.condit_code = condit_code{t};
        cfg.source.forward = source_forward;
        cfg.source.inverse = source;
        cfg.test    = i;
        cd (p.kmeanspath)
        HCP4_Kmeans_clustering_100f(substsk{t}, task, condit, p)
    end
end
%% HCP_Figure1_centroids
% Generate visualization figures for the K-means cluster centroids.
% Plots the average power spectra and corresponding spatial brain patterns 
% for each cluster.
% Figures are automatically closed after saving to free system memory.
for t = [1 2 3 4]
    for condit = 1:Ncondit(t)
        cfg = [];
        cfg.numtask = t;
        cfg.condit  = condit;
        cfg.task    = tasks{cfg.numtask};
        cfg.condit_label = condit_label;
        cfg.Nk      = 25;
        cfg.Nvox    = 1925;
        cfg.Nsub    = length(substsk{t});
        cfg.fig_powsp = 1;
        cfg.stats     = 0;
        cfg.fig_brainpatt     = 1;
        cfg.locgenerators_fig = 0;
        % cfg.locgenerators_vox = [2479 2472 2925 4816 4636 3578 3929 2066
        % 1717 4637 4634 5053 4634 5126 5117 4261 4633 5070 5047 5863 5477
        % 4459 4451 3722 4091];     % list of generators of brain rhythms,
        % Not used in HCP
        cfg.fig_gamma = 0;          % Supplementary Figure 1

        HCP_Figure1_centroids_100f (cfg, p)
        close all
    end
end

%% Compute natural frequencies for individual subjects
% Compute natural frequency maps for each subject individually. Each voxel
% is assigned a frequency based on the similarity of it's spectral profile
% to each of the clusters.
% Reloads spatial templates, iterates through all subjects/tasks/conditions, 
% and estimates subject-specific intrinsic brain oscillation Nat. freq.
cfg = [];
%cfg.Ntr     = 100;
cd('Z:\HCP\data\100307\Restin3')
load source_forward % only for template
load source_inverse % only for template

for s = 1: Nsub
    sub = deblank(subs(s,:));
    for t = [1 2 3 4]
        task = tasks{t};
        cd (p.kmeanspath)
        for condit = 1:Ncondit(t)
            %for i = 1:length(badk{t}{condit}) % HCP does not have badk
                cfg.Nk      = 25;
                cfg.Nvox    = 1925;
                cfg.Nsub    = length(substsk{t});
                cfg.numtask = t;
                cfg.condit  = condit;
                cfg.task    = tasks{cfg.numtask};
                cfg.condit_label = condit_label;
                cfg.condit_code = condit_code{t};
                cfg.source.forward = source_forward;
                cfg.source.inverse = source;
                cfg.test    = i;
                %cfg.badk    = badk{t}{condit}{i}; % HCP does not have badk
                HCP5_Naturalfreq_singlesubject_100f(cfg, substsk, task, condit, p, sub) % Uses Omega_neighbors_ly
            %end
        end
    end
end
%% Compute natural frequencies for all subjects (group map)
% Aggregate individual natural frequency estimates to generate group-level 
% average maps.
% Computes the spatial average of intrinsic oscillation frequencies across 
% all subjects for each task-condition combination.
cfg = [];
%cfg.Ntr     = 100;
cd('Z:\HCP\data\100307\Restin3')
load source_forward % only for template
load source_inverse % only for template
% Define the root directory where the subject folders are located
rootDir = 'Z:\HCP\data'; % datapath
for t = [1 2 3 4]
    task = tasks{t};
    cd (p.kmeanspath)
    for condit = 1:Ncondit(t)
        cfg.Nk      = 25;
        cfg.Nvox    = 1925;
        cfg.Nsub    = length(substsk{t});
        cfg.numtask = t;
        cfg.condit  = condit;
        cfg.task    = tasks{cfg.numtask};
        cfg.condit_label = condit_label;
        cfg.condit_code = condit_code{t};
        cfg.source.forward = source_forward;
        cfg.source.inverse = source;
        HCP6_Naturalfreq_groupmap_100f (cfg, substsk, task, condit, p, rootDir) % Unchanged for now
    end
end

%% Structure fnat dataframe
% Extract and compile natural frequency values and associated power spectra 
% for all subjects.
% Uses a try-catch block to gracefully handle subjects with missing or 
% invalid data, concatenates successful results into a unified 
% multi-dimensional array, and saves it to disk.
cfg = [];
%cfg.Ntr     = 100; % If you want to use just 100 tr.
cd('Z:\HCP\data\100307\Restin3')
load source_forward % only for template
load source_inverse % only for template

for t = 1:4
    fnat_full{t} = [];
    task = tasks{t};
    cd (p.kmeanspath)        
    for condit = 1:Ncondit(t)
        cts = 1;
        for s = 1: Nsub
            sub = deblank(subs(s,:));
            fnat_sub = [];
            cfg.Nk      = 25;
            cfg.Nvox    = 1925;
            cfg.Nsub    = length(substsk{t});
            cfg.numtask = t;
            cfg.condit  = condit;
            cfg.task    = tasks{cfg.numtask};
            cfg.condit_label = condit_label;
            cfg.condit_code = condit_code{t};
            cfg.test    = i;
            try
                [fnat_sub, fnat_powsp] = HCP7_Naturalfreq_singlesubject_fnat(cfg, substsk, task, condit, p, sub);
            end
            if ~isempty(fnat_sub)
                fnat_full{t}(condit,cts,:) = fnat_sub; 
                cts = cts+1;
            end
        end
    end
end

cd([p.kmeanspath '\FNAT'])
save fnat_full fnat_full
%% Condition subdivision (fnat_struct)
% Reorganize the compiled dataset into a structured MATLAB object for 
% downstream analysis.
% Dynamically generates field names based on task-condition labels, removes 
% singleton dimensions using squeeze, and saves the final structured 
% dataset (ready for Part 2: statistical analysis).

% Load natural frequencies in single dataframe and give it structure
cd([p.kmeanspath '\FNAT'])
load fnat_full
% empty struct object
fnat_struct = struct();
% for loop to iterate over tasks and conditions
for task_number = 1:length(tasks)
    task_name = tasks{task_number};
    num_conditions = Ncondit(task_number);
    
    for condition = 1:num_conditions
        condition_name = condit_label{task_number}{condition};
        % Construir el nombre del campo dinámicamente
        field_name = sprintf('%s_%s', task_name, condition_name);
        % Aplicar squeeze y almacenar el resultado
        fnat_struct.(field_name) = squeeze(fnat_full{task_number}(condition,:,:));

    end
end

% Mostrar resultados
disp(fnat_struct);
save fnat_struct fnat_struct

%% HCP (Part 2)
% ------------------------------------------------------------------------%
% Part 2
% ------------------------------------------------------------------------%
% This script (2) implements a complete statistical and visualization work-
% flow for Human Connectome Project (HCP) MEG data. It begins by configuring 
% data paths, loading preprocessed natural frequency datasets, and defining 
% frequency bins of interest. The pipeline then executes voxel-wise 
% condition comparisons (computing T-statistics and Cohen's d and other 
% effect sizes), followed by non-parametric permutation testing to establish 
% statistical thresholds for multiple comparison correction. Finally, it 
% generates and exports spatial brain maps (via MRIcroN) and condition-
% specific spectral plots for significant clusters and individual voxels, 
% preparing the results for downstream analysis and publication.
% If Part 1 has been run already, you can run part 2 directly.
% -----------------------------------------------------------------------
%% Paths and files
% Configure MATLAB paths for FieldTrip, and custom HCP scripts.
% Define directory structure for raw/processed data, scan the data folder 
% to list subject IDs, and map HCP tasks to their corresponding sessions 
% and experimental conditions. Validates data availability per subject 
% and creates structured labels/codes for downstream processing.
% Add required toolboxes and custom script directories to MATLAB's search path
addpath('Z:\Toolbox\fieldtrip-20230118\') %fieldtrip              
addpath(genpath('Z:\HCP\Juancho\')) % Work folder
addpath(genpath('Z:\HCP\Juancho\github\subfxs')) % for function: Omega_neighbors_ly 
% % rmpath('C:\toolbox\fieldtrip-20210912\external\lagextraction\')      % if incompatibility with savefig

% Define main directories for raw data, processed results, clustering outputs, and figures
p.rawpath    = 'Z:\HCP\raw_preproc\';           % preprocessed data
p.datapath   = 'Z:\HCP\data\';                  % data results
p.kmeanspath = 'Z:\HCP\Juancho\Kmeans';     % kmeans results: Kmeans (per_condition)
p.outpath    = 'Z:\HCP\Juancho\Kmeans\Figures'; % Figure directory
p.figpath    = 'Z:\HCP\Juancho\Kmeans\Figures';

cd (p.datapath)
temp = ls;
subs = deblank(temp(3:end,:));                      % all but  .  and  ..  
tasks    = {'Restin', 'StoryM', 'Wrkmem', 'Motort'};
sessions = {[3 4 5],[8 9],[6 7],[10 11]};
Ncondit  = [1, 2, 5, 5]; %[1, 2, 4, 5]

Nsub = size(subs,1);
Ntask = length(tasks);

substsk = {};
for t = 1:Ntask
    task = tasks{t};
    ct = 1;
    for s = 1: Nsub
        sub = deblank(subs(s,:));
        if exist([p.rawpath sub '\MEG\' tasks{t}])==7
            substsk{t}{ct} = sub;
            ct = ct+1;
        end
    end
end

% Resting-state
condit_label{1}{1} = 'Restin';

% StoryM
condit_label{2}{1} = 'Story';
condit_label{2}{2} = 'Math';

% Wrkmem
condit_label{3}{1} = '0_back'; %'0-back'
condit_label{3}{2} = '2_back'; %'2-back'
condit_label{3}{3} = 'Face';
condit_label{3}{4} = 'Tools';
condit_label{3}{5} = 'WM_Fixation';

% Motort
condit_label{4}{1} = 'Left_hand';
condit_label{4}{2} = 'Left_foot';
condit_label{4}{3} = 'Right_hand';
condit_label{4}{4} = 'Right_foot';
condit_label{4}{5} = 'Fixation';

% Codes for different conditions
condit_code={};
condit_code{2} = [1 2];     % StoryM             % 1: Story, 2: Math % TO DO
condit_code{3} = [1 2 1 2 0];     % Wrkmem             % TO DO
condit_code{4} = [1 2 4 5 6];     % Motort     % 1-Left Hand,  2 - Left Foot, 4 - Right Hand. 5 - Right Foot, 6 - Fixation
%% Templates
% Load precomputed forward and inverse source models from a reference subject (100307).
% These serve as spatial templates for mapping results back to source space across all subjects.
cd('Z:\HCP\data\100307\Restin3')
load source_forward % only for template
load source_inverse % only for template
%% Frequencies for difference analysis
% Load the precomputed natural frequency dataset (fnat_struct) from the K-means pipeline.
% Define the exponential frequency grid (foi) for spectral analysis and compute a symmetric adjacency matrix
% that defines spatial neighbors for each voxel, enabling spatial smoothing and cluster-based statistics.
cd('Z:\HCP\Juancho\Kmeans\FNAT_max')
load fnat_struct

f1         = 0.55:0.05:0.9;          % frequencies are limited by trial length (1.2 s -> min freq 2.6 Hz with 3 cycles; lower frequencies analyzed with 2 up to 3 cycles)
foi1      = exp(f1);

f2         = 0.95:0.05:4.6; %3.55;
foi2      = exp(f2);

foi = [foi1 foi2];

% For visualization purposes:
f1 = 0.5;
f2 = 4.6; %3.50;
d  = 0.2;
foi1 = exp(f1:d:f2);
foi1x = exp(f1+d/2:d:f2-d/2);
foi2 = interp(foi1(1:end-1),10);
foi2 = foi2(1:end-3)

% Symetric matrix for all the neighbours of each voxel.
neighbor_matrix = find_neighbours(1.5); % Available in subfx folder
neighbor_matrix = neighbor_matrix.connmat;

%% CONDITION COMPARISON AREA ----------------------------------------------

%% Condition comparison loop (FNAT)
% Compute voxel-wise statistical differences between predefined task conditions.
% Iterates through comparison sets (Motor, Working Memory, Story, Rest), calculates Cohen's d,
% T-sum statistics, and effect sizes (Cliff's delta: delta_M) for each of the 1925 source voxels, and stores 
% results in a structured array for downstream thresholding and visualization.
warning('off'); % Disable warnings
tic; % for timming
% Create an empty struct to store the results
cohen_d_struct = struct();

% Get all field names
fields = fieldnames(fnat_struct);

% Compare Motort_Fixation against all conditions containing Motort_
motort_fields = fields(contains(fields, 'Motort_'));
reference_condition_Motort = 'Motort_Fixation';

for i = 1:length(motort_fields)
    current_condition = motort_fields{i};
    
    if strcmp(current_condition, reference_condition_Motort)
        continue;
    end
    
    cohen_d_per_voxel = NaN(1, 1925);
    T_sum_per_voxel = NaN(1, 1925);
    normality_p_per_voxel = NaN(length(foi2), 1925);
    
    for vox = 1:1925
        [normality_p, T_sum, cohen_d_M, delta_M] = HCP9_Naturalfreq_cohen_d_max(fnat_struct, reference_condition_Motort, current_condition, vox, foi1, foi1x, foi2, neighbor_matrix, substsk, 'rep_measures');
        cohen_d_per_voxel(vox) = cohen_d_M;
        T_sum_per_voxel(vox) = T_sum;
        delta_M_per_voxel(vox) = delta_M;
        normality_p_per_voxel(:,vox) = normality_p;
    end
    
    comparison_name = [reference_condition_Motort '_vs_' current_condition];
    cohen_d_struct.(comparison_name).T_sum = T_sum_per_voxel;
    cohen_d_struct.(comparison_name).cohen_D = cohen_d_per_voxel;
    cohen_d_struct.(comparison_name).delta_M = delta_M_per_voxel;
    cohen_d_struct.(comparison_name).normality_p = normality_p_per_voxel;
end

% Compare Wrkmem_WM_Fixation against all conditions containing Wrkmem_
wrkmem_fields = fields(contains(fields, 'Wrkmem_'));
reference_condition_Wrkmem = 'Wrkmem_WM_Fixation';

for i = 1:length(wrkmem_fields)
    current_condition = wrkmem_fields{i};
    
    if strcmp(current_condition, reference_condition_Wrkmem)
        continue;
    end
    
    cohen_d_per_voxel = NaN(1, 1925);
    T_sum_per_voxel = NaN(1, 1925);
    normality_p_per_voxel = NaN(length(foi2), 1925);
    
    for vox = 1:1925
        [normality_p, T_sum, cohen_d_M, delta_M] = HCP9_Naturalfreq_cohen_d_max(fnat_struct, reference_condition_Wrkmem, current_condition, vox, foi1, foi1x, foi2, neighbor_matrix, substsk, 'rep_measures');
        cohen_d_per_voxel(vox) = cohen_d_M;
        T_sum_per_voxel(vox) = T_sum;
        delta_M_per_voxel(vox) = delta_M;
        normality_p_per_voxel(:,vox) = normality_p;
    end
    
    comparison_name = [reference_condition_Wrkmem '_vs_' current_condition];
    cohen_d_struct.(comparison_name).T_sum = T_sum_per_voxel;
    cohen_d_struct.(comparison_name).cohen_D = cohen_d_per_voxel;
    cohen_d_struct.(comparison_name).delta_M = delta_M_per_voxel;
    cohen_d_struct.(comparison_name).normality_p = normality_p_per_voxel;
end

% Compare Wrkmem_0_back VS Wrkmem_2_back and Wrkmem_Face VS Wrkmem_Tools
comparisons = {
    'Wrkmem_0_back', 'Wrkmem_2_back';
    'Wrkmem_Face', 'Wrkmem_Tools'
};

for i = 1:size(comparisons, 1)
    condition1 = comparisons{i, 1};
    condition2 = comparisons{i, 2};

    cohen_d_per_voxel = NaN(1, 1925);
    T_sum_per_voxel = NaN(1, 1925);
    normality_p_per_voxel = NaN(length(foi2), 1925);
    
    for vox = 1:1925
        [normality_p, T_sum, cohen_d_M, delta_M] = HCP9_Naturalfreq_cohen_d_max(fnat_struct, condition1, condition2, vox, foi1, foi1x, foi2, neighbor_matrix, substsk, 'rep_measures');
        cohen_d_per_voxel(vox) = cohen_d_M;
        T_sum_per_voxel(vox) = T_sum;
        delta_M_per_voxel(vox) = delta_M;
        normality_p_per_voxel(:,vox) = normality_p;
    end
    
    comparison_name = [condition1 '_vs_' condition2];
    cohen_d_struct.(comparison_name).T_sum = T_sum_per_voxel;
    cohen_d_struct.(comparison_name).cohen_D = cohen_d_per_voxel;
    cohen_d_struct.(comparison_name).delta_M = delta_M_per_voxel;
    cohen_d_struct.(comparison_name).normality_p = normality_p_per_voxel;
end

% Compare StoryM_Story against StoryM_Math
reference_condition_StoryM = 'StoryM_Story';
current_condition_StoryM = 'StoryM_Math';

cohen_d_per_voxel = NaN(1, 1925);
T_sum_per_voxel = NaN(1, 1925);
normality_p_per_voxel = NaN(length(foi2), 1925);

for vox = 1:1925
    [normality_p, T_sum, cohen_d_M, delta_M] = HCP9_Naturalfreq_cohen_d_max(fnat_struct, reference_condition_StoryM, current_condition_StoryM, vox, foi1, foi1x, foi2, neighbor_matrix, substsk, 'rep_measures');    
    cohen_d_per_voxel(vox) = cohen_d_M;
    T_sum_per_voxel(vox) = T_sum;
    delta_M_per_voxel(vox) = delta_M;
    normality_p_per_voxel(:,vox) = normality_p;
end

comparison_name = [reference_condition_StoryM '_vs_' current_condition_StoryM];
cohen_d_struct.(comparison_name).T_sum = T_sum_per_voxel;
cohen_d_struct.(comparison_name).cohen_D = cohen_d_per_voxel;
cohen_d_struct.(comparison_name).delta_M = delta_M_per_voxel;
cohen_d_struct.(comparison_name).normality_p = normality_p_per_voxel;

% Compare StoryM_Story VS Restin_Restin and StoryM_Math VS Restin_Restin
comparisons = {
    'Restin_Restin', 'StoryM_Story';
    'Restin_Restin', 'StoryM_Math'
};

for i = 1:size(comparisons, 1)
    condition1 = comparisons{i, 1};
    condition2 = comparisons{i, 2};

    cohen_d_per_voxel = NaN(1, 1925);
    T_sum_per_voxel = NaN(1, 1925);
    normality_p_per_voxel = NaN(length(foi2), 1925);
    
    for vox = 1:1925
        [normality_p, T_sum, cohen_d_M, delta_M] = HCP9_Naturalfreq_cohen_d_max(fnat_struct, condition1, condition2, vox, foi1, foi1x, foi2, neighbor_matrix, substsk, 'rep_measures');
        cohen_d_per_voxel(vox) = cohen_d_M;
        T_sum_per_voxel(vox) = T_sum;
        delta_M_per_voxel(vox) = delta_M;
        normality_p_per_voxel(:,vox) = normality_p;
    end
    
    comparison_name = [condition1 '_vs_' condition2];
    cohen_d_struct.(comparison_name).T_sum = T_sum_per_voxel;
    cohen_d_struct.(comparison_name).cohen_D = cohen_d_per_voxel;
    cohen_d_struct.(comparison_name).delta_M = delta_M_per_voxel;
    cohen_d_struct.(comparison_name).normality_p = normality_p_per_voxel;
end
toc;
% Save the result into a single file
if ~exist('Z:\HCP\Juancho\Kmeans\FNAT_max', 'dir') % FNAT_max: max effect size method for Natural Frequency Analysis
    mkdir('Z:\HCP\Juancho\Kmeans\FNAT_max');
end
cd('Z:\HCP\Juancho\Kmeans\FNAT_max')
save('cohen_d_struct.mat', 'cohen_d_struct');
%% Save image to MRIcroN
% Load the computed statistics, organize them by task type, and generate NIfTI-formatted brain maps
% for visualization in MRIcroN. Processes Cohen's d and Cliff's delta metrics for each comparison and
% exports spatial maps into task-specific directories for neuroimaging review.
load('Z:\HCP\Juancho\Kmeans\FNAT_max\cohen_d_struct.mat')

% Get all field names in cohen_d_struct
field_names = fieldnames(cohen_d_struct);

% Define the conditions and their respective folders
conditions = {'Motort_', 'Wrkmem_', 'StoryM_'};
base_dir = 'Z:\HCP\Juancho\Kmeans\FNAT_max\img\MRICroN_nii\FNAT'; % New folder for images

% Create folders if they do not exist
for j = 1:length(conditions)
    condition = conditions{j};
    condition_dir = fullfile(base_dir, condition);
    if ~exist(condition_dir, 'dir')
        mkdir(condition_dir);
    end
end

% Iterate over each field and apply the save_cohen_d_image function
for i = 1:length(field_names)
    % Get current name
    current_name = field_names{i};
    
    % Determine the corresponding folder
    if contains(current_name, 'Motort_')
        condition_dir = fullfile(base_dir, 'Motort_');
    elseif contains(current_name, 'Wrkmem_')
        condition_dir = fullfile(base_dir, 'Wrkmem_');
    elseif contains(current_name, 'StoryM_')
        condition_dir = fullfile(base_dir, 'StoryM_');
    else
        continue; % If it does not match any condition, skip to the next
    end
    
    % Change to the corresponding directory
    cd(condition_dir);
    
    % Call the save_cohen_d_image function with the current name
    %HCP9_Naturalfreq_Plot_MRIcroN(cohen_d_struct, current_name, source, source_forward, 'T_sum');
    HCP9_Naturalfreq_Plot_MRIcroN(cohen_d_struct, current_name, source, source_forward, 'cohen_D');
    HCP9_Naturalfreq_Plot_MRIcroN(cohen_d_struct, current_name, source, source_forward, 'delta_M');
end
%% PERMUTATION ANALYSIS ---------------------------------------------------

%% Condition comparison loop with permutations (FNAT)
% Run non-parametric permutation testing (1000 iterations) to generate null distributions.
% Computes maximum statistic thresholds across all voxels for multiple comparison correction,
% storing the max T-sum, Cohen's d, and delta_M values for each comparison to establish 
% statistically robust significance cutoffs.
warning('off'); % Disable warnings
n_perm = 1000;
cohen_d_struct_1000_perm_max = struct();

tic;
% Get all field names
fields = fieldnames(fnat_struct);

% Compare Motort_Fixation against all conditions containing Motort_
motort_fields = fields(contains(fields, 'Motort_'));
reference_condition_Motort = 'Motort_Fixation';

for i = 1:length(motort_fields)
    condition2 = motort_fields{i};
    
    if strcmp(condition2, reference_condition_Motort)
        continue;
    end

    [max_t_values, max_cohen_d_values, max_delta_M_values] = HCP9_Naturalfreq_permuted_all_vox_max(fnat_struct, reference_condition_Motort, condition2, foi1, foi1x, foi2, neighbor_matrix, substsk, n_perm, 'rep_measures');
    
    comparison_name = sprintf('%s_vs_%s_%dperm_max_values', reference_condition_Motort, condition2, n_perm);
    cohen_d_struct_1000_perm_max.(comparison_name).T_sum_max = max_t_values;
    cohen_d_struct_1000_perm_max.(comparison_name).D_mean_max = max_cohen_d_values;
    cohen_d_struct_1000_perm_max.(comparison_name).delta_mean_max = max_delta_M_values;
end

% Compare Wrkmem_WM_Fixation against all conditions containing Wrkmem_
wrkmem_fields = fields(contains(fields, 'Wrkmem_'));
reference_condition_Wrkmem = 'Wrkmem_WM_Fixation';

for i = 1:length(wrkmem_fields)
    condition2 = wrkmem_fields{i};
    
    if strcmp(condition2, reference_condition_Wrkmem)
        continue;
    end

    [max_t_values, max_cohen_d_values, max_delta_M_values] = HCP9_Naturalfreq_permuted_all_vox_max(fnat_struct, reference_condition_Wrkmem, condition2, foi1, foi1x, foi2, neighbor_matrix, substsk, n_perm, 'rep_measures');
    
    comparison_name = sprintf('%s_vs_%s_%dperm_max_values', reference_condition_Wrkmem, condition2, n_perm);
    cohen_d_struct_1000_perm_max.(comparison_name).T_sum_max = max_t_values;
    cohen_d_struct_1000_perm_max.(comparison_name).D_mean_max = max_cohen_d_values;
    cohen_d_struct_1000_perm_max.(comparison_name).delta_mean_max = max_delta_M_values;
end

% Compare Wrkmem_0_back VS Wrkmem_2_back and Wrkmem_Face VS Wrkmem_Tools
comparisons = {
    'Wrkmem_0_back', 'Wrkmem_2_back';
    'Wrkmem_Face', 'Wrkmem_Tools'
};

for i = 1:size(comparisons, 1)
    condition1 = comparisons{i, 1};
    condition2 = comparisons{i, 2};
    [max_t_values, max_cohen_d_values, max_delta_M_values] = HCP9_Naturalfreq_permuted_all_vox_max(fnat_struct, condition1, condition2, foi1, foi1x, foi2, neighbor_matrix, substsk, n_perm, 'rep_measures');
    
    comparison_name = sprintf('%s_vs_%s_%dperm_max_values', condition1, condition2, n_perm);
    cohen_d_struct_1000_perm_max.(comparison_name).T_sum_max = max_t_values;
    cohen_d_struct_1000_perm_max.(comparison_name).D_mean_max = max_cohen_d_values;
    cohen_d_struct_1000_perm_max.(comparison_name).delta_mean_max = max_delta_M_values;
end

% Compare StoryM_Story against StoryM_Math
reference_condition_StoryM = 'StoryM_Story';
current_condition_StoryM = 'StoryM_Math';

[max_t_values, max_cohen_d_values, max_delta_M_values] = HCP9_Naturalfreq_permuted_all_vox_max(fnat_struct, reference_condition_StoryM, current_condition_StoryM, foi1, foi1x, foi2, neighbor_matrix, substsk, n_perm, 'rep_measures');

comparison_name = sprintf('%s_vs_%s_%dperm_max_values', reference_condition_StoryM, current_condition_StoryM, n_perm);
cohen_d_struct_1000_perm_max.(comparison_name).T_sum_max = max_t_values;
cohen_d_struct_1000_perm_max.(comparison_name).D_mean_max = max_cohen_d_values;
cohen_d_struct_1000_perm_max.(comparison_name).delta_mean_max = max_delta_M_values;


% Compare StoryM_Story VS Restin_Restin and StoryM_Math VS Restin_Restin
comparisons = {
    'Restin_Restin', 'StoryM_Story';
    'Restin_Restin', 'StoryM_Math'
};

for i = 1:size(comparisons, 1)
    condition1 = comparisons{i, 1};
    condition2 = comparisons{i, 2};
    [max_t_values, max_cohen_d_values, max_delta_M_values] = HCP9_Naturalfreq_permuted_all_vox_max(fnat_struct, condition1, condition2, foi1, foi1x, foi2, neighbor_matrix, substsk, n_perm, 'rep_measures');
    
    comparison_name = sprintf('%s_vs_%s_%dperm_max_values', condition1, condition2, n_perm);
    cohen_d_struct_1000_perm_max.(comparison_name).T_sum_max = max_t_values;
    cohen_d_struct_1000_perm_max.(comparison_name).D_mean_max = max_cohen_d_values;
    cohen_d_struct_1000_perm_max.(comparison_name).delta_mean_max = max_delta_M_values;
end

toc;

% Save the result into a single file
cd('Z:\HCP\Juancho\Kmeans\FNAT_max')
save('cohen_d_struct_1000_perm_max.mat', 'cohen_d_struct_1000_perm_max');
%% GRAPH SECTION ----------------------------------------------------------

%% Plotting loop for saved comparisons
% Generate whole-brain statistical maps by applying permutation-derived significance thresholds.
% Uses the 95th percentile of the max-statistic null distribution as the cutoff for Cohen's d and delta_M,
% then exports thresholded spatial maps for each condition comparison.
% Load the differences and permutations
load('Z:\HCP\Juancho\Kmeans\FNAT_max\cohen_d_struct.mat')
load('Z:\HCP\Juancho\Kmeans\FNAT_max\cohen_d_struct_1000_perm_max.mat')

fields = fieldnames(cohen_d_struct)
comparisons = fields%(contains(fields, 'StoryM_'))

fields2 = fieldnames(cohen_d_struct_1000_perm_max)
comparisons2 = fields2%(contains(fields2, 'StoryM_'))

% Generate and save images for cohen_d_struct
for i = 1:length(comparisons)
    comparison_name = comparisons{i};
    comparison_name2 = comparisons2{i};

    if ~exist('Z:\HCP\Juancho\Kmeans\FNAT_max\img\MRICroN_nii\sig', 'dir')
        mkdir('Z:\HCP\Juancho\Kmeans\FNAT_max\img\MRICroN_nii\sig');
    end

    cd('Z:\HCP\Juancho\Kmeans\FNAT_max\img\MRICroN_nii\sig')
    cutoff_point_D = quantile(cohen_d_struct_1000_perm_max.(comparison_name2).D_mean_max, 0.95);
    HCP9_Naturalfreq_Plot_MRIcroN_2(cohen_d_struct, comparison_name, source, source_forward, 'cohen_D', cutoff_point_D);
    cutoff_point_delta = quantile(cohen_d_struct_1000_perm_max.(comparison_name2).delta_mean_max, 0.95);
    HCP9_Naturalfreq_Plot_MRIcroN_2(cohen_d_struct, comparison_name, source, source_forward, 'delta_M', cutoff_point_delta);
end
%% Plotting loop for saved comparisons local max
% Identify and visualize local maxima within significant clusters for each condition comparison.
% Applies the same permutation-based thresholds but restricts output to peak voxels only,
% saving focused spatial maps that highlight the strongest condition-specific effects.
fields = fieldnames(cohen_d_struct)
comparisons = fields%(contains(fields, 'StoryM_'))

fields2 = fieldnames(cohen_d_struct_1000_perm_max)
comparisons2 = fields2%(contains(fields2, 'StoryM_'))

% Generate and save images for cohen_d_struct
for i = 1:length(comparisons)
    comparison_name = comparisons{i};
    comparison_name2 = comparisons2{i};

    if ~exist('Z:\HCP\Juancho\Kmeans\FNAT_max\img\MRICroN_nii\sig\localmax', 'dir')
        mkdir('Z:\HCP\Juancho\Kmeans\FNAT_max\img\MRICroN_nii\sig\localmax');
    end

    cd('Z:\HCP\Juancho\Kmeans\FNAT_max\img\MRICroN_nii\sig\localmax')
    cutoff_point_D = quantile(cohen_d_struct_1000_perm_max.(comparison_name2).D_mean_max, 0.95);
    HCP9_Naturalfreq_Plot_MRIcroN_LocalMax(cohen_d_struct, comparison_name, source, source_forward, 'cohen_D', cutoff_point_D);
    cutoff_point_delta = quantile(cohen_d_struct_1000_perm_max.(comparison_name2).delta_mean_max, 0.95);
    HCP9_Naturalfreq_Plot_MRIcroN_LocalMax(cohen_d_struct, comparison_name, source, source_forward, 'delta_M', cutoff_point_delta);
end
%% PLOT specific voxels
%% Plotting loop for saved comparisons local max 
% Extract voxel indices corresponding to local maxima, sort them by effect size (Cohen's d),
% and prepare a structured list of significant voxels for targeted spectral difference plotting.
% Load the differences and permutations
load('Z:\HCP\Juancho\Kmeans\FNAT_max\cohen_d_struct.mat')
load('Z:\HCP\Juancho\Kmeans\FNAT_max\cohen_d_struct_1000_perm_max.mat')

fields = fieldnames(cohen_d_struct);
comparisons = fields;  

fields2 = fieldnames(cohen_d_struct_1000_perm_max);
comparisons2 = fields2;  

% Initialize structure to store only the numbers of the selected voxels
selected_voxels_struct = struct();

% Create folder if it does not exist
output_dir = 'Z:\HCP\Juancho\Kmeans\FNAT_max\img\MRICroN_nii\sig\localmax';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

cd(output_dir)

% Generate and save images for each comparison
for i = 1:length(comparisons)
    comparison_name = comparisons{i};
    comparison_name2 = comparisons2{i};

    % Initialize substructure for this comparison
    selected_voxels_struct0.(comparison_name) = [];

    % Define thresholds based on quantiles
    cutoff_point_D_max = quantile(cohen_d_struct_1000_perm_max.(comparison_name2).D_mean_max, 0.95);
    cutoff_point_delta_max = quantile(cohen_d_struct_1000_perm_max.(comparison_name2).delta_mean_max, 0.95);

    % Local maxima
    [~, ~, locmx_D_max, voxel_numbers_D_max] = HCP9_Naturalfreq_Plot_MRIcroN_LocalMax3(...
        cohen_d_struct, comparison_name, source, source_forward, 'cohen_D', cutoff_point_D_max, 6);
    
    [~, ~, locmx_delta_max, voxel_numbers_delta_max] = HCP9_Naturalfreq_Plot_MRIcroN_LocalMax3(...
        cohen_d_struct, comparison_name, source, source_forward, 'delta_M', cutoff_point_delta_max, 6);

    % Save only the voxel number within the range 1-1925
    selected_voxels_struct0.(comparison_name).cohen_D_max = voxel_numbers_D_max;  
    selected_voxels_struct0.(comparison_name).delta_max = unique(voxel_numbers_delta_max(:));  
end

% Display the selected voxels structure
disp('Selected Voxels Structure:');
disp(selected_voxels_struct0);

% Initialize a structure to store the selected voxels
selected_voxels_struct1 = struct();
field_names = fieldnames(cohen_d_struct);
field_names2 = fieldnames(cohen_d_struct_1000_perm_max);

% Iterate over each field of cohen_d_struct
for i = 1:length(field_names)
    % Get the current field name
    current_name = field_names{i};
    current_name2 = field_names2{i};
    
    % Extract data from the current field (Cohen's d for each voxel)
    data = cohen_d_struct.(current_name).cohen_D;
    
    % Select the localmax voxels
    selected_voxels = selected_voxels_struct0.(current_name).cohen_D_max
    
    % Extract the Cohen's d values corresponding to the selected voxels
    selected_d_values = data(selected_voxels);
    
    % Sort the selected voxels in descending order of Cohen's d
    [~, sorted_indices] = sort(selected_d_values, 'descend');
    selected_voxels_sorted = selected_voxels(sorted_indices);
    
    % Save the sorted results in the new structure
    selected_voxels_struct1.(current_name) = selected_voxels_sorted;
end

% Show the results for verification
disp(selected_voxels_struct1);

%% Manually selecting voxels for visualization
% Override automated voxel selection with manually curated indices for precise, publication-ready figures.
% Maps key brain regions to specific experimental contrasts based on anatomical and functional relevance.
% The following voxel indices where manually selected through MRICroN image
% inspection. Selected voxel where all significant and either local maxima 
% values in effect size or symetrical to those and/or representative of a
% given area.
selected_voxels_struct0.Motort_Fixation_vs_Motort_Left_hand = [1731 1593 1482 1644 635];
selected_voxels_struct0.Motort_Fixation_vs_Motort_Left_foot = [156 1068 93 717];
selected_voxels_struct0.Motort_Fixation_vs_Motort_Right_hand = [1853 1762 1482 1644 65 60 383 201 7 20];
selected_voxels_struct0.Motort_Fixation_vs_Motort_Right_foot = [1042 62 1357 45];

selected_voxels_struct0.Wrkmem_0_back_vs_Wrkmem_2_back = [1524 1332 1840 11];
selected_voxels_struct0.Wrkmem_Face_vs_Wrkmem_Tools = [724 1368 669 975 39];
selected_voxels_struct0.Wrkmem_WM_Fixation_vs_Wrkmem_0_back = [1743 384 1060 331 22];
selected_voxels_struct0.Wrkmem_WM_Fixation_vs_Wrkmem_2_back = [1743 1756 1795 379 162 65];
selected_voxels_struct0.Wrkmem_WM_Fixation_vs_Wrkmem_Face = [1743 1761 850 331 383];
selected_voxels_struct0.Wrkmem_WM_Fixation_vs_Wrkmem_Tools = [1743 1719 715 1875 401 1060];

selected_voxels_struct0.StoryM_Story_vs_StoryM_Math = [1563 761 1412 1134 63];
selected_voxels_struct0.Restin_Restin_vs_StoryM_Story = [649 867 1581 1586 351 153 68];
selected_voxels_struct0.Restin_Restin_vs_StoryM_Math = [1494 1306 649 867];
disp(selected_voxels_struct0);

%save selected_voxels_struct0.mat 'selected_voxels_struct0'

%% TEST IMAGES
% Comprehensive testing loop: generates condition-specific plots, applies permutation thresholds,
% exports local maxima visualizations, and organizes outputs into structured task/condition subfolders.
% Ensures all plotting pipelines run correctly across the full dataset before manuscript preparation.
selected_voxels_struct1 = selected_voxels_struct0; % Optional
output_folder = 'Z:\HCP\Juancho\Kmeans\FNAT_max\img\test_img';
%conditions = {'Motort_', 'Wrkmem_', 'StoryM_'};
tasks_ = {'Restin', 'StoryM', 'Wrkmem', 'Motort'};
comparisons = fieldnames(selected_voxels_struct1);
comparisons2 = fieldnames(cohen_d_struct_1000_perm_max);

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

for i = 1:length(comparisons)
    comparison_name = comparisons{i};
    
    % Initialize task found flag
    found_task = false;
    
    % Search which task is at the beginning of the comparison name
    for j = 1:length(tasks_)
        task_candidate = tasks_{j};
        if startsWith(comparison_name, [task_candidate '_'])
            found_task = true;
            break;
        end
    end
    
    if ~found_task
        continue; % Skip if no valid task is found
    end
    
    % Extract conditions: part before and after '_vs_'
    parts = strsplit(comparison_name, '_vs_');
    cond1 = parts{1};  % E.g., 'Motort_Fixation'
    cond2 = parts{2};  % E.g., 'Motort_Left_hand'

    % Remove the task prefix from each condition
    cond1_clean = strrep(cond1, [task_candidate '_'], '');
    cond2_clean = strrep(cond2, [task_candidate '_'], '');

    % Folder name: Fixation_vs_Left_hand
    condition_subfolder = sprintf('%s_vs_%s', cond1_clean, cond2_clean);
    
    % Full path: output_folder / task / condition_subfolder
    condition_dir = fullfile(output_folder, task_candidate, condition_subfolder);
    
    % Create folder if it does not exist
    if ~exist(condition_dir, 'dir')
        mkdir(condition_dir);
    end
    
    % Get the comparison name
    comparison_name = comparisons{i};
    comparison_name2 = comparisons2{i};
    
    % Split the name to extract the two conditions
    parts = split(comparison_name, '_vs_');
    
    if length(parts) ~= 2
        warning('Unexpected format in comparison name: %s', comparison_name);
        continue;
    end
    
    condition1 = parts{1};
    condition2 = parts{2};

    for j = 1:length(selected_voxels_struct1.(comparison_name))
        % Call the function with the condition names
        HCP9_Naturalfreq_Plot_Conditions_max_final(fnat_struct, ...
            condition1, condition2, selected_voxels_struct1.(comparison_name)(j), ...
            foi1, foi1x, foi2, neighbor_matrix, substsk, true, 'rep_measures');

        % Save the two generated plots
        fig1 = gcf; % Get the first figure (assuming it is the active one)
        % Save the figure with a name based on `comparison_name` and the first element of the structure
        %saveas(fig1, fullfile(condition_dir, sprintf('%s_%d_D.png', comparison_name, selected_voxels_struct1.(comparison_name)(j))));
        close(fig1); % Close the figure to free memory

        fig2 = gcf; % The second figure should now be active
        %saveas(fig2, fullfile(condition_dir, sprintf('%s_%d_VS.png', comparison_name, selected_voxels_struct1.(comparison_name)(j))));
        print(fig2, fullfile(condition_dir, sprintf('%s_%d_VS.tiff', comparison_name, selected_voxels_struct1.(comparison_name).cohen_D_max(j))), '-dtiff', '-r300');
        close(fig2); % Close the second figure

        cutoff_point_D = quantile(cohen_d_struct_1000_perm_max.(comparison_name2).D_mean_max, 0.95);
        plot_and_save_images_max_vistas_separadas(comparison_name, cohen_d_struct, 'cohen_D', source, source_forward, cutoff_point_D, condition_dir)

        voxel_index = selected_voxels_struct1.(comparison_name).cohen_D_max(j);
        plot_and_save_vox_max(comparison_name, cohen_d_struct, 'cohen_D', source, source_forward, cutoff_point_D, voxel_index, condition_dir)

    end
end
