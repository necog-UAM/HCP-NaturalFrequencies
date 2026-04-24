function data2 = HCP2_Freqsource (sub, task, p, onlysource)

% sub: subject code (e.g.,102816)
% task: 'resting' for resting, 'motor' for motor ......
% p: path information
%
% For example, HCP2_Freqsource (102816, 'resting', p) computes power
% spectrum (freq_allvox) in source-space (1 result per 1.2-s-interval, adapting 
% the length of the time window to each frequency) 

sub      = num2str(sub);
rawpath  = p.rawpath;   
datapath = p.datapath; 

%% 2.1. Read data, source_forward and source_inverse

datas = {};
% if strcmp(task,'Restin')
%     cd([rawpath sub '\MEG\' 'Restin\rmegpreproc'])
%     a = ls('*Restin*');
%     for i = 1:size(a,1)
%         disp(['Reading ' a(i,:)])
%         load(a(i,:))
%         datas{i} = data;
%     end
% elseif strcmp(task,'StoryM')
%     cd([rawpath sub '\MEG\' 'StoryM\tmegpreproc'])
%     a = ls('*StoryM_tmegpreproc_BU*');
%     for i = 1:size(a,1)
%         disp(['Reading ' a(i,:)])
%         load(a(i,:))
%         datas{i} = data;
%     end
% elseif strcmp(task,'Wrkmem')
%     cd([rawpath sub '\MEG\' 'Wrkmem\tmegpreproc'])
%     a = ls('*Wrkmem_tmegpreproc_TIM*');
%     for i = 1:size(a,1)
%         disp(['Reading ' a(i,:)])
%         load(a(i,:))
%         datas{i} = data;
%     end
% elseif strcmp(task,'Motort')
%     cd([rawpath sub '\MEG\' 'Motort\tmegpreproc'])
%     a = ls('*Motort_tmegpreproc_TEMG*');
%     for i = 1:size(a,1)
%         disp(['Reading ' a(i,:)])
%         load(a(i,:))
%         datas{i} = data;
%     end
% end
% Define the base path and task-specific subdirectories
base_dir = [rawpath sub '\MEG\'];

% Use a switch-case structure for better readability
switch task
    case 'Restin'
        subdir = 'Restin\rmegpreproc';
        pattern = '*Restin*';
    case 'StoryM'
        subdir = 'StoryM\tmegpreproc';
        pattern = '*StoryM_tmegpreproc_BU*';
    case 'Wrkmem'
        subdir = 'Wrkmem\tmegpreproc';
        pattern = '*Wrkmem_tmegpreproc_TIM*';
    case 'Motort'
        subdir = 'Motort\tmegpreproc';
        pattern = '*Motort_tmegpreproc_TEMG*';
    otherwise
        error('Invalid task specified.');
end

% Construct the full directory path
full_dir = fullfile(base_dir, subdir);

% Try to change to the directory and list files
try
    % Check if the directory exists before attempting to cd
    if ~isfolder(full_dir)
        warning('Directory does not exist: %s', full_dir);
        datas = {}; % Return an empty cell array if the directory does not exist
        return;
    end
    
    % Change to the directory
    cd(full_dir);
    
    % List files matching the pattern
    a = ls(pattern);
    
    % Check if any files were found
    if isempty(a)
        warning('No files found in directory: %s', full_dir);
        datas = {}; % Return an empty cell array if no files are found
        return;
    end
    
    % Load data from each file
    datas = {}; % Initialize the data cell array
    for i = 1:size(a, 1)
        disp(['Reading ' a(i,:)]);
        load(a(i,:)); % Load the file
        datas{end+1} = data; % Store the loaded data in the cell array
    end

catch ME
    % Handle unexpected errors
    warning('An unexpected error occurred: %s', ME.message);
    datas = {}; % Return an empty cell array on error
end
%%
cd([datapath sub])
if strcmp(task,'Restin')
    a = ls('Restin*');
elseif strcmp(task,'StoryM')
    a = ls('StoryM*');
elseif strcmp(task,'Wrkmem')
    a = ls('Wrkmem*');
elseif strcmp(task,'Motort')
    a = ls('Motort*');
end

for i = 1:size(a,1)
    cd([datapath sub '\' a(i,:)])
    load source_forward
    load source_inverse
    sourcef{i} = source_forward;
    sourcei{i} = source;
end


%% 2.2. Reconstruct source-level data

tmax = 1.2;
Nses  = length(datas);
datax = data;             % source reconstructed data
datax.trialinfo = [];
ct = 1;

for ses = 1:Nses
    disp(['Source reconstruction: Subject ' sub ' - Session ' num2str(ses) '/' num2str(Nses)])
    data = datas{ses};

    cfg         = [];
    cfg.channel = 'MEG';               % remove EMG channels in motor task
    data        = ft_preprocessing(cfg,data);

    if strcmp(task,'Restin') | strcmp(task,'Wrkmem') | strcmp(task,'Motort')
        cfg        = [];
        cfg.toilim = [0.002 tmax-0.001]; %%
        data      = ft_redefinetrial(cfg,data);

    elseif strcmp(task,'StoryM')                 % split entire blocks into 1.2-s trials
        for i = 1:20                             % up to 20 seconds
            cfg        = [];
            cfg.toilim = [0.002 tmax-0.001] + tmax*(i-1); %%
            datai{i}   = ft_redefinetrial(cfg,data);
            if size(datai{i}.time{1},2) == 609
                cfg        = [];
                cfg.toilim = [0.001 tmax-0.001] + tmax*(i-1);
                datai{i}   = ft_redefinetrial(cfg,data);
            end
        end
        data = ft_appenddata([],datai{:});
    end
    
    trls = [];
    trlnan = [];
    for i = 1:length(data.trial)
        trls(i) = size(data.trial{i},2);   
        trlnan(i) = sum(sum(isnan(data.trial{i})));
    end

    cfg = [];
    cfg.trials = trls==610 & trlnan==0;           % remove shorter trials and trials with NaNs
    data       = ft_redefinetrial(cfg,data);

    source_forward = sourcef{ses};
    source = sourcei{ses};

    load ('aal_voxel_label')
    voxel_inside = find(source.inside==1);
    time         = data.time{1};

    Nvox = length(voxel_inside);
    Ntime = length(time);
    Ntrial = length(data.trial);

    datasource = zeros(Nvox,Ntime,Ntrial);
    addpath(genpath('Z:\Toolbox\fieldtrip-20230118\'))                   % addpath gives problems in ft_read_mri
    [~, selch] = match_str(data.label, source.avg.label);
    rmpath(genpath('Z:\Toolbox\fieldtrip-20230118\'))                   
    addpath('Z:\Toolbox\fieldtrip-20230118')                 

    for i = 1:Nvox
        sourcefilt = source.avg.filter{voxel_inside(i)}(selch);
        for tr = 1:Ntrial
            datasource(i,:,tr) = sourcefilt * data.trial{tr};
        end
    end
    datasource = datasource./repmat(std(datasource,0,2),[1,Ntime,1]);    % correction of centre of the head bias
        
    datax.label     = {};
    for i = 1:Nvox
        datax.label{i} = ['V' num2str(i)];
    end

    for tr = 1:Ntrial
        datax.trial{ct} = single(datasource(:,:,tr));
        if strcmp(task,'StoryM')
            datax.trialinfo(ct,1) = data.trialinfo(tr,2);     % 1: Story, 2: Math
        elseif strcmp(task,'Wrkmem')
            datax.trialinfo(ct,1) = data.trialinfo(tr,4);     % imgType: 1- Face, 2- Tools  0- Fixation
            datax.trialinfo(ct,2) = data.trialinfo(tr,5);     % memoryType :  1: 0-Back   2: 2-Back
        elseif strcmp(task,'Motort')
            datax.trialinfo(ct,1) = data.trialinfo(tr,2);     % 1-Left Hand,  2 - Left Foot, 4 - Right Hand. 5 - Right Foot, 6 - Fixation
        else
            datax.trialinfo(ct,1) = NaN;
        end
     
        datax.time{ct}  = time;
        ct = ct+1;
    end
    clear datasource
end

Ntrial     = length(datax.trial);
trialinfo  = datax.trialinfo;

% Change data structure to continuous mode to speed freqanalysis up and for padding
data2          = datax;
data2.trial    = {};
data2.trial{1} = [];
data2.time     = {};
data2.time{1}  = [];
data2.sampleinfo = [1 length(data2.trial{1})];

for tr = 1:Ntrial
    data2.trial{1} = [data2.trial{1} datax.trial{tr}];
    data2.time{1}  = [data2.time{1} datax.time{tr}+tmax*(tr-1)];
end
clear datax

%% 2.3. Frequency analysis

% Frequency analysis parameters: foi, toi, t_ftimwin

if onlysource == 0
    f1         = 0.55:0.05:0.9;          % frequencies are limited by trial length (1.2 s -> min freq 2.6 Hz with 3 cycles; lower frequencies analyzed with 2 up to 3 cycles)
    foi1      = exp(f1);
    t_ftimwin1 = [2:0.125:2.875]./foi1;             
    
    f2         = 0.95:0.05:4.6;     %3.55    
    foi2      = exp(f2);
    t_ftimwin2 = 3./foi2;                 % time-window = 3 cicles of corresponding frequency
    
    foi = [foi1 foi2];
    t_ftimwin = [t_ftimwin1 t_ftimwin2];
    
    % Frequency analysis computation
    cfg              = [];
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.output       = 'pow';
    cfg.foi          = foi;
    cfg.toi          = tmax/2:tmax:data2.time{1}(end);    % one spectrum per trial (every 1.2-s), at the midpoint of each trial (0.6 s)
    cfg.t_ftimwin    = t_ftimwin;
    cfg.keeptrials   = 'yes';
    cfg.pad          ='nextpow2';
    freq             = ft_freqanalysis(cfg, data2);
    
    powsp = single(squeeze(freq.powspctrm));
    
    cd([datapath sub])
    mkdir(['_' task])
    cd(['_' task])
    save freqsource_100f powsp foi trialinfo
end
