function HCP1_Source_reconstruction(sub, task, ses, p)

% sub: subject code (e.g.,102816)
% task: 'resting' for resting, 'motor' for motor ......
% ses: ses number (e.g., 3, 4, or 5 for resting)
% p: path information
% 
% For example, HCP1_Source_reconstruction (102816, 'resting', 3, p) prepares
% source_forward and source_inverse for subject 102816, session 3 of resting 

sub      = num2str(sub);
rawpath = p.rawpath;   
datapath = p.datapath; 

%% 2.1. Coregistration of MEG-MRI spaces
% mri.transform is the transformation matrix to go from mri space to sensor space
    
cd([rawpath sub '\MEG\anatomy\'])
addpath(genpath('C:\toolbox\fieldtrip-20210912\'))                   % addpath gives problems in ft_read_mri
mri  =  ft_read_mri('T1w_acpc_dc_restore.nii');    % open .TAR.GZ file into MRIcroN and save it as .nii to read it here
rmpath(genpath('C:\toolbox\fieldtrip-20210912\'))                   
addpath('C:\toolbox\fieldtrip-20210912\')                 

fid  = fopen([sub '_MEG_anatomy_transform.txt']);
txt  = textscan(fid,'%s');
txt  = txt{1};
fclose(fid);

for tt = 1:length(txt)
    if strcmp(txt{tt},'transform.vox07mm2bti')
        t1 = tt;
    end
end

transform = [];
for i = 1:16
    transform(i) = str2num(txt{t1+i+2});
end
mri.transform = reshape(transform,[4 4])';

cd(datapath)
mkdir(sub)
cd([datapath sub])

%% 2.2. MRI normalization
% mri normalized (mrin) and transformation matrix (normtrans)
% check coordsys!: yes, a, l, s, i

cd([datapath sub])
if exist([datapath sub '\' 'mri.mat'])==2
    load mri
else
    cfg            = [];
    cfg.nonlinear  = 'no';
    cfg.spmversion = 'spm8';
    mrin           = ft_volumenormalise(cfg, mri);
    
    % determine the affine source->template coordinate transformation (from fieldtrip-20180405) 
    normtrans = mrin.params.VG.mat * inv(mrin.params.Affine) * inv(mrin.params.VF.mat) * mrin.initial;

    save mri mri mrin normtrans
end

%% 2.3. Head model
% semi-realistic singleshell head model based on the implementation from Guido Nolte
% check coordsys!: yes, a, l, s, i

if exist([datapath sub '\' 'vol.mat'])==2
    load vol
else
    cfg             = [];
    cfg.spmversion  = 'spm8';
    segment         = ft_volumesegment(cfg,mri);
    segment.anatomy = mri.anatomy;
    
    cfg        = [];
    cfg.method = 'singleshell';
    vol        = ft_prepare_headmodel(cfg, segment);
    save vol vol
end

%% 2.4. Forward model
% output: grid2 (normalized grid and leadfields)

% Read data to extract grad information
if strcmp(task,'Restin')
    cd([rawpath sub '\MEG\' 'Restin\rmegpreproc'])
    a = ls('*Restin*');
    for i = 1:size(a,1)
        temp = strfind(a(i,:),['MEG_' num2str(ses)]);
        if temp > 1
            disp(['Reading ' a(i,:)])
            load(a(i,:))
        end
    end
elseif strcmp(task,'StoryM')
    cd([rawpath sub '\MEG\' 'StoryM\tmegpreproc'])
    a = ls('*StoryM_tmegpreproc_BU*');
    for i = 1:size(a,1)
        temp = strfind(a(i,:),['MEG_' num2str(ses) '-StoryM_tmegpreproc_BU' ]);  
        if temp > 1
            disp(['Reading ' a(i,:)])
            load(a(i,:))
        end
    end
elseif strcmp(task,'Wrkmem')
    cd([rawpath sub '\MEG\' 'Wrkmem\tmegpreproc'])
    a = ls('*Wrkmem_tmegpreproc_TIM*');
    for i = 1:size(a,1)
        temp = strfind(a(i,:),['MEG_' num2str(ses) '-Wrkmem_tmegpreproc_TIM']);   
        if temp > 1
            disp(['Reading ' a(i,:)])
            load(a(i,:))
        end
    end
elseif strcmp(task,'Motort')
    cd([rawpath sub '\MEG\' 'Motort\tmegpreproc'])
    a = ls('*Motort_tmegpreproc_TEMG*');
    for i = 1:size(a,1)
        temp = strfind(a(i,:),['MEG_' num2str(ses) '-Motort_tmegpreproc_TEMG']);    
        if temp > 1
            disp(['Reading ' a(i,:)])
            load(a(i,:))
        end
    end
end

cfg         = [];
cfg.channel = 'MEG';               % remove EMG channels in motor task
data        = ft_preprocessing(cfg,data);

grad = data.grad;


% Load normalized template grid (10mm) (in fieldtrip/template/sourcemodel)
load standard_sourcemodel3d10mm         
grid = sourcemodel; 

% Load normalized mri and head model
cd([datapath sub])
load mri
load vol

% Adapt the normalized grid to each individual's brain space
posmni    = grid.pos;
pos       = ft_warp_apply(inv(normtrans), grid.pos*10, 'homogenous')/10;
grid.pos  = pos;
grid.unit = 'cm';

% Convert grad, vol and grid to common units (mm)
grad = ft_convert_units(grad, vol.unit);
grid = ft_convert_units(grid, vol.unit);

% Select only voxels within cortical mask (e.g., cerebellum is excluded)
cd(datapath)
load ('correccion_vox_inside_10mm.mat')
grid.inside = inside;

% Compute leadfields for each grid's voxel
cfg             = [];
cfg.grid        = grid;
cfg.grad        = grad;
cfg.vol         = vol;
cfg.channel     = {'MEG'};
cfg.normalize   = 'yes';
cfg.reducerank  = 2;
grid2           = ft_prepare_leadfield(cfg);

% Check that grad, vol and grid are coregistered
plot3 (grad.chanpos(:,1), grad.chanpos(:,2), grad.chanpos(:,3), '.','MarkerEdgeColor',[0.8 0 0],'MarkerSize',25), hold on
plot3 (vol.bnd.pos(:,1), vol.bnd.pos(:,2), vol.bnd.pos(:,3), '.','MarkerEdgeColor',[0 0 0.8]), hold on
plot3 (grid2.pos(grid2.inside,1), grid2.pos(grid2.inside,2), grid2.pos(grid2.inside,3), '+k'), hold off

% Save grad, vol, grid and mri in source_forward structure to be used later
source_forward      = [];
source_forward.vol  = vol;
source_forward.mri  = mrin;
source_forward.grad = grad;
source_forward.grid = grid2;

cd([datapath sub])
mkdir([task num2str(ses)])
cd([task num2str(ses)])
savefig('coreg')
save source_forward source_forward   

%% 2.5. Source reconstruction
% source_forward and source_inverse should be computed specifically for each ses (different grad coordinates)

cfg              = [];
cfg.covariance   = 'yes';
cfg.vartrllength = 2;
datacov          = ft_timelockanalysis(cfg, data);       % covariance matrix

% Compute spatial filter (in source.avg.filter)
cfg                   = [];
cfg.method            = 'lcmv';
cfg.grad              = source_forward.grad;
cfg.headmodel         = source_forward.vol;
cfg.grid              = source_forward.grid;
cfg.lcmv.fixedori     = 'yes';
cfg.lcmv.normalize    = 'yes';
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.keepfilter   = 'yes';          % save filter for using it later
cfg.lcmv.lambda       = '10%';          
cfg.lcmv.reducerank   = 2;
source                = ft_sourceanalysis(cfg, datacov);

load standard_sourcemodel3d10mm
source.avg.ori = {};
source.avg.mom = {};
source.avg.noisecov = {};
source.pos     = sourcemodel.pos;               % standard grid positions

cd(datapath)
load ('correccion_vox_inside_10mm.mat')
source.inside=inside;

cd([datapath sub])
mkdir([task num2str(ses)])
cd([task num2str(ses)])
save source_inverse source


