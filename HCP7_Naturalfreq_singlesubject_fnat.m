function [fnat_sub, fnat_powsp] = HCP7_Naturalfreq_singlesubject_fnat (cfg, substsk, task, condit, p, sub, OMEGA, ses)

if nargin < 6
    OMEGA = 1;
    sess = false;
end

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
%source = cfg.source.inverse;
%source_forward = cfg.source.forward;

cfg2 = cfg;

%voxel_inside = find(source.inside==1);
%%
if OMEGA == 0
    route = fullfile([datapath 'sub-' sub '\ses-' ses]);
else
    route = fullfile([datapath sub '\_' task]); % cambiar 1 por indice de task (t)
end
if exist(route)
    cd(route)
    fil = sprintf('load singlesub_fnat_%s_%d_new',task,condit);
    eval(fil)
    fnat_sub = fnat.fnatsig;
    fnat_powsp = fnat.powsp;
end
end