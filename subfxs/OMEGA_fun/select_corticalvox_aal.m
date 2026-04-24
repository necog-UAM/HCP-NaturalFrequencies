%% Create cortical voxels based on AAL atlas
addpath(genpath('I:\Matems\toolbox_MEG\spm8'))

cd('I:\Mis documentos\PROYECTOS\OMEGA\scripts')
aal = spm_vol('aal.nii');           % this is to change to radiological convention (the T1_big_LR has been modified manually to improve it, so it is not necessary to invert it again, only to smooth it)

load aal_voxel_label

dat = aal.private.dat(:,:,:);
dat_aal = zeros(size(dat));     % keep aal label
dat_mask = zeros(size(dat));    % 0/1 mask 

aal_removed = [0,41,42,71:78,91:116];      % 91:116 is cerebellum            

for i = 1:size(dat,1)
    for j = 1:size(dat,2)
        for k = 1:size(dat,3)
            if sum(aal_removed == dat(i,j,k)) == 0
                dat_aal(i,j,k) = dat(i,j,k);
                dat_mask(i,j,k) = 1;
            else
                dat_aal(i,j,k) = 0;
                dat_mask(i,j,k) = 0;
            end
        end
    end
end

aal.fname = 'aal_cortical_subcortical.nii';
spm_write_vol(aal,dat_aal)

aal.fname = 'aal_cortical_subcortical_mask.nii';
spm_write_vol(aal,dat_mask)

% addpath spm8 
spm_reslice({'T1h.img','aal_cortical_subcortical.nii'})               % images in 91x109x91 scale
spm_reslice({'T1h.img','aal_cortical_subcortical_mask.nii'})


%% Load data
sub = '0001';
ses = '0001';
dpath = 'I:\Mis documentos\PROYECTOS\OMEGA\';              % folder for processed data
cd([dpath 'OMEGA_data\\sub-' sub '\ses-' ses])
load source_forward_10mm
load source_inverse_10mm

%% Check voxels-areas
source2 = source;
voxel_inside = find(source.inside==1);
nvox = length(voxel_inside);
source2.avg.pow(voxel_inside) = [1:length(voxel_inside)];

cfg              = [];
cfg.parameter    = 'pow';
cfg.downsample   = 2;
cfg.interpmethod = 'nearest';
temp             = ft_sourceinterpolate(cfg, source2, source_forward.mri);
template         = reshape(temp.pow,[91 109 91]);      

cd('I:\Mis documentos\PROYECTOS\OMEGA\scripts')
% aal = spm_vol('aal.nii');           % this is to change to radiological convention (the T1_big_LR has been modified manually to improve it, so it is not necessary to invert it again, only to smooth it
% aalLR = zeros(aal.dim);
% Ni = aal.dim(1);
% 
% change L/R to radiological convention
% for i = 1:Ni
%     aalLR(i,:,:) = aal.private.dat(Ni-i+1,:,:);
% end
% aal.fname = 'aal_LR.nii';
% spm_write_vol(aal,aalLR)

% now L/R is in radiological convention: left is right

% aal = ft_read_atlas('scortical_mask_T1_big_LR2.nii');    % LR2 is 91x109x91
aal = ft_read_atlas('raal_cortical_mask.nii');    % same for kmeans than for ICA [big for kmeans; normal for ICA}
aal_voxel = NaN(size(template));
aal_label = NaN(size(template));

for i = 1:size(template,1)
    for j = 1:size(template,2)
        for k = 1:size(template,3)      
            if aal.brick0(i,j,k) == 1 
                aal_voxel(i,j,k) = template(i,j,k);
            end
        end
    end
end

voxel_inside_aal = [];
label_inside_aal = [];
Nvox = [];

for i = 1:nvox
    Nvox(i) = sum(aal_voxel(:) == i);
end

[N,ind] = sort(Nvox);            
% ind = ind (N>8);            % for ICA analysis: only MEG voxels that have at least 3*3*3 voxel inside (125 = 5*5*5 is completely inside) (RM voxel size 2 mm; MEG voxel size 10 mm)
ind = ind (N>8);            % for kmeans 10mm (1 for kmeans 12mm)

% source2.avg.pow(voxel_inside(ind))
source_forward.grid.inside = zeros(length(source.inside),1);
source_forward.grid.inside(voxel_inside(ind))=1;
source_forward.grid.inside=source_forward.grid.inside==1;
inside = source_forward.grid.inside;

% load correccion_vox_inside_10mm
% inside(voxin)=1;         % these should be inside as well

cd('I:\Mis documentos\PROYECTOS\OMEGA\scripts')
save aal_cortical_mask_inside_10mm inside

rmpath(genpath('I:\Matems\toolbox_MEG\spm8'))
