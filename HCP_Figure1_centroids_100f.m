function HCP_Figure1_centroids (cfg, p)

% outpath: output folder
% figpath: folder to save localization of generators in cortical surface (F1.4.)
% replicate sample: 'tot' (N=128), 'rep1' (1st replicate sample N/2), 'rep2' (2nd replicate sample N/2)
% cfg: specification of Nk (number of clusters), Nvox (number of voxels), and Nsub (number of subjects)
% cfg.fig_powsp: 1 to plot centroid power spectra 
% cfg.stats: 1 to compute stats;
% cfg.fig_brainpatt: 1 to plot brain distribution associated to each cluster
% cfg.locgenerators_vox: list of peak voxels generating different brain rhythms
% cfg.locgenerators_fig: 1 to plot the localization of generators in cortical surface
% cfg.fig_gamma: 1 to plot (artefactual) generators of gamma-band activity

% -------------------------------------------------------- %
% F2.1. Plot centroid power spectra for each cluster (1-30Hz)
% F2.2. Brain pattern of each cluster: stats
% F2.3. Brain pattern of each cluster: figure (masked with T-stat)
% F2.4. Localize generators in cortical surface
% F2.5. Supplementary Figure 1: plot high-frequency (>30 Hz) generators in slices
% -------------------------------------------------------- %

%% F2.1. Plot centroid power spectra for each cluster (1-30Hz)

Nk   = cfg.Nk;
Nsub = cfg.Nsub;
Nvox = cfg.Nvox;
task = cfg.task;
numtask = cfg.numtask;
condit = cfg.condit;
condit_label = cfg.condit_label;
cfg2 = cfg;

outpath = p.outpath;
figpath = p.figpath; 
kmeanspath = p.kmeanspath; 

cd(kmeanspath)
eval(sprintf('load kmeans_%s_%d_Nk%d_3c',cfg.task,cfg.condit,Nk));   

f         = 0.55:0.05:4.6;%3.55;          % frequencies are limited by trial length (1.2 s -> min freq 2.6 Hz with 3 cycles; lower frequencies analyzed with 2 up to 3 cycles)
foi      = exp(f);
ff  = [];
for k = 1:Nk
    sp = C(k,:);
    [pks,locs] = findpeaks(sp);
    fx  = round(foi(locs(find(pks == max(pks)))),1);
    % f  = round(foi(find(sp==max(sp)))*10)/10;
    ff(k) = fx;
end
[~,idf]=sort(ff);

% Nk2   = sum(ff < 30);   % select only clusters with peak centroids <30Hz (>30Hz are most probably artefacts) 
% foi   = foi(1:findbin(foi,35));
c     = jet_omega_mod; % jet_omega_mod is 64x3 --> This causes the error when foi is 1x82.
fqlog =  exp(0.65:0.048:4.55);%exp(0.7:0.044:3.5);

cd(outpath)
mkdir(task)
cd(task)
if strcmp(task,'Restin')==0
    mkdir(condit_label{numtask}{condit})
    cd(condit_label{numtask}{condit})
end

for k = 1:Nk
    % sp = C(idf(k),1:findbin(foi,35));
    sp = C(idf(k),:);
    [pks,locs] = findpeaks(sp,'MinPeakHeight',0.1);

    if length(pks) > 0
        f = round(foi(locs),1);
        fc = [];
        for i = 1:length(locs)
            fc(i) = findbin(fqlog,f(i));
        end
        % Clamp fc values to the valid range [1, 64]
        fc = min(max(fc, 1), 64);
        if cfg2.fig_powsp == 1
            figure('Color', [1 1 1])
            set(gcf,'Position',[100 100 850 350])

            plot(foi,sp,'Color',c(fc(find(pks == max(pks))),:),'LineWidth',6) % error is caused here
            if length(f) == 1
                title(['k ' num2str(k),' - ',num2str(f),' Hz'])
            else
                title(['k ' num2str(k),' - ',num2str(f(1)),' / ',num2str(f(2)),' Hz'])
            end
            set (gca,'XLim',[1.5 90],'YLim',[0 .4],'FontName','Calibri','FontSize',30,'XTick',[10:10:90],'YTickLabel',{'','','','',''})
            box(gca,'off');
            print('-dtiff','-r300',['Nk' num2str(Nk) '_Power_Spectrum_k' num2str(k),' - ',num2str(f(1)), 'Hz.tiff']);
        end
    end
end


%% F2.2. Brain pattern of each cluster: stats

rng('shuffle')
if cfg2.stats == 1
    cd(kmeanspath)
    load source_forward
    load source_inverse
    voxel_inside = find(source.inside==1);

    propkz = (propk-mean(propk,2))./std(propk,0,2);        % Z-normalize
    propkzm = mean(propkz,3);

    cd(outpath)
    for k = 1:Nk
        disp(['Cluster ' num2str(k) '/' num2str(Nk)])
        d_obs  = squeeze(propkz(k,:,:))';
        d_obsm = repmat(mean(d_obs,2),[1,Nvox]);
        [h,p,ci,stats] = ttest(d_obs,d_obsm);
        t_obs(k,:) = stats.tstat;
        for it = 1:1000
            for s = 1:Nsub
                vox_shuf = randperm(Nvox);
                dx = squeeze(propkz(k,vox_shuf,s))';
                d_shuf(s,:)= dx;
            end
            d_shufm = repmat(mean(d_shuf,2),[1,Nvox]);
            [h,p,ci,stats] = ttest(d_shuf,d_shufm);
            t_max = max(stats.tstat);
            t_shuf(it) = t_max;
        end
        t_thr05(k)  = prctile(t_shuf,95);
        t_thr01(k)  = prctile(t_shuf,99);
        t_thr001(k) = prctile(t_shuf,99.9);
    end
    
    cd([outpath '\' task '\'  condit_label{numtask}{condit}])
    save kmeans_stats_DEF idf t_obs t_thr05 t_thr01 t_thr001
end


%% F2.3. Brain pattern of each cluster: figure (masked with T-stat)

if cfg2.fig_brainpatt == 1
    cd(kmeanspath)
    load source_forward
    load source_inverse
    voxel_inside = find(source.inside==1);
    
    propkz = (propk-mean(propk,2))./std(propk,0,2);        % Z-normalize
    propkzm = mean(propkz,3);
    
    if strcmp(task,'Restin')==0
        cd([outpath '\' task '\'  condit_label{numtask}{condit}])
    else
        cd([outpath '\' task ])
    end
    % load kmeans_stats_DEF
    
    for k = 1:Nk
        sp = C(idf(k),:);
        [pks,locs] = findpeaks(sp,'MinPeakHeight',0.1);
        if length(pks) > 0
            if length(locs) > 1
                locs = locs(find(pks == max(pks)));
            end
            f = round(foi(locs)*10)/10;

            source2 = source;
            source2.avg = rmfield(source2.avg,'mom');
            % source2.avg.pow(voxel_inside) = propkzm(idf(k),:).*t_sig;        % Z-value masked by stat
            source2.avg.pow(voxel_inside) = propkzm(idf(k),:);

            cfg = [];
            cfg.parameter  = 'avg.pow';
            cfg.downsample = 2;
            cfg.interpmethod  =  'linear';
            source_interp  = ft_sourceinterpolate (cfg, source2, source_forward.mri);

            % Write MRIcroN image
            cfg = [];
            cfg.filetype  = 'nifti';
            cfg.parameter = 'pow';
            cfg.filename  = ['Nk' num2str(Nk) '_k' num2str(k) '_' num2str(round(f)) 'Hz'];
            ft_sourcewrite(cfg, source_interp)

            % Plot and save results in cortical surface
            figure('Color',[1 1 1]);
            set(gcf,'Position',[100 100 1100 900])

            cfg               = [];
            cfg.method        = 'surface';
            cfg.funparameter  = 'pow';
            cfg.maskparameter = cfg.funparameter;
            cfg.funcolorlim   = [0  2];
            cfg.funcolormap   = 'hot_omega_mod';
            cfg.projmethod    = 'nearest';
            cfg.opacitymap    = 'rampup';
            cfg.camlight      = 'no';
            cfg.colorbar      = 'yes';
            cfg.figure        = 'gcf';
            cfg.surffile      = 'surface_pial_left.mat';
            cfg.surfinflated  = 'surface_inflated_left_caret_white.mat';
            subplot(2,2,1), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('left')
            subplot(2,2,3), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('left')

            cfg.surffile      = 'surface_pial_right.mat';
            cfg.surfinflated  = 'surface_inflated_right_caret_white.mat';
            subplot(2,2,2), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('right')
            subplot(2,2,4), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('right')

            print('-dtiff','-r300',['Nk' num2str(Nk) '_k' num2str(k) '_' num2str(round(f)) 'Hz.tiff']);
        end
    end
end


%% F2.4. Localize generators in cortical surface

if cfg2.locgenerators_fig == 1
    cd(figpath)
    vx = cfg2.locgenerators_vox;
    for i = 1:length(vx)
        v(i) = find(voxel_inside==vx(i));
    end
    
    for i = 1:length(v)
        figure('Color',[1 1 1]);
        set(gcf,'Position',[100 100 1100 900])
        
        source2 = source;
        source2.avg.pow(voxel_inside) = 0;        
        source2.avg.pow(voxel_inside(v(i))) = 1;    
        
        cfg=[];
        cfg.parameter  = 'avg.pow';
        cfg.downsample = 2;
        cfg.interpmethod  = 'nearest';
        source_interp = ft_sourceinterpolate (cfg, source2, source_forward.mri);
        
        cfg               = [];
        cfg.method        = 'surface';
        cfg.funparameter  = 'pow';
        cfg.maskparameter = cfg.funparameter;
        cfg.funcolormap   = 'hot_omega_mod';   % parula
        cfg.projmethod    = 'nearest';
        cfg.opacitymap    = 'rampup';
        cfg.camlight      = 'no';
        cfg.colorbar      = 'yes';
        cfg.surffile      = 'surface_pial_left.mat';
        cfg.surfinflated  = 'surface_inflated_left_caret_white.mat';
        subplot(2,2,1), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('left')
        subplot(2,2,3), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('left')
        
        cfg.surffile      = 'surface_pial_right.mat';
        cfg.surfinflated  = 'surface_inflated_right_caret_white.mat';
        subplot(2,2,2), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('right')
        subplot(2,2,4), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('right')
        
        print('-dtiff','-r300',['AAL_voxel' num2str(vx(i)) '.tiff']);
        close(gcf)
    end
end


%% F2.5. Supplementary Figure 1: plot high-frequency (>30 Hz) generators in slices

if cfg2.fig_gamma == 1
    
    Nk   = cfg.Nk;
    Nsub = cfg.Nsub;
    Nvox = cfg.Nvox;
    cfg2 = cfg;
    
    cd([outpath  'Nk25_10mm_' rep])
    fil = sprintf('load kmeans_10mm_Nk%d_%s',Nk,rep);
    eval(fil)
    
    f   = 0.55:0.05:4.6;
    foi = exp(f);
    ff  = [];
    for k = 1:Nk
        sp = C(k,:);
        f  = round(foi(find(sp==max(sp)))*10)/10;
        ff(k) = f;
    end
    [~,idf]=sort(ff);
    
    cd(outpath)
    load source_forward_10mm
    load source_inverse_10mm
    voxel_inside = find(source.inside==1);
    
    propkz = (propk-mean(propk,2))./std(propk,0,2);        % Z-normalize
    propkzm = mean(propkz,3);
    
    cd([outpath  'Nk25_10mm_' rep])
    load kmeans_stats_DEF
    
    for k = 20:Nk
        t_sig = t_obs(idf(k),:);
        t_sig(t_sig < t_thr05(idf(k)))  = 0;
        t_sig(t_sig >= t_thr05(idf(k))) = 1;
        
        sp = [C(idf(k),:) 0];
        [pks,locs] = findpeaks(sp,'MinPeakHeight',0.1);
        if length(locs) == 2          % if mu rhythm, then assign the frequency of the beta peak
            locs = locs(2);
        end
        f = round(foi(locs)*10)/10;
        
        source2 = source;
        source2.avg.pow(voxel_inside) = propkzm(idf(k),:).*t_sig;        % Z-value masked by stat
        
        temp = source2.avg.pow(voxel_inside);                 % check minimum Z-value significant
        temp(temp ==0) = NaN;
        zthr(k)=nanmin(temp);
        
        cfg = [];
        cfg.parameter  = 'avg.pow';
        cfg.downsample = 2;
        cfg.interpmethod  =  'linear';
        source_interp  = ft_sourceinterpolate (cfg, source2, source_forward.mri);
        
        figure('Color',[1 1 1]);
        set(gcf,'Position',[100 100 1100 900])
        cfg               = [];
        cfg.method        = 'slice';
        cfg.nslices       = 10;
        cfg.funparameter  = 'pow';
        cfg.funcolorlim   = [0 2.5];
        cfg.maskparameter = cfg.funparameter;
        cfg.funcolormap   = 'hot';
        cfg.projmethod    = 'nearest';
        % cfg.interactive   = 'yes';
        cfg.inputcoord    = 'mni';
        ft_sourceplot(cfg,source_interp)
        
        print('-dtiff','-r300',['SupplFig1_Nk' num2str(Nk) '_k' num2str(k) '_' num2str(round(f)) 'Hz.tiff']);
        
    end
end

