function HCP5_Naturalfreq_singlesubject (cfg, substsk, task, condit, p, sub)

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

cfg2 = cfg; % Since cfg will be overwritten later

voxel_inside = find(source.inside==1);

%% 4.1. Reconstruction of source-level time-series
    
%sub = deblank(subs(s,:));
cd([datapath sub '\' ])
% Control de condiciones
a = ls('*_*');
if length(strmatch(['_' task],a)) > 0
    cd(['_' task])
    load freqsource_100f %freqsource % Contiene trialinfo
    cd(kmeanspath)
    if exist([kmeanspath '\' sprintf('kmeans_%s_%d_Nk%d_3c.mat',task,condit,Nk)])==2
        eval(sprintf('load kmeans_%s_%d_Nk%d_3c',task,condit,Nk));

        if strmatch(condit_label{t}, 'Restin')
            powspsubc = powsp(:,:,:);
        elseif ismember(condit_label{t}{t}, {'0-back', '2-back'})
            powspsubc = powsp(:,:,trialinfo(:,2) == condit_code(condit));
        elseif ismember(condit_label{t}{t}, {'Face', 'Tools', 'WM_Fixation'})
            powspsubc = powsp(:,:,trialinfo(:,1) == condit_code(condit));
        else
            powspsubc = powsp(:,:,trialinfo == condit_code(condit));
        end
        ct = 1;
        powspsubc2  = single(zeros(size(powspsubc,1)*size(powspsubc,3),size(powspsubc,2)));
        kvox   = [];
        ktrial = [];
        for i = 1:size(powspsubc,1)
            for tr = 1:size(powspsubc,3)
                powspsubc2(ct,:) = powspsubc(i,:,tr);
                kvox(ct)    = i;
                ktrial(ct)  = tr;
                ct = ct+1;
            end
        end

        bl       = sum(powspsubc2,2);                                % save kmeans_10mm_baseline bl
        powspsubc2 = powspsubc2./repmat(bl,[1 size(powspsubc2,2)]);      % compute relative power to correct the center of the head bias

        [D2,I] = pdist2(C,powspsubc2, 'cosine','Smallest',1); % D2 porque ya existe D.

    % 6.2. Proportion of power spectra categorized in each cluster

        if strmatch(condit_label{t}, 'Restin')
            Ntr = length(trialinfo);
        elseif ismember(condit_label{t}{t}, {'0-back', '2-back'})
            trial_logic = trialinfo(:,2) == condit_code(condit);
            Ntr = length(trial_logic(trial_logic==1)); % Number of trials
        elseif ismember(condit_label{t}{t}, {'Face', 'Tools', 'WM_Fixation'})
            trial_logic = trialinfo(:,1) == condit_code(condit);
            Ntr = length(trial_logic(trial_logic==1)); % Number of trials
        else
            trial_logic = trialinfo == condit_code(condit);
            Ntr = length(trial_logic(trial_logic==1)); % Number of trials
        end

        propk = NaN(Nk,Nvox);
        for k = 1:Nk
            disp(['Cluster ' num2str(k) '/' num2str(Nk)])
            ctk = find(I==k);
            ctkvox = kvox(ctk)';
            ctkvoxfilt = ctkvox==[1:Nvox];
            propk(k,:)=sum(ctkvoxfilt)./Ntr;
        end
        
        label = condit_label{numtask}{condit};
        cd([datapath sub '\_' task])
        fil = sprintf('save kmeans_10mm_%s_%d_Nk%d_%dtr D D2 I C propk -v7.3',task,condit,Nk,Ntr)
        eval(fil)

    %% 7.1. Compute natural frequency of each voxel (output: natfreq)

        f   = 0.55:0.05:4.6;%3.55;
        foi = exp(f);

%         f   = 0.9:0.05:3.5;
%         fx = exp(f);
%         foi2 = interp(fx(1:end-1),10);

        ff  = [];
        for k = 1:Nk
            sp = C(k,:);
            [pks,locs] = findpeaks(sp);
            fx  = round(foi(locs(find(pks == max(pks)))),1);
            ff(k) = fx;
        end
        [ffsort,idf]=sort(ff); 
        ffsort = unique(ffsort);

        ff = NaN(Nk,2);
        badk = [];
        for k = 1:Nk
            sp = C(idf(k),:);
            [pks,locs] = findpeaks(sp,'MinPeakHeight',0.1);
            % if it does not find any peak, go to the next centroid
            if length(pks) == 1                                % unimodal spectrum
                ff(idf(k),1) = round(foi(locs(1)),1);
                ff(idf(k),2) = round(foi(locs(1)),1);      % if unimodal spectrum, repeat 1st peak
            elseif length(pks) == 2                           % bimodal spectrum
                if (foi(locs(1)-1)*2 > foi(locs(2)-1) &  foi(locs(1)-1)*2 < foi(locs(2)+1)) | ...
                        (foi(locs(1)+1)*2 > foi(locs(2)-1) &  foi(locs(1)+1)*2 < foi(locs(2)+1))
                    ff(idf(k),1) = round(foi(locs(1)),1);     % harmonics: only 1st peak (fundamental frequency)
                    ff(idf(k),2) = round(foi(locs(1)),1);
                else  
                    ff(idf(k),1) = round(foi(locs(1)),1);
                    ff(idf(k),2) = round(foi(locs(2)),1);
                end
            elseif length(pks) > 2                         % in case of a third residual peak
                [~,id] = sort(pks,'descend');
                locs = locs(id(1:2));
                if (foi(locs(1)-1)*2 > foi(locs(2)-1) &  foi(locs(1)-1)*2 < foi(locs(2)+1)) | ...
                        (foi(locs(1)+1)*2 > foi(locs(2)-1) &  foi(locs(1)+1)*2 < foi(locs(2)+1))
                    ff(idf(k),1) = round(foi(locs(1)),1);     % harmonics: only 1st peak (fundamental frequency)
                    ff(idf(k),2) = round(foi(locs(1)),1);
                else     
                    ff(idf(k),1) = round(foi(locs(1)),1);     % no harmonics: take both peaks
                    ff(idf(k),2) = round(foi(locs(2)),1);
                end
            elseif length(pks) == 0
                badk = [badk idf(k)];
            end
        end
    
        propk2 = zeros(length(ffsort),Nvox);

        for i=1:length(ffsort)
            [r,c] = find(ff==ffsort(i));
            r = unique(r);
            for j=1:length(r)
                if ff(r(j),1) ~= ff(r(j),2)
                    propk2(i,:) = propk2(i,:) + 1/2.*sum(propk(r(j),:),1);    % bimodal
                elseif ff(r(j),1) == ff(r(j),2)
                    propk2(i,:) = propk2(i,:) + sum(propk(r(j),:),1);    % unimodal
                end
            end
        end

        emptycol = find(sum(propk2,2)==0);
        propk2(emptycol,:) = [];
        ffsort(emptycol) = [];

        propkz = (propk2-mean(propk2,2))./std(propk2,0,2);        % z-normalize
    
        [dim,xx,yy,zz,connmat, dtempl] = Omega_neighbors_ly(source); % pers function
    
        propkz2sm = NaN(size(propkz));
        f2 = interp(ffsort(1:end-1),10);
        fnatsm=[];
        fnatsm_T=[];
        fnatsm_p=[];
        fnatsm_powsp = zeros(Nvox,length(foi));
    
        for vx=1:Nvox
            vxneigh = connmat(vx,:)==1;
            propkzneig= propkz(:,vxneigh);
            propinterp = [];
            for i=1:size(propkzneig,2)
                propinterp(:,i) = interp1(ffsort,propkzneig(:,i),foi,'pchip');
            end
            [h,p,ci,stats] = ttest(propinterp');
    
            [tmax,idmax] = max(stats.tstat);
            fnatsm(vx) = foi(idmax);
            fnatsm_T(vx) = tmax;
            fnatsm_p(vx) = p(idmax);
    
            tval = stats.tstat;
            tval(p>.05) = NaN;
            tval(tval < 0) = NaN;
            [pks,locs] = findpeaks(tval);
    
            fnat.pks{vx}=foi(locs);
            fnat.w{vx}=pks*100./sum(pks);

            fnatsm_powsp(vx,:) = mean(propinterp,2);    
        end
    
        fnatsm2 = fnatsm;
        fnatsm2(fnatsm_p>.05) = NaN;
    
        fnat.fnat = fnatsm;
        fnat.fnatsig = fnatsm2;
        fnat.tval = fnatsm_T;
        fnat.pval = fnatsm_p;
        fnat.powsp = fnatsm_powsp;

        cd([datapath sub '\_' task])
        fil = sprintf('save singlesub_fnat_%s_%d_new fnat propkz -v7.3',task,condit)
        eval(fil) 
    
        % Plot single-subject natural frequency maps
    
        source2 = source;
        source2.avg.pow(voxel_inside) = log(fnatsm2);
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
        cfg.funcolorlim   = [0.7 3.4];
        cfg.funcolormap   = 'jet_omega_mod';
        cfg.projmethod    = 'nearest';
        cfg.opacity       = 0.8;
        cfg.camlight      = 'no';
        cfg.colorbar      = 'no';
        cfg.surffile     = 'surface_pial_left.mat';
        cfg.surfinflated  = 'surface_inflated_left_caret_white.mat';
        subplot(2,2,1), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('left')
        subplot(2,2,3), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('left')
    
        cfg.surffile     = 'surface_pial_right.mat';
        cfg.surfinflated  = 'surface_inflated_right_caret_white.mat';
        subplot(2,2,2), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('right')
        subplot(2,2,4), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('right')
        
        if ~exist([outpath '\' task '\'  condit_label{numtask}{condit}], 'dir')
            mkdir([outpath '\' task '\'  condit_label{numtask}{condit}]);
        end
        
        if strmatch(condit_label{t}, 'Restin')
            cd([outpath '\' task '\'])
        else
            cd([outpath '\' task '\'  condit_label{numtask}{condit}])
        end
        current_dir = pwd
        if ~exist([current_dir 'singlesub'], 'dir')
            mkdir('singlesub');
        end
        mkdir('singlesub');
        cd('singlesub')
        print('-dtiff','-r300',['singlesub_fnat_' sub '_new.tiff']);
        close
    end
end

end

