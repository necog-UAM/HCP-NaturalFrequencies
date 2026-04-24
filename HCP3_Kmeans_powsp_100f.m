% function HCP3_Kmeans_powsp (substsk, task, p, cfg)
% 
% % task: 'resting' for resting, 'motor' for motor ......
% % subs: list of subjects; e.g., subs = [100307, 102816, 559053];
% %
% % For example, HCP3_Kmeans_powsp ('resting') prepares data from all subjects
% % and sessions of task to subsequetly perform Kmeans clustering in HCP4
% 
% Nsub = length(substsk);
% rawpath  = p.rawpath;   
% datapath = p.datapath; 
% kmeanspath = p.kmeanspath;
% task    = cfg.task;
% numtask = cfg.numtask;
% condit  = cfg.condit;
% condit_label = cfg.condit_label;
% condit_code  = cfg.condit_code;
% % source = cfg.source.inverse;
% % source_forward = cfg.source.forward;
% %%  3.1. Selection of trials per condition
% 
% Ntr      = 100;               % select only 100 trials per subject to reduce computational load
% 
% if strcmp(task,'Restin')
%     Ncondit = 1;
% elseif strcmp(task,'StoryM')
%     Ncondit = 2;
% elseif strcmp(task,'Wrkmem')
%     Ncondit = 5; %4
% elseif strcmp(task,'Motort')
%     Ncondit = 5;
% end
% 
% %% control
% cd(kmeanspath)
% if exist([kmeanspath '\' sprintf('kmeans_%s_%d_powsp',task,condit)])==2
%     eval(sprintf('load kmeans_%s_%d_powsp',task,condit));
% else
% 
%     for c = 1:Ncondit
%         powsptotc{c} = [];
%         ksubc{c}     = [];
%         kvoxc{c}     = [];
%         ktrialc{c}   = [];
%     end
% 
%     seed = rng(sum(clock),'twister');
% 
%     for s = 1:Nsub
%         display(['Subject ' num2str(s) '/' num2str(Nsub)])
% 
%         sub      = num2str(substsk{s});
%         route = fullfile([datapath sub '\_' task]);
%         if exist(route)
%             cd(route)
%             load freqsource_100f %freqsource
% 
%             if strcmp(task,'StoryM')
%                 powsp2{1} = powsp(:,:,trialinfo==1);               % 1: Story
%                 powsp2{2} = powsp(:,:,trialinfo==2);               % 2: Math
%             elseif strcmp(task,'Wrkmem')
%                 powsp2{1} = powsp(:,:,trialinfo(:,2)==1);          % 1: 0-Back
%                 powsp2{2} = powsp(:,:,trialinfo(:,2)==2);          % 2: 2-Back
%                 powsp2{3} = powsp(:,:,trialinfo(:,1)==1);          % 3: Face
%                 powsp2{4} = powsp(:,:,trialinfo(:,1)==2);          % 4: Tools
%                 powsp2{5} = powsp(:,:,trialinfo(:,1)==0);        % 5: Fixation   % N/A
%             elseif strcmp(task,'Motort')
%                 powsp2{1} = powsp(:,:,trialinfo==1);          % 1: Left Hand
%                 powsp2{2} = powsp(:,:,trialinfo==2);          % 2: Left Foot
%                 powsp2{3} = powsp(:,:,trialinfo==4);          % 3: Right Hand
%                 powsp2{4} = powsp(:,:,trialinfo==5);          % 4: Right Foot
%                 powsp2{5} = powsp(:,:,trialinfo==6);          % 5: Fixation
%             else
%                 powsp2{1} = powsp;
%             end
% 
%             % Count how many trials per condition
%             for c = 1:Ncondit
%                 ntrials(s,c) = size(powsp2{c},3);
%             end
% 
%             for c = 1:Ncondit
%                 if size(powsp2{c},3) > Ntr
%                     rndtr    = randperm(size(powsp2{c},3),Ntr);
%                 else
%                     Ntr2 = 43;
%                     rndtr    = randperm(size(powsp2{c},3),Ntr2); % We need to change this section and adjust Ntr2
%                 end
%                 powsp2{c}   = powsp2{c}(:,:,rndtr);
% 
%                 ct = 1;
%                 powsp3{c}  = single(zeros(size(powsp2{c},1)*size(powsp2{c},3),size(powsp2{c},2)));
%                 ksub2{c}   = [];
%                 kvox2{c}   = [];
%                 ktrial2{c} = [];
%                 for i = 1:size(powsp2{c},1)
%                     for tr = 1:size(powsp2{c},3)
%                         powsp3{c}(ct,:) = powsp2{c}(i,:,tr);
%                         ksub2{c}(ct)    = s;
%                         kvox2{c}(ct)    = i;
%                         ktrial2{c}(ct)  = rndtr(tr);
%                         ct = ct+1;
%                     end
%                 end
%                 powsptotc{c} = [powsptotc{c};powsp3{c}];
%                 ksubc{c}     = [ksubc{c} ksub2{c}];
%                 kvoxc{c}     = [kvoxc{c} kvox2{c}];
%                 ktrialc{c}   = [ktrialc{c} ktrial2{c}];
%             end
%         end
% 
%         clear powsp powsp2 powsp3 bl
%         for c = 1:Ncondit
%             bl = sum(powsptotc{c},2);                   % mean power across time for each frequency
%             powsptotc{c} = powsptotc{c}./repmat(bl,[1 size(powsptotc{c},2)]);
%         end
% 
%         cd(kmeanspath)
%         eval(sprintf('save ntrials_%s ntrials',task));
% 
%         for c = 1:Ncondit
%             powsptot = powsptotc{c};
%             ksub     = ksubc{c};
%             kvox     = kvoxc{c};
%             ktrial   = ktrialc{c};
%             eval(sprintf('save kmeans_%s_%d_powsp powsptot ksub kvox ktrial seed -v7.3',task,c));
%         end
%     end
% end


function HCP3_Kmeans_powsp_OMEGA (substsk, task, p, cfg, OMEGA, sess)

% task: 'resting' for resting, 'motor' for motor ......
% subs: list of subjects; e.g., subs = [100307, 102816, 559053];
%
% For example, HCP3_Kmeans_powsp ('resting') prepares data from all subjects
% and sessions of task to subsequetly perform Kmeans clustering in HCP4
if nargin < 5
    OMEGA = 1;
    sess = false;
end
Nsub = length(substsk);
rawpath  = p.rawpath;   
datapath = p.datapath; 
kmeanspath = p.kmeanspath;
task    = cfg.task;
numtask = cfg.numtask;
condit  = cfg.condit;
condit_label = cfg.condit_label;
condit_code  = cfg.condit_code;
% source = cfg.source.inverse;
% source_forward = cfg.source.forward;
%%  3.1. Selection of trials per condition

Ntr      = 100;               % select only 100 trials per subject to reduce computational load

if strcmp(task,'Restin')
    Ncondit = 1;
elseif strcmp(task,'StoryM')
    Ncondit = 2;
elseif strcmp(task,'Wrkmem')
    Ncondit = 5; %4
elseif strcmp(task,'Motort')
    Ncondit = 5;
end

%% control
cd(kmeanspath)
if exist([kmeanspath '\' sprintf('kmeans_%s_%d_powsp',task,condit)])==2
    eval(sprintf('load kmeans_%s_%d_powsp',task,condit));
else

    for c = 1:Ncondit
        powsptotc{c} = [];
        ksubc{c}     = [];
        kvoxc{c}     = [];
        ktrialc{c}   = [];
    end

    seed = rng(sum(clock),'twister');

    for s = 1:Nsub
        display(['Subject ' num2str(s) '/' num2str(Nsub)])

        sub      = num2str(substsk{s});
        if OMEGA == 0
            ses = char(sess(s));
            route = fullfile([datapath 'sub-' sub '\ses-' ses]);
        else
            route = fullfile([datapath sub '\_' task]);
        end
        if exist(route)
            cd(route)
            if OMEGA == 0
                load freq_allvox_10mm 
            else
                load freqsource_100f %freqsource
            end

            if strcmp(task,'StoryM')
                powsp2{1} = powsp(:,:,trialinfo==1);               % 1: Story
                powsp2{2} = powsp(:,:,trialinfo==2);               % 2: Math
            elseif strcmp(task,'Wrkmem')
                powsp2{1} = powsp(:,:,trialinfo(:,2)==1);          % 1: 0-Back
                powsp2{2} = powsp(:,:,trialinfo(:,2)==2);          % 2: 2-Back
                powsp2{3} = powsp(:,:,trialinfo(:,1)==1);          % 3: Face
                powsp2{4} = powsp(:,:,trialinfo(:,1)==2);          % 4: Tools
                powsp2{5} = powsp(:,:,trialinfo(:,1)==0);        % 5: Fixation   % N/A
            elseif strcmp(task,'Motort')
                powsp2{1} = powsp(:,:,trialinfo==1);          % 1: Left Hand
                powsp2{2} = powsp(:,:,trialinfo==2);          % 2: Left Foot
                powsp2{3} = powsp(:,:,trialinfo==4);          % 3: Right Hand
                powsp2{4} = powsp(:,:,trialinfo==5);          % 4: Right Foot
                powsp2{5} = powsp(:,:,trialinfo==6);          % 5: Fixation
            else
                powsp2{1} = powsp;
            end

            % Count how many trials per condition
            for c = 1:Ncondit
                ntrials(s,c) = size(powsp2{c},3);
            end

            for c = 1:Ncondit
                if size(powsp2{c},3) > Ntr
                    rndtr    = randperm(size(powsp2{c},3),Ntr);
                else
                    Ntr2 = 43;
                    rndtr    = randperm(size(powsp2{c},3),Ntr2); % We need to change this section and adjust Ntr2
                end
                powsp2{c}   = powsp2{c}(:,:,rndtr);

                ct = 1;
                powsp3{c}  = single(zeros(size(powsp2{c},1)*size(powsp2{c},3),size(powsp2{c},2)));
                ksub2{c}   = [];
                kvox2{c}   = [];
                ktrial2{c} = [];
                for i = 1:size(powsp2{c},1)
                    for tr = 1:size(powsp2{c},3)
                        powsp3{c}(ct,:) = powsp2{c}(i,:,tr);
                        ksub2{c}(ct)    = s;
                        kvox2{c}(ct)    = i;
                        ktrial2{c}(ct)  = rndtr(tr);
                        ct = ct+1;
                    end
                end
                powsptotc{c} = [powsptotc{c};powsp3{c}];
                ksubc{c}     = [ksubc{c} ksub2{c}];
                kvoxc{c}     = [kvoxc{c} kvox2{c}];
                ktrialc{c}   = [ktrialc{c} ktrial2{c}];
            end
        end

        clear powsp powsp2 powsp3 bl
        for c = 1:Ncondit
            bl = sum(powsptotc{c},2);                   % mean power across time for each frequency
            powsptotc{c} = powsptotc{c}./repmat(bl,[1 size(powsptotc{c},2)]);
        end

        cd(kmeanspath)
        eval(sprintf('save ntrials_%s ntrials',task));

        for c = 1:Ncondit
            powsptot = powsptotc{c};
            ksub     = ksubc{c};
            kvox     = kvoxc{c};
            ktrial   = ktrialc{c};
            eval(sprintf('save kmeans_%s_%d_powsp powsptot ksub kvox ktrial seed -v7.3',task,c));
        end
    end
end