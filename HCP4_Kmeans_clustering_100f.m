function HCP4_Kmeans_clustering (substsk, task, condit, p)

% task: e.g., 'Restin' for resting
% condit: condition within task
%
% For example, HCP4_Kmeans_clustering ('StoryM', 1) performs Kmeans clustering 
% for task 'StoryM' 1st condition (Story) with data from all subjects pooled together

Nsub = length(substsk);
rawpath  = p.rawpath;   
datapath = p.datapath; 
kmeanspath = p.kmeanspath; 

Nk   = 25;
Nvox = 1925;
if (task == 'Wrkmem' & condit == 5)
    Ntr = 43;
else
    Ntr  = 100;
end

%%  4.1. K-means clustering
     
cd(kmeanspath)
if exist([kmeanspath '\' sprintf('kmeans_%s_%d_NK%d_3c.mat',task,condit, Nk)])==2
    eval(sprintf('load kmeans_%s_%d_NK%d_3c',task,condit, Nk));
    eval(sprintf('load kmeans_%s_%d_powsp',task,condit));
else
    eval(sprintf('load kmeans_%s_%d_powsp',task,condit));
     
    f1         = 0.55:0.05:0.9;          % frequencies are limited by trial length (1.2 s -> min freq 2.6 Hz with 3 cycles; lower frequencies analyzed with 2 up to 3 cycles)
    foi1      = exp(f1);
    t_ftimwin1 = [2:0.125:2.875]./foi1;             
    
    f2         = 0.95:0.05:4.6; %3.55;         
    foi2      = exp(f2);
    t_ftimwin2 = 3./foi2;                 % time-window = 3 cicles of corresponding frequency
    
    foi = [foi1 foi2];
    t_ftimwin = [t_ftimwin1 t_ftimwin2];
    
    % % new!!!!! only 3 cycles
    % f         = 0.95:0.05:4.6;          % frequencies are limited by trial length (1.2 s -> min freq 2.6 Hz with 3 cycles; lower frequencies analyzed with 2 up to 3 cycles)
    % foi       = exp(f);
    % powsptot  = powsptot(:,9:end);
    % % --------------------- %
    
    [idx,C,sumd,D]  = kmeans(powsptot,Nk,'Distance','cosine','Display','iter','Replicates',5,'MaxIter',200);
    
    
    
    % check how many segments each subject has
    % sustituir por lo de Lydia:
    propk = NaN(Nk,Nvox);
        for k = 1:Nk
            disp(['Cluster ' num2str(k) '/' num2str(Nk)])
            ctk = find(idx==k); % I
            ctkvox = kvox(ctk)';
            ctkvoxfilt = ctkvox==[1:Nvox];
            propk(k,:)=sum(ctkvoxfilt)./Ntr; % Aquí hay un error. Se divide por Ntr = 100; pero creo que no debería ser así. Ntr debería ser el total (ya que para working memory fixation hay menos de 100 tr disponibles).
        end
        
       
    
    %--------------------------- Antiguo código -------------------------------
    % propk=NaN(Nk,Nvox,Nsub);
    % for k=1:Nk
    %     disp(['Cluster ' num2str(k) '/' num2str(Nk)])
    %     ctk = find(idx==k);
    %     ctksub = ksub(ctk)';
    %     ctkvox = kvox(ctk)';
    % 
    %     ctksubfilt = ctksub==[1:Nsub];
    %     ctkvoxfilt = ctkvox==[1:Nvox];
    %     for s=1:Nsub
    %         propk(k,:,s)=sum(ctksubfilt(:,s) & ctkvoxfilt)./Ntr;
    %     end
    % end 
    %--------------------------------------------------------------------------
    
    cd(kmeanspath)
    %eval(sprintf('save kmeans_%s_%d_Nk%d_3c idx C sumd D Nvox propk -v7.3', task,condit,Nk, Ntr));
    eval(sprintf('save(''kmeans_%s_%d_Nk%d_3c.mat'', ''idx'', ''C'', ''sumd'', ''D'', ''Nvox'', ''propk'', ''-v7.3'')', task, condit, Nk));
end
