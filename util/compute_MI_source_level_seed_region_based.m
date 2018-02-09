function compute_MI_source_level_seed_region_based(mode,files,seedIx,hemi)

%% set path defs
[path2files] = '/home/rouxf/prj/TC/matFilesRev/';

%%
for zt = 1:size(files,2)
        
    %% hilbert transform the virtual channel data    
    if mode ==1
        lfDat = load([path2files,files(1,zt).name]);
        hfDat = load([path2files,files(2,zt).name]);
        
        cfg = [];
        cfg.channel = lfDat.VC.label(seedIx);%extract alpha phase       
        cfg.hilbert = 'angle';% phase in rad        
        cfg.continuous = 'no';
        
        [lfDat] = ft_preprocessing(cfg, lfDat.VC);% single trial lfDat activity        
        
        cfg = [];        
        cfg.channel = {'all'};
        cfg.hilbert = 'abs'; % compute the power        
        cfg.continuous = 'no';
        
        [hfDat] = ft_preprocessing(cfg,hfDat.VC);
    else
        dat = load([path2files,files(zt).name]);
        
        cfg = [];
        cfg.channel = dat.VC.label(seedIx);%extract alpha phase
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [9 11];
        cfg.bpfilttype = 'fir';
        cfg.hilbert = 'angle';% phase in rad
        cfg.padding = max(dat.VC.time{1});%
        cfg.padtype = 'mirror';
        cfg.continuous = 'no';
        
        [lfDat] = ft_preprocessing(cfg, dat.VC);% single trial lfDat activity
                
        cfg = [];
        cfg.channel = {'all'};
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [65 75];
        cfg.bpfilttype = 'fir';
        cfg.hilbert = 'abs'; % compute the power
        cfg.padding = max(dat.VC.time{1});%
        cfg.padtype = 'mirror';
        cfg.continuous = 'no';
        
        [hfDat] = ft_preprocessing(cfg,dat.VC);
        clear dat;
    end;
    
    %% concatenate
    [concat1] = zeros(length(lfDat.label),length(lfDat.time{1})*length(lfDat.trial));
    [concat2] = zeros(length(hfDat.label),length(hfDat.time{1})*length(hfDat.trial));
    
    idx = 1:length(lfDat.time{1});
    for it = 1:length(lfDat.trial)
        concat1(:,idx) = lfDat.trial{it};
        concat2(:,idx) = hfDat.trial{it}.^2;% square the amp to obtain power
        idx = idx + length(lfDat.time{it});
    end;
    clear lfDat hfDat;
    
    %% start parallel pool
    if isempty(gcp('nocreate'))
       parpool(36,'SpmdEnabled', false);
    end;
        
    %% phase bins
    [pbins] = -pi:pi/8:pi;
    
    [phi] = concat1;clear concat1;% extract the lfDat phase from seed region
    
    %preallocate
    [PAH] = zeros(size(concat2,1),length(pbins));
     
    for it = 1:size(concat2,1)%loop over virtual channels
        fprintf([num2str(it),'/',num2str(size(concat2,1))]);
        [amp] = concat2(it,:);% extract whole brain hfDat power for ind VC
        
        X = zeros(length(pbins),1);        
        parfor kt = 1:length(pbins)-1% loop over phase bins  
            amp2 = amp;
            pbins2 =pbins;
            tmp = zeros(1,size(phi,1));
            for yt = 1:size(phi,1)
                [idx] = find(phi(yt,:)  >= pbins2(kt) & phi(yt,:) < pbins2(kt+1));            
                 tmp(yt)= mean(amp2(idx));  
            end;
            X(kt) = mean(tmp);
        end;
        X(end) = X(1);
        
        PAH(it,:) = X;
        fprintf('\n');
    end;
    
    %% close parallel pool
    clear concat*;
    %delete(gcp);
    
    %% compute MI
    PAH = PAH./repmat(sum(PAH,2),[1 size(PAH,2)]);
    H = -sum(log(PAH).*PAH,2);
    n = length(pbins);
    MI = (log(n)-H)./log(n);
    
    %% save data
    savepath = path2files;
    if mode ==1
        save([savepath,'whole_brain_alphaGammaPAC_seed_based_',num2str(length(seedIx)),'_',hemi,'H_',files(1,zt).name],'MI','pbins','PAH');
    else
        save([savepath,'whole_brain_alphaGammaPAC_seed_based_',num2str(length(seedIx)),'_',hemi,'H_',files(zt).name],'MI','pbins','PAH');
    end;
    
end;
%% clear and exit
clear;
%exit;
