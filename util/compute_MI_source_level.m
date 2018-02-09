function compute_MI_source_level(mode,files)

%% set path defs
[path2files] = '/home/rouxf/prj/TC/matFilesRev/';

%%
for zt = 1:size(files,2)
        
    %% hilbert transform the virtual channel data    
    if mode ==1
        
        lfDat = load([path2files,files(1,zt).name]);
        hfDat = load([path2files,files(2,zt).name]);
        
        cfg = [];        
        cfg.hilbert = 'angle';% phase in rad
        cfg.continuous = 'no';
        
        [lfDat] = ft_preprocessing(cfg, lfDat.VC);% single trial lfDat activity        
        
        cfg = [];        
        cfg.hilbert = 'abs'; % compute the power        
        cfg.continuous = 'no';
        
        [hfDat] = ft_preprocessing(cfg,hfDat.VC);
        
    else
        
        dat = load([path2files,files(zt).name]);
        
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [9 11];
        cfg.bpfilttype = 'fir';
        cfg.hilbert = 'angle';% phase in rad
        cfg.padding = max( dat.VC.time{1});%
        cfg.padtype = 'mirror';
        cfg.continuous = 'no';
        
        [lfDat] = ft_preprocessing(cfg, dat.VC);% single trial lfDat activity
                
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [60 80];
        cfg.bpfilttype = 'fir';
        cfg.hilbert = 'abs'; % compute the power
        cfg.padding = max( dat.VC.time{1});%
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
    
    %% phase bins
    [pbins] = -pi:pi/8:pi;
    %preallocate
    [PAH] = zeros(size(concat1,1),length(pbins));
    
    %%
    if isempty(gcp('nocreate'))
        parpool(36,'SpmdEnabled', false);
    end;
    
    %% loop over virtual channels
    parfor it = 1:size(concat1,1)
        
        [phi] = concat1(it,:);
        [amp] = concat2(it,:);
        X = zeros(length(pbins),1);
        for kt = 1:length(pbins)-1% loop over phase bins                        
            X(kt) = mean(amp(phi  >= pbins(kt) & phi < pbins(kt+1)));            
        end;
        X(end) = X(1);
        PAH(it,:) = X;
        
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
        save([savepath,'whole_brain_alphaGammaPAC_local_',files(1,zt).name],'MI','pbins','PAH');
    else
        save([savepath,'whole_brain_alphaGammaPAC_local_',files(zt).name],'MI','pbins','PAH');
    end;
    
end;

%% clear and exit
clear all;
%exit;



















