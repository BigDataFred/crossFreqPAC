function compute_MI_source_level_reversed_seed_region_based(mode,files,seedIx,hemi)

%% set path defs
[path2files] = '/home/rouxf/prj/TC/matFilesRev/';

%%
for zt = 1:size(files,2)

    %% hilbert transform the virtual channel data    
    if mode ==1
        lfDat = load([path2files,files(1,zt).name]);
        hfDat = load([path2files,files(2,zt).name]);
        
        %low frequency data
        cfg = [];
        cfg.channel = {'all'};%extract alpha phase       
        cfg.hilbert = 'angle';% phase in rad        
        cfg.continuous = 'no';
        
        [lfDat] = ft_preprocessing(cfg, lfDat.VC);% single trial lfDat activity        
        
        %high frequency data
        cfg = [];        
        cfg.channel = hfDat.VC.label(seedIx);%extract hfDat amp
        cfg.hilbert = 'abs'; % compute the power        
        cfg.continuous = 'no';
        
        [hfDat] = ft_preprocessing(cfg,hfDat.VC);
    else
        dat = load([path2files,files(zt).name]);
        
        %low frequency data
        cfg = [];
        cfg.channel = {'all'};
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [9 11];
        cfg.bpfilttype = 'fir';
        cfg.hilbert = 'angle';% phase in rad
        cfg.padding = max(dat.VC.time{1});%
        cfg.padtype = 'mirror';
        cfg.continuous = 'no';
        
        [lfDat] = ft_preprocessing(cfg, dat.VC);% single trial lfDat activity
            
        %high frequency data    
        cfg = [];
        cfg.channel =dat.VC.label(seedIx);%extract hfDat amp
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [60 80];
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
    %preallocate
    [PAH] = zeros(size(concat2,1),size(concat1,1),length(pbins));
    
    %% loop over virtual channels    
    %for it = 1:size(concat2,1)%loop over virtual channels
        %fprintf([num2str(it),'/',num2str(size(concat2,1))]);
        
        %[amp] = concat2(it,:);% extract the seed region hfDat power for ind VC
        [amp] = mean(concat2,1);
        [PAH] = zeros(1,size(concat1,1),length(pbins));
        [dum] = zeros(size(concat1,1),length(pbins));
        for jt = 1:size(concat1,1)
            fprintf([num2str(jt),'/',num2str(size(concat1,1))]);
            amp2=amp;
            [phi] = concat1(jt,:);% extract the whole brain lfDat phase             
            
            [X] = zeros(1,length(pbins),1);
            parfor kt = 1:length(pbins)-1% loop over phase bins                                
                X(kt) = mean(amp2(phi  >= pbins(kt) & phi < pbins(kt+1)));                
            end;
            X(end) = X(1);            
            dum(jt,:) = X;   
            fprintf('\n');
        end;   
        PAH(1,:,:) = dum;
        clear amp phi dum;
        %PAH(it,:,:) = dum;
        %fprintf('\n');
    %end;
    
    %% close parallel pool
    clear concat*;
    %delete(gcp);
    
    %% compute MI
    PAH = PAH./repmat(sum(PAH,3),[1 1 size(PAH,3)]);
    H = -sum(log(PAH).*PAH,3);
    n = length(pbins);
    MI = (log(n)-H)./log(n);
    
    %% save data
    savepath = path2files;
    if mode ==1
        save([savepath,'whole_brain_alphaGammaPAC_reversed_seed_based_',num2str(length(seedIx)),'_',hemi,'H_',files(1,zt).name],'MI','pbins','PAH');
    else
        save([savepath,'whole_brain_alphaGammaPAC_reversed_seed_based_',num2str(length(seedIx)),'_',hemi,'H_',files(zt).name],'MI','pbins','PAH');
    end;
    
end;

%% clear and exit
clear all;
%exit;
