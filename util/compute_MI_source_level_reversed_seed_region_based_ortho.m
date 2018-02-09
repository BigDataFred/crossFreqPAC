function compute_MI_source_level_reversed_seed_region_based_ortho(mode,files,seedIx,hemi)


%% set path defs
[path2files] = '/home/rouxf/prj/TC/matFilesRev/';

%%
for zt = 1:size(files,2)
         
    %%  loading the data & some prepcoessing  
    if mode ==1
        
        lfDat = load([path2files,files(1,zt).name]);
        hfDat = load([path2files,files(2,zt).name]);        
        
        lfDat = lfDat.VC;
        hfDat = hfDat.VC;                
    
        %hilbert transform the virtual channel data
        cfg = [];
        cfg.channel = hfDat.label(seedIx);%extract hfDat amp for the reversed-seed grid points
        cfg.hilbert             = 'abs'; % compute the envelope
        
        [hfDat] = ft_preprocessing( cfg , hfDat );
    else
        dat = load([path2files,files(zt).name]);
        
        % DO NOT use hilbert transform here
        cfg = [];
        cfg.channel = {'all'};
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [9 11];
        cfg.bpfilttype = 'fir';
        cfg.padding = max(dat.VC.time{1});%
        cfg.padtype = 'mirror';
        cfg.continuous = 'no';
        
        [lfDat] = ft_preprocessing(cfg, dat.VC);% single trial lfDat activity
        
        % DO use hilbert transform here
        cfg = [];
        cfg.channel = dat.VC.label(seedIx);%extract hfDat amp for the reversed-seed grid points
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [60 80];
        cfg.bpfilttype = 'fir';
        cfg.padding = max(dat.VC.time{1});%
        cfg.padtype = 'mirror';
        cfg.continuous = 'no';
        cfg.hilbert             = 'abs'; % compute the envelope
        
        [hfDat] = ft_preprocessing(cfg,dat.VC);
        clear dat;
    end;    
       
    %% start parallel pool
    if isempty(gcp('nocreate'))
       parpool(36,'SpmdEnabled', false);
    end;        
    
    %% reverse seed voxel orthogonalization
    
    %idx = setdiff(1:length(lfDat.label),seedIx);%do orthogonalization for all except ref channel
    idx = 1:length(lfDat.label);%do orthogonalization for all except ref channel
    
    for nt = 1:length(lfDat.trial) 
        fprintf([num2str(nt),'/',num2str(length(lfDat.trial))]);
        X = lfDat.trial{nt}(seedIx,:);
        Y = lfDat.trial{nt}(idx,:);
        
        lfDat.trial{nt} = zeros(size(lfDat.trial{nt}));
        
        parfor xt = 1:length(idx)   
            Yo = zeros(size(X,1),size(X,2));
            for yt = 1:size(X,1)
                [Yo(yt,:)] = orthogonalize_time_domain( X(yt,:) , Y(xt,:) );            
            end;
            if any(any(isnan(Yo)))
                error('NaNs detected');
            end;
            Y(xt,:) = mean(Yo);
        end;      
        
        %lfDat.trial{nt}(seedIx,:) = X;% keep the original
        lfDat.trial{nt}(idx,:) = Y;% replace the original with the orthogonalized part
         
        clear X Y:
        fprintf('\n');
    end;
    %delete(gcp);
    
    %% hilbert transform to obtain phase time series on orthogonalized data
    cfg = [];
    cfg.channel = {'all'};%extract alpha phase
    cfg.hilbert = 'angle';% phase in rad
    cfg.continuous = 'no';
    
    [lfDat]                 = ft_preprocessing(cfg,lfDat);
    
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
    
%     if sign(size(concat2,1)-1) == 1
%         error('number of channels is out of range');
%     end;
    
    %% phase bins
    [pbins] = -pi:pi/8:pi;       
    
    %% loop over virtual channels
    [amp] = concat2;% extract the seed region hfDat power for ind VC
    [PAH] = zeros(size(concat1,1),length(pbins));
    
    for it = 1:size(concat1,1)%loop over virtual channels
        
        fprintf([num2str(it),'/',num2str(size(concat1,1))]);        
        [phi] = concat1(it,:);% extract the whole brain lfDat phase     
        
        [X] = zeros(1,length(pbins));
        parfor kt = 1:length(pbins)-1% loop over phase bins     
            amp;pbins;
            X(kt) = mean(mean(amp(:,phi  >= pbins(kt) & phi < pbins(kt+1))));           
        end;
        X(end) = X(1);        
        PAH(it,:) = X;
        
        fprintf('\n'); 
    end;
    clear amp;
    
    %% close parallel pool
    clear concat*;
    
    %% compute MI
    PAH = PAH./repmat(sum(PAH,2),[1 size(PAH,2)]);
    H = -sum(log(PAH).*PAH,2);
    n = length(pbins);
    MI = (log(n)-H)./log(n);
    
    %% save data
    savepath = path2files;
    if mode ==1
        save([savepath,'whole_brain_alphaGammaPAC_reversed_seed_based_ortho_',num2str(length(seedIx)),'_',hemi,'H_',files(1,zt).name],'MI','pbins','PAH');
    else
        save([savepath,'whole_brain_alphaGammaPAC_reversed_seed_based_ortho_',num2str(length(seedIx)),'_',hemi,'H_',files(zt).name],'MI','pbins','PAH');
    end;

end;

%% clear and exit
clear;
%exit;