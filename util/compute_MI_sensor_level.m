function compute_MI_sensor_level(files)

%%
p2df = '~rouxf/prj/TC/matFiles/';
files = dir([p2df,files]);

%%
for yt = 1:length(files)
    
    load([p2df,files(yt).name]);% load data
    
    %%
    for xt = 1:length(save_data1)
        
        raw1 = save_data1{xt};
        
        %% format the data into rpt_chan_time format
        cfg = [];
        cfg.keeptrials = 'yes';
        
        [raw1] = ft_timelockanalysis(cfg,raw1);
        
        %% parameters for PAC analysis
        % frequencies of interest
        pfoi = 4:2:20;
        afoi = 45:5:140;
        % phase bins
        [pbins] = -pi:pi/8:pi;
        
        %% preallocate
        [PAH] = zeros(length(raw1.label),length(pfoi),length(afoi),length(pbins));
        
        %% PAC for MEG channels
        for it = 1:length(pfoi) % loop over frequency for phase
            cfg = [];
            cfg.bpfilter = 'yes';
            cfg.bpfreq = [pfoi(it)-1 pfoi(it)+1];
            cfg.bpfilttype = 'fir';
            cfg.hilbert = 'angle';% phase in rad
            cfg.padding = 3*max(raw1.time);%
            cfg.padtype = 'mirror';
            cfg.continuous = 'no';
            cfg.channel = raw1.label;
            
            [phi] = ft_preprocessing(cfg,raw1);% single trial alpha activity
            concat1 = zeros(length(phi.label),length(phi.time)*size(phi.trial,1));
            idx = 1:length(phi.time);
            for zt = 1:size(phi.trial,1)
                concat1(:,idx) = squeeze(phi.trial(zt,:,:));
                idx = idx + length(phi.time);
            end;
            clear phi;
            
            for jt = 1:length(afoi) % loop over frequency for amplitude
                cfg = [];
                cfg.bpfilter = 'yes';
                cfg.bpfreq = [afoi(jt)-5 afoi(jt)+5];
                cfg.bpfilttype = 'fir';
                cfg.hilbert = 'abs'; % compute the power
                cfg.padding = 3*max(raw1.time);%
                cfg.padtype = 'mirror';
                cfg.continuous = 'no';
                cfg.channel = raw1.label;
                
                [amp] = ft_preprocessing(cfg,raw1);
                
                % concatenate
                concat2 = zeros(length(amp.label),length(amp.time)*size(amp.trial,1));
                
                idx = 1:length(amp.time);
                for zt = 1:size(amp.trial,1)
                    concat2(:,idx) = squeeze(amp.trial(zt,:,:)).^2;% square the amp to obtain power
                    idx = idx + length(amp.time);
                end;
                clear amp;
                
                %% compute the phase-amplitude histogram
                parfor zt = 1:size(concat1,1)%% loop over MEG channels
                    
                    X = zeros(length(pbins),1);
                    for kt = 1:length(pbins)-1% loop over phase bins
                        
                        [p] = concat1(zt,:);
                        [a] = concat2(zt,:);
                        
                        idx = [];
                        [idx] = find(p  >= pbins(kt) & p < pbins(kt+1));
                        
                        X(kt) = mean(a(idx));
                        
                    end;
                    X(end) = X(1);
                    PAH(zt,it,jt,:) = X;
                end;                
                
            end;
            clear concat*;
        end;        
        
        %% compute MI
        PAH = PAH./repmat(sum(PAH,4),[1 1 1 size(PAH,4)]);
        H = -sum(log(PAH).*PAH,4);
        n = length(pbins);
        MI = (log(n)-H)./log(n);
        
        %% save results
        savepath = p2df;
        save([savepath,'MI_sensorLevel_',files(yt).name],'-v7.3','raw1','PAH','H','n','MI','pbins','pfoi','afoi');
    end;
end;

%%
%exit;