%%
restoredefaultpath;
addpath('/bcbl/home/home_a-f/froux/fieldtrip-20151020/');
addpath(genpath('~froux/froux/dipole_simulation_thalamicPAC/'));
ft_defaults;
savepath = '~froux/froux/dipole_simulation_thalamicPAC/matFiles/';
%%
SIunit = 'mm';
spm_template = '/bcbl/home/home_a-f/froux/fieldtrip-20151020/external/spm8/templates/T1.nii';
%%
cfg = [];
cfg.dataset = '/bcbl/home/home_a-f/froux/dipole_simulation_thalamicPAC/THA07_DEVELOPMENT_20090924_12.ds';
cfg.continous = 'yes';
cfg.trl = [240*1200+1 480*1200 0];
cfg.demean = 'yes';
cfg.detrend = 'yes';
cfg.channel = {'MEG' '-MLP12' '-MLT41' '-MRC12'  '-MRC14' '-MRC25' '-MRP56' '-MRT12' '-MRT21' '-MRT23' '-MLT57' '-MLT52' '-MRF22' '-MRF13' '-MRF24' '-MRF43' '-MLO22' '-MLF25' '-MRO31' '-MRO21'};

[meg_data] = ft_preprocessing(cfg);

meg_data.grad = ft_convert_units(meg_data.grad,SIunit);
%%
cfg = [];
cfg.channel = {'MEG'};
cfg.length = 7;

[meg_data] = ft_redefinetrial(cfg,meg_data);
%%
p2d = '~froux/froux/dipole_simulation_thalamicPAC/matFiles/';
load([p2d,'individual_MNIwarped_grid.mat']);
%%
p2d = '~froux/froux/dipole_simulation_thalamicPAC/matFiles/';
load([p2d,'individual_headmodel.mat']);
%% compute the coordinates in CTF sensor space for a dipole in the
% left and right occipital cortex
%MNIcoord =[-30 -66 49];% left parietal cortex MNI coordinates
%MNIcoord =[30 -66 49];% right parietal cortex MNI coordinates
%MNIcoord =[-30 -66 49;31 -66 49];% left and right parietal cortex MNI coordinates
%MNIcoord = [-9 -26 9];% left thalamus coordinates
%MNIcoord = [-9 -26 9;9 -26 9];% left and right thalami coordinates

MNIcoord = [-30 -66 49;-9 -26 9];% left thalamus and left parietla cortex coordinates

%MNIcoord =[-30 -66 49;31 -66 49];% left and right parietal cortex MNI coordinates

[CTFcoord,Gridcoord,VCidx] = convert_MNI2Source_gridCoord(tm,MNIcoord,grid);

idx = find(grid.inside ==1);
if grid.pos(idx(VCidx),:)-Gridcoord ~= 0
    error('VC index does not lign up with output grid coordinates');
end;

[CTFcoord;Gridcoord]
%%
param= 1:.5:10;
for xt = 1:length(param)
    %Orientation of the dipole [qx qy qz] in MNI corrdinates: ARS
    %note: x-axis goes towards R, y-axis goes towards anterior, z-axis
    %goes towards superior
    %[1 0 0] points to Nasion
    %[0 1 0] points to LPA
    %[0 0 1] points to vertex
    %[-1 0 0] points to Inion
    %[0 -1 0] points to RPA
    %[0 0 -1] points away from vertex towards the feet
    
    cfg = [];
    cfg.headmodel = hdm;%forward model
    cfg.grad = meg_data.grad;% gradiometer positions
    
    %Position of dipole in cartesian space (x,y,z)
    cfg.dip.pos = [Gridcoord];
    if size(Gridcoord,1) ==1
        cfg.dip.mom = [-1 0 0]';
    elseif size(Gridcoord,1)==2
        cfg.dip.mom = [-1 0 0 0 0 0]'+[0 0 0 -1 0 0]';
    end;
    
    cfg.ntrials = length(meg_data.trial);
    cfg.triallength = round(abs(min(meg_data.time{1})) + max(meg_data.time{1}));
    cfg.fsample = meg_data.fsample;%meg_data.fsample;
    cfg.relnoise    = 0;
    
    cfg.dip.signal = cell(1,cfg.ntrials);
    t = 0:1/cfg.fsample:cfg.triallength;
    t(1) = [];
    A = 2.5;
    B = 1.5;
    
    snr1 = zeros(cfg.ntrials,1);
    for it = 1:cfg.ntrials
        
        [asig] = A*sin(2*pi*10.*t);
        [gsig] = (sin(2*pi*70.*t)).*(B*sin(2*pi*10.*t)+1);
        
        [sig] = asig+gsig;
        
        [SNR] = compute_SNR(gsig,asig);
        
        f = 0;
        while f <1
            
            n1 = mean(sig)+(std(sig).*randn(1,cfg.triallength*cfg.fsample));
            n2 = mean(sig)+(std(sig).*randn(1,cfg.triallength*cfg.fsample));
            
            [SNR1] = compute_SNR(sig,n1);
            [SNR2] = compute_SNR(sig,n2);
            
            [sig1] = sig + n1;
            
            [sig2] = hilbert(sig);
            [sig2] = sig2.*exp(1i*pi/4);
            [sig2] = real(sig2)+n2;
            
            [sig2] = sig2.*param(xt);
            
            
            r = corr(sig1',sig2');
            
            if r <.6
                f = 1;
                snr1(it) = 10*log10(sqrt(mean(abs(sig2).^2))./sqrt(mean(abs(sig1).^2)));
            end;
        end;
        
        if size(cfg.dip.pos,1)==1
            cfg.dip.signal{it} = sig1;
        elseif size(cfg.dip.pos,1) ==2
            cfg.dip.signal{it} = [sig1;sig2];
        end;
        
    end;
    
    [raw1] = ft_dipolesimulation(cfg);
    
    cfg = [];
    cfg.channel = meg_data.label;
    
    [raw1] = ft_selectdata(cfg,raw1);
    
    SNR1 = zeros(length(raw1.trial),length(raw1.label));
    for it = 1:length(raw1.trial)
        
        m = mean(raw1.trial{it},2)*ones(1,length(raw1.time{it}));
        sd = std(raw1.trial{it},0,2)*ones(1,length(raw1.time{it}));
        chan_noise = 2*m+2*sd.*randn(length(raw1.label),length(raw1.time{1}));
        
        [SNR1(it,:)] = compute_SNR(raw1.trial{it},chan_noise);
        
        raw1.trial{it} = raw1.trial{it} + chan_noise;
    end;
    %%
    cfg = [];
    cfg.headmodel = hdm;%forward model
    cfg.grad = meg_data.grad;% gradiometer positions
    
    %Position of dipole in cartesian space (x,y,z)
    cfg.dip.pos = [Gridcoord];
    if size(Gridcoord,1) ==1
        cfg.dip.mom = [-1 0 0]';
    elseif size(Gridcoord,1)==2
        cfg.dip.mom = [-1 0 0 0 0 0]'+[0 0 0 -1 0 0]';
    end;
    
    % cfg.dip.frequency = [10 10];
    % cfg.dip.phase = [-pi pi];
    
    cfg.ntrials = length(meg_data.trial);
    cfg.triallength = round(abs(min(meg_data.time{1})) + max(meg_data.time{1}));
    cfg.fsample = meg_data.fsample;%meg_data.fsample;
    %cfg.relnoise    = 1.25;
    cfg.relnoise    = 1.75;
    
    cfg.dip.signal = cell(1,cfg.ntrials);
    
    A = 2.5;
    
    snr2 = zeros(cfg.ntrials,1);
    for it = 1:cfg.ntrials
        
        f = 0;
        while f <1
            
            [sig1] = A*(.5+(1.*randn(1,cfg.triallength*cfg.fsample)));
            [sig2] = A*(.5+(1.*randn(1,cfg.triallength*cfg.fsample)));
            [sig2] = sig2.*param(xt);
            
            [sig2] = hilbert(sig2);
            [sig2] = sig2.*exp(1i*pi/4);
            [sig2] = real(sig2);
            
            
            r = corr(sig1',sig2');
            
            if r <.6
                f = 1;
                snr2(it) = 10*log10(sqrt(mean(abs(sig2).^2))./sqrt(mean(abs(sig1).^2)));
                
            end;
        end;
        
        if size(cfg.dip.pos,1)==1
            cfg.dip.signal{it} = sig1;
        elseif size(cfg.dip.pos,1) ==2
            cfg.dip.signal{it} = [sig1;sig2];
        end;
        
    end;
    
    [raw2] = ft_dipolesimulation(cfg);
    
    cfg = [];
    cfg.channel = meg_data.label;
    
    [raw2] = ft_selectdata(cfg,raw2);
    
    SNR2 = zeros(length(raw2.trial),length(raw2.label));
    for it = 1:length(raw2.trial)
        
        m = mean(raw2.trial{it},2)*ones(1,length(raw2.time{it}));
        sd = std(raw2.trial{it},0,2)*ones(1,length(raw2.time{it}));
        chan_noise = 2*m+2*sd.*randn(length(raw2.label),length(raw2.time{1}));
        
        [SNR2(it,:)] = compute_SNR(raw2.trial{it},chan_noise);
        
        raw2.trial{it} = raw2.trial{it} + chan_noise;
    end;
end;
%%
clear sig* t;
savename5 = 'left_thalamic_and_left_parietal_dipoles_alphagammaPAC.mat';
savename6 = 'left_thalamic_and_left_parietal_dipoles_nonrhythmic.mat';
save([savepath,savename5],'raw1','CTFcoord','Gridcoord','MNIcoord','VCidx','SNR1','snr1');
save([savepath,savename6],'raw2','CTFcoord','Gridcoord','MNIcoord','VCidx','SNR2','snr2');
exit;

