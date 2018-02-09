%%
restoredefaultpath;
addpath('/bcbl/home/home_a-f/froux/fieldtrip-20151020/');
addpath(genpath('~froux/froux/project_reaction_times/mcode/'));
ft_defaults;
savepath = '~froux/froux/project_reaction_times/matFiles/';
%%
SIunit = 'mm';
spm_template = '/bcbl/home/home_a-f/froux/fieldtrip-20151020/external/spm8/templates/T1.nii';
%%
cfg = [];
cfg.dataset = '~froux/froux/project_reaction_times/THA07_DEVELOPMENT_20090924_12.ds';
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
p2d = '~froux/froux/project_reaction_times/matFiles/';
load([p2d,'individual_MNIwarped_grid.mat']);
%%
p2d = '~froux/froux/project_reaction_times/matFiles/';
load([p2d,'individual_headmodel.mat']);
%% compute the coordinates in CTF sensor space for a dipole in the
% left and right occipital cortex
%MNIcoord =[-30 -66 49];% left parietal cortex MNI coordinates
%MNIcoord =[30 -66 49];% right parietal cortex MNI coordinates
%MNIcoord =[-30 -66 49;31 -66 49];% left and right parietal cortex MNI coordinates
%MNIcoord = [-9 -26 9];% left thalamus coordinates
%MNIcoord = [-9 -26 9;9 -26 9];% left and right thalami coordinates

%MNIcoord = [-30 -66 49;-9 -26 9];% left thalamus and left parietal cortex coordinates

%MNIcoord = [10 -24 12;-14 -24 12];%[-30 -66 49;31 -66 49];% left and right parietal cortex MNI coordinates
%MNIcoord = [-9 -26 9;12 -26 9];%[-30 -66 49;31 -66 49];% left and right parietal cortex MNI coordinates
MNIcoord =[-30 -66 49;31 -66 49;-9 -26 9;12 -26 9];% left and right parietal cortex MNI coordinates

[CTFcoord,Gridcoord,VCidx] = convert_MNI2Source_gridCoord(tm,MNIcoord,grid);

idx = find(grid.inside ==1);
if grid.pos(idx(VCidx),:)-Gridcoord ~= 0
    error('VC index does not lign up with output grid coordinates');
end;

[CTFcoord;Gridcoord]
%%
param1 = -pi/2;%-pi:pi/4:pi;%phase lag param
param2 = 1;%0.2:.3:2.1;% SNR param
raw1 = cell(length(param1),length(param2));
raw2 = cell(length(param1),length(param2));

for xt = 1:length(param1)
    for yt = 1:length(param2)
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
        elseif size(Gridcoord,1)==4
            cfg.dip.mom = [-1 0 0 0 0 0 0 0 0 0 0 0]'+[0 0 0 -1 0 0 0 0 0 0 0 0]'+[0 0 0 0 0 0 -1 0 0 0 0 0]'+[0 0 0 0 0 0 0 0 0 -1 0 0]';
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
        
        %snr1 = zeros(cfg.ntrials,1);
        for it = 1:cfg.ntrials
            
            [asig] = A*sin(2*pi*10.*t);
            [gsig] = (sin(2*pi*70.*t)).*(B*sin(2*pi*10.*t)+1);
            
            [sig] = asig+gsig;
            
            %[SNR] = compute_SNR(gsig,asig);
            
            f = 0;
            while f <1
                
                n1 = mean(sig)+(std(sig).*randn(1,cfg.triallength*cfg.fsample));
                n2 = mean(sig)+(std(sig).*randn(1,cfg.triallength*cfg.fsample));
                n3 = mean(sig)+(std(sig).*randn(1,cfg.triallength*cfg.fsample));
                n4 = mean(sig)+(std(sig).*randn(1,cfg.triallength*cfg.fsample));
                
                %[SNR1] = compute_SNR(sig,n1);
                %[SNR2] = compute_SNR(sig,n2);
                
                [sig1] = sig + n1;
                [sig2] = sig + n2;
                
                %[sig1] = sig1.*param2(yt);
                
                [sig3] = hilbert(sig);
                [sig3] = sig3.*exp(1i*param1(xt));
                [sig3] = real(sig3)+n3;
                
                [sig3] = sig3.*param2(yt);
                
                [sig4] = hilbert(sig);
                [sig4] = sig4.*exp(1i*param1(xt));
                [sig4] = real(sig4)+n4;
                
                [sig4] = sig4.*param2(yt);
                [sig4] = sig4.*0.5;
                
                
                r = [];
                r(1) = corr(sig1',sig2');
                r(2) = corr(sig3',sig4');
                r(3) = corr(sig1',sig3');
                r(4) = corr(sig1',sig4');
                r(5) = corr(sig2',sig3');
                r(6) = corr(sig2',sig4');
                
                if sum(r<.6) == length(r)
                    f = 1;
                    %snr1(it) = 10*log10(sqrt(mean(abs(sig2).^2))./sqrt(mean(abs(sig1).^2)));
                end;
            end;
            
            if size(cfg.dip.pos,1)==1
                cfg.dip.signal{it} = sig1;
            elseif size(cfg.dip.pos,1) ==2
                cfg.dip.signal{it} = [sig1;sig2];
            elseif size(cfg.dip.pos,1) ==4
                cfg.dip.signal{it} = [sig1;sig2;sig3;sig4];
            end;
            
        end;
        
        [raw1{xt,yt}] = ft_dipolesimulation(cfg);
        
        cfg = [];
        cfg.channel = meg_data.label;
        
        [raw1{xt,yt}] = ft_selectdata(cfg,raw1{xt,yt});
        
        %SNR1 = zeros(length(raw1{xt,yt}.trial),length(raw1{xt,yt}.label));
        for it = 1:length(raw1{xt,yt}.trial)
            
            m = mean(raw1{xt,yt}.trial{it},2)*ones(1,length(raw1{xt,yt}.time{it}));
            sd = std(raw1{xt,yt}.trial{it},0,2)*ones(1,length(raw1{xt,yt}.time{it}));
            chan_noise = 2*m+2*sd.*randn(length(raw1{xt,yt}.label),length(raw1{xt,yt}.time{1}));
            
            %[SNR1(it,:)] = compute_SNR(raw1{xt,yt}.trial{it},chan_noise);
            
            raw1{xt,yt}.trial{it} = raw1{xt,yt}.trial{it} + chan_noise;
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
        elseif size(Gridcoord,1)==4
            cfg.dip.mom = [-1 0 0 0 0 0 0 0 0 0 0 0]'+[0 0 0 -1 0 0 0 0 0 0 0 0]'+[0 0 0 0 0 0 -1 0 0 0 0 0]'+[0 0 0 0 0 0 0 0 0 -1 0 0]';
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
        
        %snr2 = zeros(cfg.ntrials,1);
        for it = 1:cfg.ntrials
            
            f = 0;
            while f <1
                
                
                [sig1] = A*(.5+(1.*randn(1,cfg.triallength*cfg.fsample)));
                [sig2] = A*(.5+(1.*randn(1,cfg.triallength*cfg.fsample)));
                %[sig1] = sig1.*param2(yt);

                [sig3] = A*(.5+(1.*randn(1,cfg.triallength*cfg.fsample)));                
                [sig3] = hilbert(sig3);
                [sig3] = sig3.*exp(1i*param1(xt));
                [sig3] = real(sig3);
                [sig3] = sig3.*param2(yt);
                
                [sig4] = A*(.5+(1.*randn(1,cfg.triallength*cfg.fsample)));
                [sig4] = hilbert(sig4);
                [sig4] = sig4.*exp(1i*param1(xt));
                [sig4] = real(sig4);
                [sig4] = sig4.*param2(yt);
                [sig4] = sig4.*0.5;
                
                r =[];
                r(1) = corr(sig1',sig2');
                r(2) = corr(sig3',sig4');
                r(3) = corr(sig1',sig3');
                r(4) = corr(sig1',sig4');
                r(5) = corr(sig2',sig3');
                r(6) = corr(sig2',sig4');
                
                if sum(r<.6) == length(r)
                    f = 1;
                    %snr2(it) = 10*log10(sqrt(mean(abs(sig2).^2))./sqrt(mean(abs(sig1).^2)));
                    
                end;
            end;
            
            if size(cfg.dip.pos,1)==1
                cfg.dip.signal{it} = sig1;
            elseif size(cfg.dip.pos,1) ==2
                cfg.dip.signal{it} = [sig1;sig2];
            elseif size(cfg.dip.pos,1) ==4
                cfg.dip.signal{it} = [sig1;sig2;sig3;sig4];
            end;
            
        end;
        
        [raw2{xt,yt}] = ft_dipolesimulation(cfg);
        
        cfg = [];
        cfg.channel = meg_data.label;
        
        [raw2{xt,yt}] = ft_selectdata(cfg,raw2{xt,yt});
        
        %SNR2 = zeros(length(raw2{xt,yt}.trial),length(raw2{xt,yt}.label));
        for it = 1:length(raw2{xt,yt}.trial)
            
            m = mean(raw2{xt,yt}.trial{it},2)*ones(1,length(raw2{xt,yt}.time{it}));
            sd = std(raw2{xt,yt}.trial{it},0,2)*ones(1,length(raw2{xt,yt}.time{it}));
            chan_noise = 2*m+2*sd.*randn(length(raw2{xt,yt}.label),length(raw2{xt,yt}.time{1}));
            
            %[SNR2(it,:)] = compute_SNR(raw2{xt,yt}.trial{it},chan_noise);
            
            raw2{xt,yt}.trial{it} = raw2{xt,yt}.trial{it} + chan_noise;
        end;
    end;
end;
%%
clear sig* t;
%save([savepath,savename5],'raw1','CTFcoord','Gridcoord','MNIcoord','VCidx','SNR1','snr1','-v7.3');
%save([savepath,savename6],'raw2','CTFcoord','Gridcoord','MNIcoord','VCidx','SNR2','snr2','-v7.3');
for it = 1:length(param2)
    
    save_data1 = raw1(:,it);
    save_data2 = raw2(:,it);

    savename5 = ['2x_LR_thalamic_and_2x_LR_parietal_dipoles_alphagammaPAC_SNR',num2str(param2(it)),'.mat'];
    savename6 = ['2x_LR_thalamic_and_2x_LR_parietal_dipoles_nonrhythmic_SNR',num2str(param2(it)),'.mat'];
    
    save([savepath,savename5],'save_data1','CTFcoord','Gridcoord','MNIcoord','VCidx','-v7.3');
    save([savepath,savename6],'save_data2','CTFcoord','Gridcoord','MNIcoord','VCidx','-v7.3');
end;
%%
clear all;
exit;

