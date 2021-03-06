%% set path def
restoredefaultpath;
addpath(genpath('~rouxf/projects/TC/mcode/'));
addpath('~rouxf/toolboxes/fieldtrip-20161009/');
ft_defaults;

[savepath] = '~rouxf/projects/TC/matFiles/';

%% set the unit for geometrical data
SIunit = 'mm';
spm_template = '~rouxf/toolboxes/fieldtrip-20161009/external/spm8/templates/T1.nii';

%% load single subject MEG data
load([savepath,'resting_state_recording_THA07_ICAcleanedMEGdataEO.mat']);
load([savepath,'resting_state_recording_THA07_ICAcleanedMEGdataEC.mat']);

%%
ntrl = length(dc1.trial);% number of trials
trlt = max(dc1.time{1});% duration of each trial in s
Fs = 1/(dc1.time{1}(2)-dc1.time{1}(1));% sampling freq

cfg = [];
cfg.layout = 'CTF275.lay';

[lay] = ft_prepare_layout(cfg);

meg_data1 = dc1; clear dc1;
meg_data2 = dc2; clear dc2;

meg_data1.grad = ft_convert_units(meg_data1.grad,SIunit);
meg_data2.grad = ft_convert_units(meg_data2.grad,SIunit);

meg_data1.fsample = Fs;
meg_data2.fsample = Fs;

%% load invidual grid
p2d = savepath;
load([p2d,'individual_MNIwarped_grid.mat']);

%%
p2d = savepath;
load([p2d,'individual_headmodel.mat']);

%% compute the coordinates in CTF sensor space for a dipole in the
% left and right occipital cortex
%MNIcoord =[-30 -66 49];% left parietal cortex MNI coordinates
%MNIcoord =[30 -66 49];% right parietal cortex MNI coordinates
%MNIcoord =[-30 -66 49;31 -66 49];% left and right parietal cortex MNI coordinates
%MNIcoord = [-9 -26 9];% left thalamus coordinates
%MNIcoord = [-9 -26 9;9 -26 9];% left and right thalami coordinates

%MNIcoord = [-30 -66 49;-9 -26 9];% left thalamus and left parietal cortex coordinates


%MNIcoord =[-30 -66 49;31 -66 49];% left and right parietal cortex MNI coordinates

MNIcoord =[-30 -66 49;... %left parietal cortex MNI coordinates
    -44 8 43;...%left FEF MNI coordinates
    -10 -18 10];  % left thalamus

% MNIcoord =[-30 -66 49;... % left and right parietal cortex MNI coordinates
%            31 -66 49;...
%           -14 -24 12;... % left and right thalamus MNI coordinates
%            10 -24 12];

[CTFcoord,Gridcoord,VCidx] = convert_MNI2Source_gridCoord(tm,MNIcoord,grid);

idx = find(grid.inside ==1);
if grid.pos(idx(VCidx),:)-Gridcoord ~= 0
    error('VC index does not lign up with output grid coordinates');
end;

[CTFcoord;Gridcoord]

%%
% figure;
% hold on;
% ft_plot_vol(hdm,'facealpha',.7,'surfaceonly',0);
% ft_plot_mesh(grid.pos(grid.inside,:),'facecolor','cortex','surfaceonly',0);
% 
% plot3(Gridcoord(1,1),Gridcoord(1,2),Gridcoord(1,3),'ro','MarkerFaceColor','r');
% plot3(Gridcoord(2,1),Gridcoord(2,2),Gridcoord(2,3),'ro','MarkerFaceColor','r');
% plot3(Gridcoord(3,1),Gridcoord(3,2),Gridcoord(3,3),'ro','MarkerFaceColor','r');

%%
param1 = -pi/2;%-pi:pi/4:pi;%phase lag param
param2 = 1;%0.2:.3:2.1;% SNR param
param3 = 10; % number of random sources

raw1 = cell(length(param1),length(param2));
raw2 = cell(length(param1),length(param2));

%%
rand_mom = [];
rand_pos = [];

idx(VCidx) = [];

mom = [-1 0 1];
for zt = 1:param3
    idx = idx( randperm(length(idx)) );
    ix = idx(1);
    idx(1) = [];
    rand_pos = [rand_pos;grid.pos(ix,:)];
    fl = zeros(1,3);
    for yt = 1:3
        x = randperm(3);
        fl(yt) = x(1);
    end;
    rand_mom = [rand_mom;mom(fl)'];
end;

%%
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
        cfg.grad = meg_data2.grad;% gradiometer positions
        
        %Position of dipole in cartesian space (x,y,z)
        cfg.dip.pos = [Gridcoord;rand_pos];
        if size(Gridcoord,1) ==1
            cfg.dip.mom = [-1 0 0]';
        elseif size(Gridcoord,1)==2
            cfg.dip.mom = [-1 0 0 0 0 0 zeros(1,param3*3)]'+[0 0 0 -1 0 0  zeros(1,param3*3)]';
        elseif size(Gridcoord,1)==3
            cfg.dip.mom = [-1 0 0 0 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 -1 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 0 0 0 -1 0 0  zeros(1,param3*3)]';
        elseif size(Gridcoord,1)==4
            cfg.dip.mom = [-1 0 0 0 0 0 0 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 -1 0 0 0 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 0 0 0 -1 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 0 0 0 0 0 0 -1 0 0  zeros(1,param3*3)]';
        end;
        cfg.dip.mom(size(Gridcoord,1)*3+1:end) = rand_mom;
        
        cfg.ntrials = ntrl;
        cfg.triallength = trlt;
        cfg.fsample = Fs;%meg_data.fsample;
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
                
                n1 = mean(sig)+1*(std(sig).*randn(1,cfg.triallength*cfg.fsample));
                n2 = mean(sig)+1*(std(sig).*randn(1,cfg.triallength*cfg.fsample));
                n3 = mean(sig)+1*(std(sig).*randn(1,cfg.triallength*cfg.fsample));
                n4 = mean(sig)+1*(std(sig).*randn(1,cfg.triallength*cfg.fsample));
                
                % 1st cortical dipole
                %[sig1] = sig + n1;
                [sig1] = gsig + n1;
                
                %apply phase shift to 2nd cortical dipole
                %[sig3] = hilbert(gsig);
                %[sig3] = sig3.*exp(1i*param1(xt));
                %[sig3] = asig + real(sig3) + n3;
                
                [sig3] = hilbert(gsig);
                [sig3] = sig3.*exp(1i*param1(xt));
                [sig3] = real(sig3) + n3;
                
                %apply phase shift and signal boost to 1st thalamic dipole
                %                 [sig2] = hilbert(asig);
                %                 [sig2] = sig2.*exp(1i*param1(xt));
                %                 [sig2] = real(sig2) + n2;
                %                 [sig2] = sig2.*param2(yt);
                
                [sig2] = hilbert(asig);
                [sig2] = sig2.*exp(1i*param1(xt));
                [sig2] = real(sig2) + n2;
                [sig2] = sig2.*param2(yt);
                
                %apply phase shift and signal boost to 2nd thalamic dipole
                %                 [sig4] = hilbert(asig);
                %                 [sig4] = sig4.*exp(1i*param1(xt));
                %                 [sig4] = real(sig4) + n4;
                %                 [sig4] = sig4.*param2(yt);
                
                [sig4] = hilbert(asig);
                [sig4] = sig4.*exp(1i*param1(xt));
                [sig4] = real(sig4) + n4;
                [sig4] = sig4.*param2(yt);
                
                r = [corr(sig1',sig2') corr(sig1',sig3') corr(sig1',sig4') corr(sig2',sig3') corr(sig2',sig4') corr(sig3',sig4')];
                
                if max(abs(r)) <.6
                    f = 1;
                    %snr1(it) = 10*log10(sqrt(mean(abs(sig2).^2))./sqrt(mean(abs(sig1).^2)));
                end;
            end;
            
            if size(cfg.dip.pos,1)-param3==1
                cfg.dip.signal{it} = sig1;% 1x cortex
            elseif size(cfg.dip.pos,1)-param3 ==2
                cfg.dip.signal{it} = [sig1;sig2];% 1 x cortex + 1 x thalamus
            elseif size(cfg.dip.pos,1)-param3 ==3
                cfg.dip.signal{it} = [sig1;sig3;sig2];% 2 x cortex + 1 x thalamus
            elseif size(cfg.dip.pos,1)-param3 ==4
                cfg.dip.signal{it} = [sig1;sig3;sig2;sig4];% 2 x cortex + 2 x thalamus
            end;
            
            for jt = 1:size(rand_pos,1)
                rand_sig = mean(sig)+1*(std(sig).*randn(1,cfg.triallength*cfg.fsample));
                cfg.dip.signal{it} = [cfg.dip.signal{it};rand_sig];
            end;
            
        end;
        
        [raw1{xt,yt}] = ft_dipolesimulation(cfg);
        
        cfg = [];
        cfg.channel = {'MEG'};
        
        [raw1{xt,yt}] = ft_selectdata(cfg,raw1{xt,yt});
        
        cfg = [];
        cfg.channel = setdiff(meg_data2.label,setdiff(meg_data2.label,raw1{xt,yt}.label))
        
        [meg_data2] = ft_selectdata(cfg,meg_data2);
        
        %SNR1 = zeros(length(raw1{xt,yt}.trial),length(raw1{xt,yt}.label));
        for it = 1:length(raw1{xt,yt}.trial)
            
            m = mean(raw1{xt,yt}.trial{it},2)*ones(1,length(raw1{xt,yt}.time{it}));
            sd = std(raw1{xt,yt}.trial{it},0,2)*ones(1,length(raw1{xt,yt}.time{it}));
            chan_noise = 2*m+2*sd.*randn(length(raw1{xt,yt}.label),length(raw1{xt,yt}.time{1}));
            
            %[SNR1(it,:)] = compute_SNR(raw1{xt,yt}.trial{it},chan_noise);
            
            raw1{xt,yt}.trial{it} = raw1{xt,yt}.trial{it} + chan_noise + meg_data2.trial{it};
        end;
        %%
        cfg = [];
        cfg.headmodel = hdm;%forward model
        cfg.grad = meg_data2.grad;% gradiometer positions
        
        %Position of dipole in cartesian space (x,y,z)
        cfg.dip.pos = [Gridcoord;rand_pos];
        if size(Gridcoord,1) ==1
            cfg.dip.mom = [-1 0 0]';
        elseif size(Gridcoord,1)==2
            cfg.dip.mom = [-1 0 0 0 0 0 zeros(1,param3*3)]'+[0 0 0 -1 0 0  zeros(1,param3*3)]';
        elseif size(Gridcoord,1)==3
            cfg.dip.mom = [-1 0 0 0 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 -1 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 0 0 0 -1 0 0  zeros(1,param3*3)]';
        elseif size(Gridcoord,1)==4
            cfg.dip.mom = [-1 0 0 0 0 0 0 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 -1 0 0 0 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 0 0 0 -1 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 0 0 0 0 0 0 -1 0 0  zeros(1,param3*3)]';
        end;
        cfg.dip.mom(size(Gridcoord,1)*3+1:end) = rand_mom;
        
        % cfg.dip.frequency = [10 10];
        % cfg.dip.phase = [-pi pi];
        
        cfg.ntrials = ntrl;
        cfg.triallength = trlt;
        cfg.fsample = Fs;%meg_data.fsample;
        %cfg.relnoise    = 1.25;
        cfg.relnoise    = 1.75;
        
        cfg.dip.signal = cell(1,cfg.ntrials);
        
        A = 2.5;
        
        %snr2 = zeros(cfg.ntrials,1);
        for it = 1:cfg.ntrials
            
            f = 0;
            while f <1
                
                [sig1] = A*(.5+(1.*randn(1,cfg.triallength*cfg.fsample)));
                [sig3] = A*(.5+(1.*randn(1,cfg.triallength*cfg.fsample)));
                
                %[sig1] = sig1.*param2(yt);
                
                [sig2] = A*(.5+(1.*randn(1,cfg.triallength*cfg.fsample)));
                [sig2] = hilbert(sig2);
                [sig2] = sig2.*exp(1i*param1(xt));
                [sig2] = real(sig2);
                [sig2] = sig2.*param2(yt);
                
                [sig4] = A*(.5+(1.*randn(1,cfg.triallength*cfg.fsample)));
                [sig4] = hilbert(sig4);
                [sig4] = sig4.*exp(1i*param1(xt));
                [sig4] = real(sig4);
                [sig4] = sig4.*param2(yt);
                
                r = [corr(sig1',sig2') corr(sig1',sig3') corr(sig1',sig4') corr(sig2',sig3') corr(sig2',sig4') corr(sig3',sig4')];
                
                if max(abs(r)) <.6
                    f = 1;
                    %snr2(it) = 10*log10(sqrt(mean(abs(sig2).^2))./sqrt(mean(abs(sig1).^2)));
                    
                end;
            end;
            
            if size(cfg.dip.pos,1)-param3==1
                cfg.dip.signal{it} = sig1;% 1x cortex
            elseif size(cfg.dip.pos,1)-param3 ==2
                cfg.dip.signal{it} = [sig1;sig2];% 1 x cortex + 1 x thalamus
            elseif size(cfg.dip.pos,1)-param3 ==3
                cfg.dip.signal{it} = [sig1;sig3;sig2];% 2 x cortex + 1 x thalamus
            elseif size(cfg.dip.pos,1)-param3 ==4
                cfg.dip.signal{it} = [sig1;sig3;sig2;sig4];% 2 x cortex + 2 x thalamus
            end;
            
            for jt = 1:size(rand_pos,1)
                rand_sig = mean(sig)+1*(std(sig).*randn(1,cfg.triallength*cfg.fsample));
                cfg.dip.signal{it} = [cfg.dip.signal{it};rand_sig];
            end;
            
        end;
        
        [raw2{xt,yt}] = ft_dipolesimulation(cfg);
        
        cfg = [];
        cfg.channel = meg_data2.label;
        
        [raw2{xt,yt}] = ft_selectdata(cfg,raw2{xt,yt});
        
        %SNR2 = zeros(length(raw2{xt,yt}.trial),length(raw2{xt,yt}.label));
        for it = 1:length(raw2{xt,yt}.trial)
            
            m = mean(raw2{xt,yt}.trial{it},2)*ones(1,length(raw2{xt,yt}.time{it}));
            sd = std(raw2{xt,yt}.trial{it},0,2)*ones(1,length(raw2{xt,yt}.time{it}));
            chan_noise = 2*m+2*sd.*randn(length(raw2{xt,yt}.label),length(raw2{xt,yt}.time{1}));
            
            %[SNR2(it,:)] = compute_SNR(raw2{xt,yt}.trial{it},chan_noise);
            
            raw2{xt,yt}.trial{it} = raw2{xt,yt}.trial{it} + chan_noise + meg_data2.trial{it};
        end;
    end;
end;

%%
clear sig* t;
for it = 1:length(param2)
    save_data1 = raw1(:,it);
    save_data2 = raw2(:,it);
    
    savename5 = ['1x_L_thalamic_and_2x_L_cortical_dipoles_alphagammaPAC_SNR',num2str(param2(it)),'_nRandSources',num2str(param3),'_realMEG.mat'];
    savename6 = ['1x_L_thalamic_and_2x_L_cortical_dipoles_nonrhythmic_SNR',num2str(param2(it)),'_nRandSources',num2str(param3),'_realMEG.mat'];
    
    save([savepath,savename5],'save_data1','CTFcoord','Gridcoord','MNIcoord','VCidx','-v7.3');
    save([savepath,savename6],'save_data2','CTFcoord','Gridcoord','MNIcoord','VCidx','-v7.3');
end;

%%
clear all;
exit;