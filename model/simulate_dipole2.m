function simulate_dipole2(param1,param2,param3)

if nargin ==0
    param1 = -pi/2;%-pi:pi/4:pi;%phase lag param
    param2 = 1;%0.2:.3:2.1;% SNR param
    param3 = 0;%;10e1; % number of random sources
end;

%%
[savepath] = '~rouxf/prj/TC/matFiles/';

%% set the unit for geometrical data
SIunit = 'mm';

%% load single subject MEG data
load([savepath,'resting_state_recording_THA07_ICAcleanedMEGdataEO.mat']);
load([savepath,'resting_state_recording_THA07_ICAcleanedMEGdataEC.mat']);

load([savepath,'independent_alphaGenerators_timeSeries.mat']);
SigDat = dat;
for it = 1:size(r,1)
    r(it,it) = NaN;
end;

%%
ntrl = params.ntrl;%length(dc1.trial);% number of trials
trlt = max(params.t);%max(dc1.time{1});% duration of each trial in s
Fs = params.Fs;% sampling freq

%%
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
dat = load([p2d,'individual_MNIwarped_grid.mat']);
Sgrid = dat.grid;
tm = dat.tm;

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

%MNIcoord = [-30 -66 49;-10 -18 10];% left thalamus and left parietal cortex coordinates


%MNIcoord =[-30 -66 49;31 -66 49];% left and right parietal cortex MNI coordinates

%MNIcoord =[-30 -66 49;... %left parietal cortex MNI coordinates
%    -44 8 43;...%left FEF MNI coordinates
%    -10 -18 10];  % left thalamus

% MNIcoord =[-30 -66 49;... % left and right parietal cortex MNI coordinates
%             31 -66 49;...
%            -14 -24 12;... % left and right thalamus MNI coordinates
%             10 -24 12];

MNIcoord =[-14 -24 12;... % left and right thalamus MNI coordinates
    10 -24 12;...
    -18 -38 46;... % left and right parietal cortex MNI coordinates
    28 -28 46;...
    -19 -79 20;... % left and right occipital cortex MNI coordinates
    21 -83 19];

[CTFcoord,Gridcoord,VCidx] = convert_MNI2Source_gridCoord(tm,MNIcoord,Sgrid);

[CTFcoord;Gridcoord]

%%
C = [.9 0 0; 0 0 .9; .9 .75 .75; 0 .9 0; 0 0.9 .9; 0 0.5 0.5; 0.9 0.25 0];
figure;
hold on;
ft_plot_vol(hdm,'facealpha',.7,'surfaceonly',0);
%ft_plot_mesh(Sgrid.pos(Sgrid.inside,:),'facecolor','cortex','surfaceonly',0);

for it = 1:size(Gridcoord,1)
    plot3(Gridcoord(it,1),Gridcoord(it,2),Gridcoord(it,3),'o','Color',C(it,:),'MarkerFaceColor',C(it,:),'MarkerSize',6);
end;
view(153,52);

%%
raw1 = cell(length(param1),length(param2));
raw2 = cell(length(param1),length(param2));

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
        cfg.dipoleunit = 'nA*m';
        cfg.channelunit = 'T/m';
        cfg.headmodel = hdm;%forward model
        cfg.grad = meg_data2.grad;% gradiometer positions                        
        cfg.ntrials = ntrl;
        cfg.triallength = trlt;
        cfg.fsample = Fs;%meg_data.fsample;
        cfg.relnoise    = 0;
        
        cfg.dip.signal = cell(1,cfg.ntrials);
        
        %snr1 = zeros(cfg.ntrials,1);
        for it = 1:cfg.ntrials
            
            %%
            rand_mom = [];
            rand_pos = [];
                       
            idx = find(Sgrid.inside ==1);
            if Sgrid.pos(idx(VCidx),:)-Gridcoord ~= 0
                error('VC index does not lign up with output grid coordinates');
            end;
            idx(VCidx) = [];
            
            mom = [-1 0 1];
            for zt = 1:param3
                idx = idx( randperm(length(idx)) );
                ix = idx(1);
                idx(1) = [];
                rand_pos = [rand_pos;Sgrid.pos(ix,:)];
                fl = zeros(1,3);
                for vt = 1:3
                    x = randperm(3);
                    fl(vt) = x(1);
                end;
                rand_mom = [rand_mom;mom(fl)'];
            end;
            
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
            elseif size(Gridcoord,1)==6
                cfg.dip.mom = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 zeros(1,param3*3)]'+[0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0  zeros(1,param3*3)]';
            end;
            cfg.dip.mom(size(Gridcoord,1)*3+1:end) = rand_mom;
            
            %[SNR] = compute_SNR(gsig,asig);
            
            if max(max(abs(r))) <.6
                f = 1;
                %snr1(it) = 10*log10(sqrt(mean(abs(sig2).^2))./sqrt(mean(abs(sig1).^2)));
            else
                error('Signal correlations must remain < 0.6');
            end;
            
            if size(cfg.dip.pos,1)-param3==1
                cfg.dip.signal{it} = squeeze(SigDat(3,it,:));% 1x cortex
            elseif size(cfg.dip.pos,1)-param3 ==2
                cfg.dip.signal{it} = squeeze(SigDat([1 3],it,:));% 1 x thalamus + 1 x cortex
            elseif size(cfg.dip.pos,1)-param3 ==3
                cfg.dip.signal{it} = squeeze(SigDat([1 3 4],it,:));% 1 x thalamus + 2 x cortex
            elseif size(cfg.dip.pos,1)-param3 ==4
                cfg.dip.signal{it} = squeeze(SigDat(1:4,it,:));%  2 x thalamus + 2 x cortex
            elseif size(cfg.dip.pos,1)-param3 ==6
                cfg.dip.signal{it} = squeeze(SigDat(:,it,:));%  2 x thalamus + 6 x cortex
            end;
            
            for jt = 1:size(rand_pos,1)
                rand_sig = mean(mean(cfg.dip.signal{it}(3:4,:)))+1*mean(std(cfg.dip.signal{it}(3:4,:),0,2)).*randn(1,cfg.triallength*cfg.fsample+1);
                rand_sig = rand_sig./10;
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
            
            raw1{xt,yt}.trial{it} = raw1{xt,yt}.trial{it} + chan_noise;%
        end;
        
        %%
        cfg = [];
        cfg.dipoleunit = 'nA*m';
        cfg.channelunit = 'T/m';
        cfg.headmodel = hdm;%forward model
        cfg.grad = meg_data2.grad;% gradiometer positions                        
        cfg.ntrials = ntrl;
        cfg.triallength = trlt;
        cfg.fsample = Fs;%meg_data.fsample;
        %cfg.relnoise    = 0.15;
        cfg.relnoise    = 1.75;
        
        A = 1;
        cfg.dip.signal = cell(1,cfg.ntrials);
        %snr2 = zeros(cfg.ntrials,1);
        for it = 1:cfg.ntrials                        
            
            %%
            rand_mom = [];
            rand_pos = [];
            
            idx = find(Sgrid.inside ==1);
            if Sgrid.pos(idx(VCidx),:)-Gridcoord ~= 0
                error('VC index does not lign up with output grid coordinates');
            end;
            idx(VCidx) = [];
            
            mom = [-1 0 1];
            for zt = 1:param3
                idx = idx( randperm(length(idx)) );
                ix = idx(1);
                idx(1) = [];
                rand_pos = [rand_pos;Sgrid.pos(ix,:)];
                fl = zeros(1,3);
                for vt = 1:3
                    x = randperm(3);
                    fl(vt) = x(1);
                end;
                rand_mom = [rand_mom;mom(fl)'];
            end;
            
            sig = zeros(size(Gridcoord,1),cfg.triallength*cfg.fsample+1);
            for jt = 1:size(Gridcoord,1)
                [sig(jt,:)] = A*(.5+(1.*randn(1,cfg.triallength*cfg.fsample+1)));
            end;
            
            for jt = 2:6
                [sig(jt,:)] = hilbert(sig(jt,:));
                [sig(jt,:)] = sig(jt,:).*exp(1i*param1(xt));
                [sig(jt,:)] = real(sig(jt,:));
                [sig(jt,:)] = sig(jt,:).*param2(yt);
            end;
            sig(3:4,:) = sig(3:4,:)./100;
            sig([1:2 5:6],:) = sig([1:2 5:6],:)./10;
            
            r = zeros(size(sig,1),size(sig,1));
            for jt = 1:size(sig,1)
                for kt = 1:size(sig,1)
                    r(jt,kt) = corr(sig(jt,:)',sig(kt,:)');
                end;
                r(jt,jt) = 0;
            end;
            
            if max(max(abs(r))) <.6
                %snr2(it) = 10*log10(sqrt(mean(abs(sig2).^2))./sqrt(mean(abs(sig1).^2)));
            else
                error('Signal correlations must remain <0.6');
            end;
            
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
            elseif size(Gridcoord,1)==6
                cfg.dip.mom = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 zeros(1,param3*3)]'+[0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0  zeros(1,param3*3)]'+[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0  zeros(1,param3*3)]';
            end;
            cfg.dip.mom(size(Gridcoord,1)*3+1:end) = rand_mom;
        
            if size(cfg.dip.pos,1)-param3==1
                cfg.dip.signal{it} = sig(3,:);% 1x cortex
            elseif size(cfg.dip.pos,1)-param3 ==2
                cfg.dip.signal{it} = sig([1 3],:,:);% 1 x thalamus + 1 x cortex
            elseif size(cfg.dip.pos,1)-param3 ==3
                cfg.dip.signal{it} = sig([1 3 4],:);%1 x thalamus + 2 x cortex
            elseif size(cfg.dip.pos,1)-param3 ==4
                cfg.dip.signal{it} = sig(1:4,:);% 2 x thalamus + 2 x cortex
            elseif size(cfg.dip.pos,1)-param3 ==6
                cfg.dip.signal{it} = sig;% 2 x thalamus + 4 x cortex
            end;
            
            for jt = 1:size(rand_pos,1)
                rand_sig = mean(mean(sig(3:4,:)))+1*(mean(std(sig(3:4,:),0,2)).*randn(1,cfg.triallength*cfg.fsample+1));
                rand_sig = rand_sig./10;
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
            
            raw2{xt,yt}.trial{it} = raw2{xt,yt}.trial{it} + chan_noise;%
        end;
    end;
end;

%%
clear sig* t;
for it = 1:length(param2)
    save_data1 = raw1(:,it);
    save_data2 = raw2(:,it);
    
    savename5 = ['2x_L_thalamic_and_4x_L_cortical_dipoles_alphagammaPAC_SNR',num2str(param2(it)),'_nRandSources',num2str(param3),'_simMEG.mat'];
    savename6 = ['2x_L_thalamic_and_4x_L_cortical_dipoles_nonrhythmic_SNR',num2str(param2(it)),'_nRandSources',num2str(param3),'_simMEG.mat'];
    
    save([savepath,savename5],'save_data1','CTFcoord','Gridcoord','MNIcoord','VCidx','-v7.3');
    save([savepath,savename6],'save_data2','CTFcoord','Gridcoord','MNIcoord','VCidx','-v7.3');
end;
