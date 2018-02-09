%% set path defs
restoredefaultpath;

%fieldtrip
addpath(genpath('~rouxf/tbx/chronux_2_11/'));
ft = dir('~rouxf/tbx/fieldtrip-*');
addpath(['~rouxf/tbx/',ft.name]);
ft_defaults;

%local scripts/functions
addpath(genpath('/home/rouxf/prj/TC/mcode/'));

%%
%simulate_independent_alpha_generators([1 1 1],1);
%simulate_dipole_ROI(0,1,0,[5 5 5]);

%%            simulate_dipole_ROI(snr,lt,nDip(kt),scf(jt,:));

load sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources100_simMEG_ROI_fieldMode1.mat
dc3 = save_data1{1};

load(['/home/rouxf/prj/TC/matFiles/','resting_state_recording_THA07_ICAcleanedMEGdataEO.mat']);
load(['/home/rouxf/prj/TC/matFiles/','resting_state_recording_THA07_ICAcleanedMEGdataEC.mat']);
 
cfg                     = [];
cfg.length              = length(dc1.time{1})/dc1.fsample;

[dc3] = ft_redefinetrial( cfg , dc3);
%[dc2] = ft_redefinetrial( cfg , dc2);

sel = randperm(length(dc3.trial));

cfg                     = [];
cfg.trials              = sel(1:length(dc1.trial));

[dc3] = ft_selectdata( cfg, dc3 );
%[dc2] = ft_selectdata( cfg, dc2 );

%% path for read and write 
[p2df] = '~rouxf/prj/TC/matFiles/';

% %%
% Fs = 1.2e3;
% cfg                     = [];
% cfg.dataset             = '/home/rouxf/prj/TC/THA07_DEVELOPMENT_20090924_12.ds';
% cfg.channel             = {'MEG'};
% 
% cfg.trl                 = [110*Fs 120*Fs 0];
% [dc1] = ft_preprocessing( cfg );
% 
% cfg.trl                 = [280*Fs 290*Fs 0];
% [dc2] = ft_preprocessing( cfg );
% 
% cfg                     = [];
% cfg.length              = 1;
% 
% [dc1] = ft_redefinetrial( cfg , dc1);
% 
% cfg                     = [];
% cfg.length              = 1;
% 
% [dc2] = ft_redefinetrial( cfg , dc2);
% 
% save(['resting_state_recording_THA07_MEGdataEOnoCleaning.mat'],'dc1');
% save(['resting_state_recording_THA07_MEGdataECnoCleaning.mat'],'dc2');

%% load anatomical volume data
% template grid (MNI template)
load([p2df,'template_volumes.mat'],'template_grid');

%load the forward model: inward shift:0, gird res: 1cm
dat = load([p2df,'individual_MNIwarped_grid.mat']);
sgrid = dat.grid;
tm = dat.tm;

dat = load([p2df,'individual_headmodel.mat']);
hdm = dat.hdm;
clear dat;

%% get the individual grid dimensions
%to align grid with anatomical template
Nx = length(template_grid.xgrid);
Ny = length(template_grid.ygrid);
Nz = length(template_grid.zgrid);

%%
ix = regexp(dc1.label,'[O]');
sel = [];
for it = 1:length(ix)
    if ~isempty(ix{it})
        sel = [sel it];
    end;
end;

for it = 1:length( dc1.time )
    
    dc1.time{it} = dc1.time{it} - min(  dc1.time{it} );
    dc2.time{it} = dc2.time{it} - min(  dc2.time{it} );
    
end;

% cfg = [];
% cfg.length = 1.0008;
% 
% [dum1] = ft_redefinetrial( cfg, dc1 );
% [dum2] = ft_redefinetrial( cfg, dc2 );

cfg                     = [];
cfg.channel             = dc2.label(sel);

[dum1] = ft_selectdata( cfg , dc1 );
[dum2] = ft_selectdata( cfg , dc2 );
[dum3] = ft_selectdata( cfg , dc3 );

Fs = dc1.fsample;

T = length(dc1.time{1})/Fs;
W = 1;%1/T;
TW = T*W;
k = floor(2*TW-1);

params                  = [];
params.Fs               = 1.2e3;
params.pad              = 0;
params.fpass            = [2 30];
params.tapers           = [TW k];
params.trialave         = 1;

S1 = {};
S2 = {};
S3 = {};
for it = 1:length(dum1.trial)
    fprintf([num2str(it),'/',num2str(size(dum1.trial,2))]);
    [S1{it},f] = mtspectrumc( (dum1.trial{it})' , params);
    [S2{it},f] = mtspectrumc( (dum2.trial{it})' , params);
    [S3{it},f] = mtspectrumc( (dum3.trial{it})' , params);
    fprintf('\n');
end;

x1 = zeros(size(S1{1}));
x2 = zeros(size(S2{1}));
x3 = zeros(size(S3{1}));
for it = 1:length( S1 )
    x1 = x1+S1{it};
    x2 = x2+S2{it};
    x3 = x3+S3{it};
end;
x1 = x1./it;
x2 = x2./it;
x3 = x3./it;

figure;
subplot(2,3,[1 2]);
hold on;
h = [];
h(1) = plot(f,20*log10(x1),'LineWidth',3);
h(2) = plot(f,20*log10(x2),'r','LineWidth',3);
h(3) = plot(f,20*log10(x3),'k','LineWidth',3);

axis tight;
xlabel('Frequency (Hz)');
ylabel('Spectral Power [dB]');
legend(h,'Eyes Open','Eyes Closed','Simulated');
legend 'boxoff';

cfg                     = [];
cfg.bpfilter            = 'yes';
cfg.bpfreq              = [8 12];

[dum1] = ft_preprocessing( cfg, dc1 );
[dum2] = ft_preprocessing( cfg, dc2 );
[dum3] = ft_preprocessing( cfg, dc3 );

for it = 1:length( dum1.trial )    
    dum1.trial{it} = dum1.trial{it}.^2;
    dum2.trial{it} = dum2.trial{it}.^2;
    dum3.trial{it} = dum3.trial{it}.^2;
end;

cfg                     = [];
cfg.avgovertime         = 'yes';

[dum1] = ft_selectdata( cfg , dum1 );
[dum2] = ft_selectdata( cfg , dum2 );
[dum3] = ft_selectdata( cfg , dum3 );

for it = 1:length( dum1.trial )
    
    dum1.trial{it} =sqrt( dum1.trial{it} );
    dum2.trial{it} =sqrt( dum2.trial{it} );
    dum3.trial{it} =sqrt( dum3.trial{it} );
    
end;

cfg                     = [];
cfg.keeptrials          = 'no';

[dum1] = ft_timelockanalysis( cfg ,dum1 );
[dum2] = ft_timelockanalysis( cfg ,dum2 );
[dum3] = ft_timelockanalysis( cfg ,dum3 );

cfg                     = [];
cfg.layout              = 'CTF275.lay';
cfg.comment             = 'no';
cfg.marker              = 'off';
%cfg.highlight           = 'on';
cfg.highlightchannel    = sel;
cfg.highlightcolor      = [1 1 1];
cfg.highlightsymbol     = 'o';
cfg.highlightsize       = 8;
cfg.highlightlinewidth  = 3;

subplot(234);
a =gca;
ft_topoplotER( cfg , dum1 );
ca1 = caxis;
cb1 = colorbar;
title('Eyes Open');
subplot(235);
a = [a gca];
ft_topoplotER( cfg , dum2 );
ca2 = caxis;
cb2 = colorbar;
title('Eyes Closed');
subplot(236);
a = [a gca];
ft_topoplotER( cfg , dum3 );
ca3 = caxis;
cb3 = colorbar;
title('Simulated');

for it = 1:length(a)
    caxis(a(it),ca3);
end;
set(cb1,'YTick',get(cb1,'YTick'));
set(cb2,'YTick',get(cb2,'YTick'));
set(cb1,'YTickLabel',get(cb1,'YTick')./1e-15);
set(cb2,'YTickLabel',get(cb2,'YTick')./1e-15);
set(cb3,'YTickLabel',get(cb3,'YTick')./1e-15);
set(get(cb1,'YLabel'),'String','RMS [fT]');
set(get(cb2,'YLabel'),'String','RMS [fT]');

%%
cfg                     = [];
cfg.bpfilter            = 'yes';
cfg.bpfreq              = [8 12];

[dum1] = ft_preprocessing( cfg, dc1 );
[dum2] = ft_preprocessing( cfg, dc2 );
[dum3] = ft_preprocessing( cfg, dc3 );

cfg                     = [];
cfg.covariance          = 'yes';

[cov1] = ft_timelockanalysis( cfg ,dum1 );
[cov2] = ft_timelockanalysis( cfg ,dum2 );
[cov3] = ft_timelockanalysis( cfg ,dum3 );
clear dum*;

cfg = [];
cfg.grid = sgrid;
cfg.headmodel = hdm;
cfg.senstype = 'meg';
cfg.method = 'lcmv';
cfg.normalize = 'yes';

%cfg.grid.pos = Gridcoord;
%cfg.grid.mom = [-1 0 0];

cfg.lcmv.fixedori = 'yes';% let dipole rotate freely = no (or not = yes)
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.projectmom = 'no';
cfg.lcmv.projectnoise = 'yes';
%cfg.lcmv.lambda = '10%';
cfg.lcmv.lambda = '1%';

cfg.grid.units = 'mm';
cfg.grid.dim=[Nx Ny Nz];

cfg.grad = cov1.grad;
[lcmv1] = ft_sourceanalysis(cfg,cov1);

cfg.grad = cov2.grad;
[lcmv2] = ft_sourceanalysis(cfg,cov2);

cfg.grad = cov3.grad;
[lcmv3] = ft_sourceanalysis(cfg,cov3);

clear cov*;

%%
lcmv1 = rmfield(lcmv1,'cfg');% get rid of large cfg-info
lcmv2 = rmfield(lcmv2,'cfg');% get rid of large cfg-info
lcmv3 = rmfield(lcmv3,'cfg');% get rid of large cfg-info

% align beamformer with grid-space
lcmv1.pos = template_grid.pos;
lcmv1.dim = template_grid.dim;
lcmv2.pos = template_grid.pos;
lcmv2.dim = template_grid.dim;
lcmv3.pos = template_grid.pos;
lcmv3.dim = template_grid.dim;

%%
load([p2df,'template_volumesII.mat']);

dum = lcmv1;
dum.avg.pow = 20*log10(lcmv2.avg.pow./lcmv1.avg.pow);

dum2 = lcmv1;
dum2.avg.pow = 20*log10(lcmv3.avg.pow./lcmv1.avg.pow);

cfg                     = [];
cfg.parameter           = 'avg.pow';

[intPow1] = ft_sourceinterpolate( cfg , dum , template_mri);
[intPow2] = ft_sourceinterpolate( cfg , dum2 , template_mri);

cfg                     = [];
cfg.funparameter        ='avg.pow';
cfg.maskparameter       = cfg.funparameter;
cfg.opacitymap          = 'rampup';
%cfg.opacitylim          = [40 100];
cfg.crosshair           = 'no';
cfg.interactive         = 'no';

figure;
ft_sourceplot( cfg , intPow1 );

figure;
ft_sourceplot( cfg , intPow2 );