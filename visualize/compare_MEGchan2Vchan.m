restoredefaultpath;
addpath('~froux/froux/fieldtrip-20151020/');
ft_defaults;
addpath(genpath('~froux/froux/dipole_simulation_thalamicPAC/mcode/'));
p2df = '~froux/froux/dipole_simulation_thalamicPAC/matFiles/';
%%
load([p2df,'parietal_dipoles_alphagammaPAC.mat']);
load([p2df,'virtual_channels_1:150Hz.mat']);

cfg = [];
cfg.keeptrials = 'yes';

[raw1] = ft_timelockanalysis(cfg,raw1);
[VC] = ft_timelockanalysis(cfg,VC);
%%
t = 0:1/1200:raw1.time(end);
dip_sig = ((sin(2*pi*10.*t)+1).*sin(2*pi*70.*t)).*0.25+sin(2*pi*10.*t)+(1.*randn(1,length(t)));

dip_sig = (dip_sig - min(dip_sig))./(max(dip_sig)-min(dip_sig));
%% load anatomical volume data
load([p2df,'individual_MRIdata.mat']);

load([p2df,'template_volumes.mat'],'template_grid');
%% load spatial filter
load([p2df,'lcmv_spatial_filter_1:150Hz.mat']);

[lcmv] = ft_convert_units(lcmv,'mm');
%% normalize individual T1 to template brain
[spm_template] = '/bcbl/home/home_a-f/froux/fieldtrip-20151020/external/spm8/templates/T1.nii';

cfg = [];
cfg.spmversion  = 'spm8';
cfg.template = spm_template;
cfg.coordinates = 'ctf';
cfg.nonlinear = 'no';

[norm] = ft_volumenormalise(cfg,mri);
%% compute the coordinates in CTF sensor space for a dipole in the
% left and right occipital cortex
idx = find(lcmv.inside==1);

[lcmv.pos(idx(VCidx),:);warp_apply(norm.cfg.final,Gridcoord,'homogenous')]

[warp_apply(inv(norm.cfg.final),lcmv.pos(idx(VCidx),:),'homogenous');Gridcoord]
%%
cfg = [];
cfg.channel = VC.label(VCidx);

[VCsig] = ft_selectdata(cfg,VC);

vc_sig = zeros(size(VCsig.trial,1),1,length(VCsig.time));
for it = 1:size(VCsig.trial,1)
    x = squeeze(VCsig.trial(it,1,:));
    vc_sig(it,1,:) = (x -min(x))./(max(x)-min(x));
end;

VCsig.trial = vc_sig;
%%
x = squeeze(vc_sig(15,1,:));

figure;
subplot(411);
hold on;
plot(t,dip_sig,'k','LineWidth',3);
plot(t,x,'r')
xlim([4 5]);
title('Unfiltered signals');
legend('Simulated dip.1','Virtual channel');
xlabel('Time [s]');
ylabel('Amplitude');

rxy = zeros(size(vc_sig,1),1);
PLV = zeros(size(vc_sig,1),length(dip_sig));

for it = 1:size(vc_sig,1)
    x = squeeze(vc_sig(it,1,:));
    rxy(it) = corr(dip_sig',x);
    PLV(it,:) = exp(1i.*(angle(hilbert(x))-angle(hilbert(dip_sig))'));
end;

PLV = 1/size(PLV,1)*abs(sum(PLV,1));

text(4.2,max(dip_sig),['rxy:',num2str(round(mean(rxy)*100)/100),', PLV:',num2str(round(mean(PLV)*100)/100)]);

dum1 = [];
dum1.trial{1} = dip_sig;
dum1.time{1} = t;
dum1.label{1} = 'dip_sig';

dum2 = [];
dum2.trial = cell(1,size(VCsig.trial,1));
dum2.time = cell(1,size(VCsig.time,1));
for it = 1:size(VCsig.trial,1)
    dum2.trial{it} = squeeze(VCsig.trial(it,1,:))';
    dum2.time{it} = t;
end;
dum2.label{1} = 'VC_sig';

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfilttype = 'fir';
cfg.bpfreq = [9 11];
cfg.padding = 100*max(dum1.time{1});
cfg.padtype = 'mirror';

[alpha1] = ft_preprocessing(cfg,dum1);
[alpha2] = ft_preprocessing(cfg,dum2);

subplot(412);
hold on;
plot(t,alpha1.trial{1},'k','LineWidth',3);
plot(t,alpha2.trial{1},'r');
xlim([4 5]);
title('Alpha-band signals (9-11Hz)');
xlabel('Time [s]');
ylabel('Amplitude');

rxy = zeros(length(alpha2.trial),1);
PLV = zeros(length(alpha2.trial),length(dip_sig));
for it = 1:length(alpha2.trial)
    rxy(it) = corr(alpha1.trial{1}',alpha2.trial{it}');
    PLV(it,:) = exp(1i.*(angle(hilbert(alpha1.trial{1}'))-angle(hilbert(alpha2.trial{it}'))));
end;

PLV = 1/size(PLV,1)*abs(sum(PLV,1));

text(4.2,max(alpha1.trial{1}),['rxy:',num2str(round(mean(rxy)*100)/100),', PLV:',num2str(round(mean(PLV)*100)/100)]);

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfilttype = 'fir';
cfg.bpfreq = [69 71];
cfg.padding = 100*max(dum1.time{1});
cfg.padtype = 'mirror';

[gamma1] = ft_preprocessing(cfg,dum1);
[gamma2] = ft_preprocessing(cfg,dum2);

subplot(413);
hold on;
plot(t,gamma1.trial{1},'k','LineWidth',3);
plot(t,gamma2.trial{1},'r');
xlim([4 5]);
title('Gamma-band signals (55-85Hz)');
xlabel('Time [s]');
ylabel('Amplitude');

rxy = zeros(length(gamma2.trial),1);
PLV = zeros(length(gamma2.trial),length(dip_sig));
for it = 1:length(gamma2.trial)
    rxy(it) = corr(gamma1.trial{1}',gamma2.trial{it}');
    PLV(it,:) =exp(1i.*(angle(hilbert(gamma1.trial{1}'))-angle(hilbert(gamma2.trial{it}'))));
end;

PLV = 1/size(PLV,1)*abs(sum(PLV,1));

text(4.2,max(gamma1.trial{1}),['rxy:',num2str(round(mean(rxy)*100)/100),' ,PLV:',num2str(round(mean(PLV)*100)/100)]);

cfg = [];
cfg.hilbert = 'abs';

[gamma1] = ft_preprocessing(cfg,gamma1);
[gamma2] = ft_preprocessing(cfg,gamma2);

cfg = [];
cfg.method = 'mtmfft';
cfg.pad = 'maxperlen';
cfg.output = 'pow';
cfg.taper = 'dpss';
cfg.tapsmofrq = .5;

[pow1] = ft_freqanalysis(cfg,gamma1);
[pow2] = ft_freqanalysis(cfg,gamma2);

subplot(414);
hold on;
plot(pow1.freq,pow1.powspctrm,'k','LineWidth',3);
plot(pow2.freq,pow2.powspctrm,'r');
xlabel('Frequency [Hz]');
ylabel('Power of gamma-env. [a.u.]');
xlim([0 40]);

set(gcf,'Color','w');
%%
dum = [];
dum.label = [raw1.label(:);{'dipole_1'}];
dum.trial = cell(1,size(raw1.trial,1));
for it = 1:size(raw1.trial,1)
    dum.trial{it} = zeros(length(raw1.label)+1,length(raw1.time));
end;

for it = 1:size(raw1.trial,1)
    dum.trial{it}(1:length(raw1.label),:) = squeeze(raw1.trial(it,:,:));
    dum.trial{it}(length(raw1.label)+1,:) = dip_sig;% + (2*std(dip_sig).*randn(1,length(dip_sig)));
    dum.time{it}  = raw1.time;
end;

cfg =[];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.pad = 'maxperlen';
cfg.taper = 'dpss';
cfg.tapsmofrq = 1;
cfg.keeptrials = 'yes';
cfg.foi = 0:100;

[phi] = ft_freqanalysis(cfg,dum);

cfg = [];
cfg.method = 'coh';
cfg.channel = dum.label;
cfg.channelcmb = [repmat(dum.label(end),[length(dum.label) 1]) dum.label];

[coh1] = ft_connectivityanalysis(cfg,phi);

dum = [];
dum.label = [VC.label(:);{'dipol_1'}];
dum.trial = cell(1,size(VC.trial,1));
for it = 1:size(VC.trial,1)
    dum.trial{it} = zeros(length(VC.label)+1,length(VC.time));
end;

for it = 1:size(VC.trial,1)
    dum.trial{it}(1:length(VC.label),:) = squeeze(VC.trial(it,:,:));
    dum.trial{it}(length(VC.label)+1,:) = dip_sig;% + (2*std(dip_sig).*randn(1,length(dip_sig)));
    dum.time{it}  = VC.time; 
end;

cfg =[];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.pad = 'maxperlen';
cfg.taper = 'dpss';
cfg.tapsmofrq = 1;
cfg.keeptrials = 'yes';
cfg.foi = 0:100;

[phi] = ft_freqanalysis(cfg,dum);

cfg = [];
cfg.method = 'coh';
cfg.channel = dum.label;
cfg.channelcmb = [repmat(dum.label(end),[length(dum.label) 1]) dum.label];

[coh2] = ft_connectivityanalysis(cfg,phi);
%%
RMS = squeeze(mean(sqrt(raw1.trial.^2),3));
[v,i1] = max(max(RMS));

%[i1,~] = find(coh1.cohspctrm == max(max(coh1.cohspctrm)));

meg_sig = squeeze(raw1.trial(15,i1,:));
meg_sig = (meg_sig - min(meg_sig))./(max(meg_sig)-min(meg_sig));
%% PAC for VC
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [9 11];
cfg.bpfilttype = 'fir';
cfg.hilbert = 'angle';% phase in rad
cfg.padding = 3*max(VC.time);%
cfg.padtype = 'mirror';
cfg.continuous = 'no';
cfg.channel = VC.label(VCidx);

[alpha] = ft_preprocessing(cfg,VC);% single trial alpha activity

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [65 75];
cfg.bpfilttype = 'fir';
cfg.hilbert = 'abs'; % compute the power
cfg.padding = 3*max(VC.time);%
cfg.padtype = 'mirror';
cfg.continuous = 'no';
cfg.channel = VC.label(VCidx);

[gamma] = ft_preprocessing(cfg,VC);

% concatenate
concat1 = zeros(length(alpha.label),length(alpha.time)*size(alpha.trial,1));
concat2 = zeros(length(gamma.label),length(gamma.time)*size(gamma.trial,1));

idx = 1:length(alpha.time);
for it = 1:size(alpha.trial,1)
    concat1(:,idx) = squeeze(alpha.trial(it,1,:));
    concat2(:,idx) = squeeze(gamma.trial(it,1,:)).^2;% square the amp to obtain power
    idx = idx + length(alpha.time);
end;
clear alpha gamma;

% phase bins
[pbins] = -pi:pi/8:pi;
%preallocate
[PAH] = zeros(size(concat1,1),length(pbins));%% loop over virtual channels    
for it = 1:size(concat1,1)
    
    X = zeros(length(pbins),1);
    for kt = 1:length(pbins)-1% loop over phase bins
        
        [phi] = concat1(it,:);
        [amp] = concat2(it,:);
        
        
        idx = [];
        [idx] = find(phi  >= pbins(kt) & phi < pbins(kt+1));
        
        X(kt) = mean(amp(idx));
        
    end;
    X(end) = X(1);
    PAH(it,:) = X;
end;
% close parallel pool
clear concat*;
% compute MI
PAH1 = PAH./repmat(sum(PAH,2),[1 size(PAH,2)]);
H = -sum(log(PAH1).*PAH1,2);
n = length(pbins);
MI1 = (log(n)-H)./log(n);
%% PAC for simulated dipole
dum = [];
dum.time = VC.time;
dum.label = {'dipol1'};
dum.trial = zeros(size(VC.trial,1),1,length(VC.time));
for it = 1:size(VC.trial,1)    
   dum.trial(it,:)=  dip_sig;
end;
dum.dimord = VC.dimord;

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [9 11];
cfg.bpfilttype = 'fir';
cfg.hilbert = 'angle';% phase in rad
cfg.padding = 3*max(dum.time);%
cfg.padtype = 'mirror';
cfg.continuous = 'no';

[alpha] = ft_preprocessing(cfg,dum);% single trial alpha activity

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [65 75];
cfg.bpfilttype = 'fir';
cfg.hilbert = 'abs'; % compute the power
cfg.padding = 3*max(dum.time);%
cfg.padtype = 'mirror';
cfg.continuous = 'no';

[gamma] = ft_preprocessing(cfg,dum);

% concatenate
concat1 = zeros(length(alpha.label),length(alpha.time)*size(alpha.trial,1));
concat2 = zeros(length(gamma.label),length(gamma.time)*size(gamma.trial,1));

idx = 1:length(alpha.time);
for it = 1:size(alpha.trial,1)
    concat1(:,idx) = squeeze(alpha.trial(it,1,:));
    concat2(:,idx) = squeeze(gamma.trial(it,1,:)).^2;% square the amp to obtain power
    idx = idx + length(alpha.time);
end;
clear alpha gamma;

% phase bins
[pbins] = -pi:pi/8:pi;
%preallocate
[PAH] = zeros(size(concat1,1),length(pbins));%% loop over virtual channels    
for it = 1:size(concat1,1)
    
    X = zeros(length(pbins),1);
    for kt = 1:length(pbins)-1% loop over phase bins
        
        [phi] = concat1(it,:);
        [amp] = concat2(it,:);
        
        
        idx = [];
        [idx] = find(phi  >= pbins(kt) & phi < pbins(kt+1));
        
        X(kt) = mean(amp(idx));
        
    end;
    X(end) = X(1);
    PAH(it,:) = X;
end;
% close parallel pool
clear concat*;
% compute MI
PAH2 = PAH./repmat(sum(PAH,2),[1 size(PAH,2)]);
H = -sum(log(PAH2).*PAH2,2);
n = length(pbins);
MI2 = (log(n)-H)./log(n);
%% PAC for MEG channels
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [9 11];
cfg.bpfilttype = 'fir';
cfg.hilbert = 'angle';% phase in rad
cfg.padding = 3*max(VC.time);%
cfg.padtype = 'mirror';
cfg.continuous = 'no';
cfg.channel = raw1.label(i1);

[alpha] = ft_preprocessing(cfg,raw1);% single trial alpha activity

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [65 75];
cfg.bpfilttype = 'fir';
cfg.hilbert = 'abs'; % compute the power
cfg.padding = 3*max(VC.time);%
cfg.padtype = 'mirror';
cfg.continuous = 'no';
cfg.channel = raw1.label(i1);

[gamma] = ft_preprocessing(cfg,raw1);

% concatenate
concat1 = zeros(length(alpha.label),length(alpha.time)*size(alpha.trial,1));
concat2 = zeros(length(gamma.label),length(gamma.time)*size(gamma.trial,1));

idx = 1:length(alpha.time);
for it = 1:size(alpha.trial,1)
    concat1(:,idx) = squeeze(alpha.trial(it,1,:));
    concat2(:,idx) = squeeze(gamma.trial(it,1,:)).^2;% square the amp to obtain power
    idx = idx + length(alpha.time);
end;
clear alpha gamma;

% phase bins
[pbins] = -pi:pi/8:pi;
%preallocate
[PAH] = zeros(size(concat1,1),length(pbins));%% loop over virtual channels    
for it = 1:size(concat1,1)
    
    X = zeros(length(pbins),1);
    for kt = 1:length(pbins)-1% loop over phase bins
        
        [phi] = concat1(it,:);
        [amp] = concat2(it,:);
        
        
        idx = [];
        [idx] = find(phi  >= pbins(kt) & phi < pbins(kt+1));
        
        X(kt) = mean(amp(idx));
        
    end;
    X(end) = X(1);
    PAH(it,:) = X;
end;
% close parallel pool
clear concat*;
% compute MI
PAH3 = PAH./repmat(sum(PAH,2),[1 size(PAH,2)]);
H = -sum(log(PAH3).*PAH3,2);
n = length(pbins);
MI3 = (log(n)-H)./log(n);
%%
dat = load([p2df,'alphaGammaPAC_source_level_fixedPosAndMom.mat']);
%%
figure;

subplot(5,2,[1 3]);

dum = [];
dum.label = raw1.label;
dum.dimord = 'chan_time';
dum.time = 1;
dum.avg = sqrt(squeeze(mean(mean(sqrt(abs(raw1.trial).^2),3),1)))';

cfg = [];
cfg.parameter = 'avg';
cfg.layout = 'CTF275.lay';
cfg.marker = 'off';
cfg.comment = 'no';
cfg.highlightchannel   = i1;
cfg.highlight          = 'on';
cfg.highlightsymbol    = 'o';
cfg.highlightcolor     = [1 1 1];

ft_topoplotER(cfg,dum);
title('RMS');
cb = colorbar;

subplot(5,2,[2 4]);
cfg = [];
cfg.layout = 'CTF275.lay';
cfg.parameter = 'cohspctrm';
cfg.xlim = [10 10];
cfg.comment = 'no';
cfg.marker = 'off';
cfg.refchannel = coh1.labelcmb(1,1);
cfg.highlightchannel   = i1;
cfg.highlight          = 'on';
cfg.highlightsymbol    = 'o';
cfg.highlightcolor     = [1 1 1]; 

ft_topoplotTFR(cfg,coh1);
title('Coherence MEG-dip.1');
cb = colorbar;

subplot(525);
MNIcoord = ft_warp_apply(norm.cfg.final,Gridcoord,'homogenous');

idx = find(lcmv.inside==1);
dum = lcmv;
dum.avg.pow = NaN(size(lcmv.avg.pow));
dum.avg.pow(idx) = sqrt(squeeze(mean(mean(sqrt(VC.trial.^2),3),1)));
%dum.avg.pow(idx(i2)) = 10e4;
dum.avg.pow = normalize_range(dum.avg.pow);

cfg = [];
cfg.parameter = 'avg.pow';

[int] = ft_sourceinterpolate(cfg,dum,norm);

Lm = max(max(max(int.pow)));

ft_plot_slice(int.anatomy,'orientation',[0 0 1]);view([0 90]);
ft_plot_slice(int.pow,'datmask',int.pow,'opacitylim',[0.1*Lm Lm],'colormap','jet','orientation',[0 0 1]);view([0 90]);
axis tight;axis off;

title('RMS');
%cb = colorbar;
axis off;

subplot(526);
cfg = [];
cfg.frequency = [9.9 10.1];

[X] = ft_selectdata(cfg,coh2);

idx = find(lcmv.inside==1);
dum = lcmv;
dum.avg.pow = NaN(size(lcmv.avg.pow));
dum.avg.pow(idx) = X.cohspctrm;
dum.avg.pow = normalize_range(dum.avg.pow);

cfg = [];
cfg.parameter = 'pow';

[int] = ft_sourceinterpolate(cfg,dum,norm);

Lm = max(max(max(int.pow)));

ft_plot_slice(int.anatomy,'orientation',[0 0 1]);view([0 90]);
ft_plot_slice(int.pow,'datmask',int.pow,'opacitylim',[0.1*Lm Lm],'colormap','jet','orientation',[0 0 1]);view([0 90]);
axis tight;axis off;

title('Coherence VCs-dip.1');
%cb = colorbar;
axis off;

subplot(527);
hold on;
plot(coh1.freq,coh1.cohspctrm,'Color',[.75 .75 .75]);
plot(coh1.freq,coh1.cohspctrm(i1,:),'r');
xlabel('Frequency [Hz]');
ylabel('Coherence [a.u.]');
title('MEG channels - dip.1');

subplot(528);
hold on;
plot(coh2.freq,coh2.cohspctrm,'Color',[.75 .75 .75]);
plot(coh2.freq,coh2.cohspctrm(VCidx,:),'r');
xlabel('Frequency [Hz]');
ylabel('Coherence [a.u.]');
title('VCs - dip.1');

subplot(529);
hold on;
plot(pbins,PAH2,'k','LineWidth',3);
plot(dat.pbins,dat.PAH,'--','Color',[.75 .75 .75],'LineWidth',3);
plot(pbins,PAH1,'r','LineWidth',3);
plot(pbins,PAH3,'b','LineWidth',3);
axis tight;

title('Phase-amplitude histogram');
xlabel('Alpha phase [rad]');
ylabel('Gamma amp. [a.u.]');
legend('Simulate dip.1','VC fixed','VC SVD','MEG-chan');

subplot(5,2,10);
Y1 = PAH1;
Y1 = (Y1-min(Y1))./(max(Y1)-min(Y1));

Y2 = PAH2;
Y2 = (Y2-min(Y2))./(max(Y2)-min(Y2));


Y3 = PAH3;
Y3 = (Y3-min(Y3))./(max(Y3)-min(Y3));

Y4 = dat.PAH;
Y4 = (Y4-min(Y4))./(max(Y4)-min(Y4));

hold on;
plot(pbins,Y2,'k','LineWidth',3);
plot(dat.pbins,Y4,'--','Color',[.75 .75 .75],'LineWidth',3);
plot(pbins,Y1,'r','LineWidth',3);
plot(pbins,Y3,'b','LineWidth',3);

axis tight;
title('Phase-amplitude histogram');
xlabel('Alpha phase [rad]');
ylabel('Gamma amp. [norm.]');

set(gcf,'Color','w');