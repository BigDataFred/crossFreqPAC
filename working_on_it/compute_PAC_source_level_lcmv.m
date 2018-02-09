addpath('~froux/froux/fieldtrip-20151020/');
ft_defaults;
%% load anatomical volume data
load('/bcbl/home/home_a-f/froux/Development/MRI/THA07_V2.mri_MRI4MNI.mat');
load('/bcbl/home/home_a-f/froux/Development/templates/template_GRID_1cmRes_0Indwardshift.mat','template_grid');
load('/bcbl/home/home_a-f/froux/Development/forwardMod/THA07_forward_model_WM_DEV_inward_shift_0_1cmRes.mat');
load('/bcbl/home/home_a-f/froux/dipole_simulation_thalamicPAC/parietal_dipoles_alphagammaPAC.mat');
%% get the individual grid dimensions
Nx = length(template_grid.xgrid);
Ny = length(template_grid.ygrid);
Nz = length(template_grid.zgrid);
% %% compute the spectrum of bandpass filtered MEG signals
% cfg = [];
% cfg.method= 'mtmfft';
% cfg.pad = 'maxperlen';
% cfg.taper = 'dpss';
% cfg.tapsmofrq = 1;
% 
% [pow] = ft_freqanalysis(cfg,raw1);
% figure;
% plot(pow.freq,squeeze(mean(pow.powspctrm,1)));
% clear pow*;
% %% bandpass filter the MEG sensor level data
% cfg = [];
% cfg.bpfilter = 'yes';
% cfg.bpfreq = [9 11];
% cfg.bpfilttype = 'fir';
% %cfg.bpfiltord = 3*fix(dum.fsample/cfg.bpfreq(1));
% cfg.hilbert = 'angle';
% cfg.padding = 15;%'maxperlen';
% cfg.padtype = 'mirror';
% cfg.continuous = 'no';
% 
% [alpha] = ft_preprocessing(cfg,raw1);% single trial alpha activity
% 
% cfg = [];
% cfg.bpfilter = 'yes';
% cfg.bpfreq = [65 75];
% cfg.bpfilttype = 'fir';
% %cfg.bpfiltord = 3*fix(dum.fsample/cfg.bpfreq(1));
% cfg.hilbert = 'abs';
% cfg.padding = 15;%'maxperlen';
% cfg.padtype = 'mirror';
% cfg.continuous = 'no';
% 
% [gamma] = ft_preprocessing(cfg,raw1);
% %%
% concat1 = zeros(length(alpha.label),length(alpha.time{1})*length(alpha.trial));
% concat2 = zeros(length(gamma.label),length(gamma.time{1})*length(gamma.trial));
% 
% idx = 1:length(alpha.time{1});
% for it = 1:length(alpha.trial)
%     concat1(:,idx) = alpha.trial{it};
%     concat2(:,idx) = gamma.trial{it}.^2;
%     idx = idx + length(alpha.time{it});
% end;
% clear alpha gamma;
% %%
% pbins = -pi:pi/4:pi;
% 
% [PAH] = zeros(size(concat1,1),length(pbins));
% 
% for it = 1:size(concat1,1)
%     
%     phi = concat1(it,:);            
%     
%     amp = concat2(it,:);
%    
%     for kt = 1:length(pbins)-1
%         
%         idx = [];
%         [idx] = find(phi  >= pbins(kt) & phi < pbins(kt+1));        
%         
%         PAH(it,kt) = mean(amp(idx));
%         
%     end;
%     PAH(it,kt+1) = PAH(it,1);
%     
% end;
% clear concat*;
% %%
% PAH = PAH./repmat(sum(PAH,2),[1 size(PAH,2)]);
% H = -sum(log(PAH).*PAH,2);
% n = length(pbins);
% MI = (log(n)-H)./log(n);
% %%
% dum = struct;
% dum.label = raw1.label;
% dum.time = 1;
% dum.avg = MI;
% dum.dimord = 'chan_time';
% 
% cfg = [];
% cfg.keeptrials = 'no';
% 
% [dum2] = ft_timelockanalysis(cfg,raw1);
% dum2.avg = abs(dum2.avg).^2;
% 
% cfg = [];
% cfg.layout = 'CTF275.lay';
% cfg.parameter = 'avg';
% cfg.comment = 'no';
% cfg.marker = 'off';
% 
% figure;
% subplot(121);
% ft_topoplotER(cfg,dum2);
% 
% subplot(122);
% ft_topoplotER(cfg,dum);
%% normalize individual T1 to template brain
cfg = [];
cfg.spmversion  = 'spm8';
cfg.template = '/bcbl/home/home_a-f/froux/spm8 2/templates/T1.nii';
cfg.coordinates = 'ctf';
cfg.nonlinear = 'no';

[norm] = ft_volumenormalise(cfg,mri);
%% compute NAI
load('~froux/froux/dipole_simulation_thalamicPAC/lcmv_spatial_filter_broadband.mat');
%% interpolate spatial filter on anatomy
cfg = [];
cfg.parameter = 'nai';

[int] = ft_sourceinterpolate(cfg,lcmv,norm);
%%
cfg = [];
cfg.method = 'slice';
cfg.funparameter = 'nai';
cfg.maskparameter = 'nai';
cfg.opacitymap = 'rampup';
cfg.funcolormap = 'jet';

ft_sourceplot(cfg,int);
clear int*;
% %%
% idx = find(lcmv1.inside==1);
% 
% cfg = [];
% cfg.keeptrials = 'no';
% 
% [tlck] = ft_timelockanalysis(cfg,VC);
% 
% 
% dum = lcmv1;
% dum.nai = NaN(size(lcmv1.avg.pow));
% dum.nai(idx) = mean(tlck.avg.^2,2);%
% 
% cfg = [];
% cfg.parameter = 'nai';
% 
% [int] = ft_sourceinterpolate(cfg,dum,norm);
% 
% cfg = [];
% cfg.method = 'slice';
% cfg.funparameter = 'nai';
% cfg.maskparameter = 'nai';
% cfg.opacitymap = 'rampup';
% cfg.funcolormap = 'jet';
% 
% ft_sourceplot(cfg,int);
% clear int;
% %% compute the spectrum of virtual channels
% cfg = [];
% cfg.method= 'mtmfft';
% cfg.pad = 'maxperlen';
% cfg.taper = 'dpss';
% cfg.tapsmofrq = 1;
% 
% [pow] = ft_freqanalysis(cfg,VC);
% 
% [i1,~] = find(pow.powspctrm == max(max(pow.powspctrm)));
% figure;
% hold on;
% plot(pow.freq,squeeze(mean(pow.powspctrm,1)));
% plot(pow.freq,squeeze(mean(pow.powspctrm(i1,:),1)),'r');
% clear pow*;
%% bandpass filter the virtual channel data
load('~froux/froux/dipole_simulation_thalamicPAC/virtual_channels_3:20Hz.mat');

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [9 11];
cfg.bpfilttype = 'fir';
%cfg.bpfiltord = 3*fix(dum.fsample/cfg.bpfreq(1));
cfg.hilbert = 'angle';
cfg.padding = max(VC.time{1});%'maxperlen';
cfg.padtype = 'mirror';
cfg.continuous = 'no';

[alpha] = ft_preprocessing(cfg,VC);% single trial alpha activity

load('~froux/froux/dipole_simulation_thalamicPAC/virtual_channels_50:90Hz.mat');

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [65 75];
cfg.bpfilttype = 'fir';
%cfg.bpfiltord = 3*fix(dum.fsample/cfg.bpfreq(1));
cfg.hilbert = 'abs';
cfg.padding = max(VC.time{1});%'maxperlen';
cfg.padtype = 'mirror';
cfg.continuous = 'no';

[gamma] = ft_preprocessing(cfg,VC);
clear VC;
%%
idx = find(lcmv.inside==1);

idxR = find(sign(lcmv.pos(:,1))==1);
idxL = find(sign(lcmv.pos(:,1))==-1);

idxR = intersect(idx,idxR);
idxL = intersect(idx,idxL);

m1 = lcmv.nai(idxR);
m2 = lcmv.nai(idxL);

m1 = idxR(find(m1 == max(m1)));
m2 = idxL(find(m2 == max(m2)));
%%
plvR = zeros(length(alpha.label),length(alpha.time{1}));
for jt = 1:length(alpha.trial)
    fprintf([num2str(jt),'/',num2str(length(alpha.trial))]);
    
    [ref1] = alpha.trial{jt}(find(idx ==m1),:);
    
    for it = 1:length(alpha.label)                
        
        [sig2] = alpha.trial{jt}(it,:);
        
        plvR(it,:) = plvR(it,:) + exp(1i*(ref1-sig2));
        
    end;
    fprintf('\n');
end;
N = length(alpha.trial);
plvR = abs(1/N.*plvR);

plvL = zeros(length(alpha.label),length(alpha.time{1}));
for jt = 1:length(alpha.trial)
    fprintf([num2str(jt),'/',num2str(length(alpha.trial))]);
    
    [ref1] = alpha.trial{jt}(find(idx == m2),:);
    
    for it = 1:length(alpha.label)                
        
        [sig2] = alpha.trial{jt}(it,:);
        
        plvL(it,:) = plvL(it,:) + exp(1i*(ref1-sig2));
        
    end;
    fprintf('\n');
end;
N = length(alpha.trial);
plvL = abs(1/N.*plvL);
%%
dum1 = lcmv;
dum1.nai = NaN(size(lcmv.avg.pow));
dum1.nai(idx) = ((mean(plvR,2))'+(mean(plvL,2))')./2;%
dum1.nai(m1) = 1000;
dum1.nai(m2) = 1000;

dum2 = lcmv;
dum2.nai = NaN(size(lcmv.avg.pow));
dum2.nai(idx) = ((mean(plvR,2))'+(mean(plvL,2))')./2;%

cfg = [];
cfg.parameter = 'nai';

[int1] = ft_sourceinterpolate(cfg,dum1,norm);
[int2] = ft_sourceinterpolate(cfg,dum2,norm);

cfg = [];
cfg.method = 'slice';
cfg.funparameter = 'nai';
cfg.maskparameter = 'nai';
cfg.opacitymap = 'rampup';
cfg.funcolormap = 'jet';
 
ft_sourceplot(cfg,int1);
ft_sourceplot(cfg,int2);
%%
idx = find(lcmv.inside==1);

x = zeros(length(lcmv.avg.filter{idx(1)}(:)),length(idx));
for it = 1:length(idx)
    x(:,it) = lcmv.avg.filter{idx(it)}(:);
end;

[ref1] = x(:,find(idx == m1));
[ref2] = x(:,find(idx == m2));

r1 = zeros(length(idx),1);
r2 = zeros(length(idx),1);
for it = 1:length(idx)
    
    [r1(it)] = corr(ref1,x(:,it),'Type','Pearson');
    [r2(it)] = corr(ref2,x(:,it),'Type','Pearson');
    
end;
%%
idx = find(lcmv.inside ==1);

dum1 = lcmv;
dum1.nai = NaN(size(lcmv.avg.pow));
dum1.nai(idx) = r1;%

dum2 = lcmv1;
dum2.nai = NaN(size(lcmv.avg.pow));
dum2.nai(idx) = r2;%

cfg = [];
cfg.parameter = 'nai';

[int1] = ft_sourceinterpolate(cfg,dum1,norm);
[int2] = ft_sourceinterpolate(cfg,dum2,norm);

cfg = [];
cfg.method = 'slice';
cfg.funparameter = 'nai';
cfg.maskparameter = 'nai';
cfg.opacitymap = 'rampup';
cfg.funcolormap = 'jet';
cfg.opacitylim = [.3 1];

ft_sourceplot(cfg,int1);caxis(cfg.opacitylim);
ft_sourceplot(cfg,int2);caxis(cfg.opacitylim);
%%
concat1 = zeros(length(alpha.label),length(alpha.time{1})*length(alpha.trial));
concat2 = zeros(length(gamma.label),length(gamma.time{1})*length(gamma.trial));

idx = 1:length(alpha.time{1});
for it = 1:length(alpha.trial)
    concat1(:,idx) = alpha.trial{it};
    concat2(:,idx) = gamma.trial{it}.^2;
    idx = idx + length(alpha.time{it});
end;
clear alpha gamma;
%%
idx = find(lcmv.inside ==1);
dum = lcmv;
dum.nai = NaN(size(lcmv.avg.pow));
dum.nai(idx) = mean(concat2.^2,2);%

cfg = [];
cfg.parameter = 'nai';

[int] = ft_sourceinterpolate(cfg,dum,norm);

cfg = [];
cfg.method = 'slice';
cfg.funparameter = 'nai';
cfg.maskparameter = 'nai';
cfg.opacitymap = 'rampup';
cfg.funcolormap = 'jet';
 
ft_sourceplot(cfg,int);
%%
pbins = -pi:pi/4:pi;

[PAH] = zeros(size(concat1,1),length(pbins));

for it = 1:size(concat1,1)
    
    phi = concat1(it,:);            
    
    amp = concat2(it,:);
   
    for kt = 1:length(pbins)-1
        
        idx = [];
        [idx] = find(phi  >= pbins(kt) & phi < pbins(kt+1));        
        
        PAH(it,kt) = mean(amp(idx));
        
    end;
    PAH(it,kt+1) = PAH(it,1);
    
end;
clear concat*;
%%
PAH = PAH./repmat(sum(PAH,2),[1 size(PAH,2)]);
H = -sum(log(PAH).*PAH,2);
n = length(pbins);
MI = (log(n)-H)./log(n);
%%
load('/bcbl/home/home_a-f/froux/dipole_simulation_thalamicPAC/alphaGammaPAC_source_level.mat');
load('/bcbl/home/home_a-f/froux/dipole_simulation_thalamicPAC/lcmv_spatial_filter_broadband.mat');

idx = find(lcmv.inside==1);

MIsource = lcmv;
MIsource.nai = NaN(size(lcmv.avg.pow));
MIsource.nai(idx) = MI;
%%
cfg = [];
cfg.parameter = 'nai';

[int3] = ft_sourceinterpolate(cfg,lcmv,norm);
%%
cfg = [];
cfg.method = 'slice';
cfg.funparameter = 'nai';
cfg.maskparameter = 'nai';
cfg.opacitymap = 'rampup';
%cfg.opacitylim = [max(MIsource.avg.MI)*0.6 max(MIsource.avg.MI)];%'maxabs';
cfg.funcolormap = 'jet';
 
ft_sourceplot(cfg,int3);
%%
[i1,~] = find(MI == max(max(MI)));
figure;
plot(pbins,PAH(i1,:),'k-o');
