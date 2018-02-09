%% set path environment
restoredefaultpath;
addpath('~froux/froux/fieldtrip-20151020/');
ft_defaults;
%% load the data
path2files = '/bcbl/home/home_a-f/froux/dipole_simulation_thalamicPAC/';
load([path2files,'PAC_source_alpha1.mat']);% spatial filters for alpha frequencies
load([path2files,'PAC_source_gamma1.mat']);% spatial filter for gamma frequencies
load([path2files,'parietal_dipoles_alphagammaPAC.mat']);% simulated MEG signals
%% load anatomical volume data
load('/bcbl/home/home_a-f/froux/Development/MRI/ABP24_V2.mri_MRI4MNI.mat');
%% normalize individual T1 to template brain
cfg = [];
cfg.spmversion  = 'spm8';
cfg.template = '/bcbl/home/home_a-f/froux/spm8 2/templates/T1.nii';
cfg.coordinates = 'ctf';
cfg.nonlinear = 'no';

[norm] = ft_volumenormalise(cfg,mri);
%% compute NAI for each single trial
for it = 1:length(source1)
    source1{it}.avg.nai = source1{it}.avg.pow./source1{it}.avg.noise;
    source3{it}.avg.nai = source3{it}.avg.pow./source3{it}.avg.noise;
end;
%% grand average of NAI across trials
cfg = [];
cfg.parameter = 'nai';

[sourceAVG1] = ft_sourcegrandaverage(cfg,source1{:});% alpha power
[sourceAVG2] = ft_sourcegrandaverage(cfg,source3{:});% gamma power
%% interpolate mean NAI to anatomy
cfg = [];
cfg.parameter = 'nai';

[int1] = ft_sourceinterpolate(cfg,sourceAVG1,norm);
[int2] = ft_sourceinterpolate(cfg,sourceAVG2,norm);
%% visualize NAI in source space
cfg = [];
cfg.method = 'slice';
cfg.funparameter = 'nai';
cfg.maskparameter = 'nai';
cfg.opacitymap = 'rampup';
cfg.funcolormap = 'jet';

cfg.opacitylim = [max(int1.nai)*0.35 max(int1.nai)];%'maxabs';
ft_sourceplot(cfg,int1);

cfg.opacitylim = [max(int2.nai)*0.35 max(int2.nai)];%'maxabs';
ft_sourceplot(cfg,int2);
%% bandpass filter the MEG sensor level data
cfg = [];
cfg.bpfilter ='yes';
cfg.bpfreq = [9 11];
cfg.bpfilttype = 'fir';
cfg.bpfiltord = 3*fix(raw1.fsample/cfg.bpfreq(1));

[alpha] = ft_preprocessing(cfg,raw1);% single trial alpha activity

cfg = [];
cfg.bpfilter ='yes';
cfg.bpfreq = [68 72];
cfg.bpfilttype = 'fir';
cfg.bpfiltord = 3*fix(raw1.fsample/cfg.bpfreq(1));% single trial gamma activity

[gamma] = ft_preprocessing(cfg,raw1);
%% compute the frequency spectrum of the raw and bp-fitlered data
cfg = [];
cfg.method = 'mtmfft';
cfg.foi = 0.5:100;
cfg.pad = 'maxperlen';
cfg.taper = 'dpss';
cfg.tapsmofrq = 1;

[freq1] = ft_freqanalysis(cfg,alpha);
[freq2] = ft_freqanalysis(cfg,gamma);
[freq3] = ft_freqanalysis(cfg,raw1);
%% visualize the spectra for the data at MEG sensor level
figure;
subplot(121);
hold on;
plot(freq1.freq,squeeze(mean(log(freq3.powspctrm),1)),'b');
plot(freq1.freq,squeeze(mean(log(freq1.powspctrm),1)),'r');
subplot(122);
hold on;
plot(freq1.freq,squeeze(mean(log(freq3.powspctrm),1)),'b');
plot(freq2.freq,squeeze(mean(log(freq2.powspctrm),1)),'r');
%% compute virtual channels in source space
%alpha= raw1;
%gamma = raw1;
idx = find(source1{1}.inside==1);

VC1 = alpha;
VC1.trial = cell(1,length(alpha.trial));
VC1.label = cell(1,length(idx));

VC2 = gamma;
VC2.trial = cell(1,length(gamma.trial));
VC2.label = cell(1,length(idx));
for it = 1:length(alpha.trial)
    fprintf([num2str(it),'/',num2str(length(alpha.trial))]);
    VC1.trial{it} = zeros(length(idx),length(alpha.time{1}));
    VC2.trial{it} = zeros(length(idx),length(gamma.time{1}));
    for jt = 1:length(idx)
        
        f = source1{it}.avg.filter{idx(jt)}(:);
        VC1.trial{it}(jt,:) = (f'*alpha.trial{it});% multiply MEG data with spatial filter
        VC1.label(jt) = {['virtual_channel',num2str(idx(jt))]};
        
        f = source3{it}.avg.filter{idx(jt)}(:);
        VC2.trial{it}(jt,:) = (f'*gamma.trial{it}); %multiply MEG data with spatial filter
        VC2.label(jt) = {['virtual_channel',num2str(idx(jt))]};

    end;
    fprintf('\n');
end;
%% compute the frequency spectra of the virtual channel data
cfg = [];
cfg.method = 'mtmfft';
cfg.foi = 1:100;
cfg.pad = 'maxperlen';
cfg.taper = 'dpss';
cfg.tapsmofrq = 1;

[pow1] = ft_freqanalysis(cfg,VC1);

[pow2] = ft_freqanalysis(cfg,VC2);

figure;
subplot(121);
plot(pow1.freq,squeeze(mean(pow1.powspctrm,1)));
axis tight;
subplot(122);
plot(pow2.freq,squeeze(mean(pow2.powspctrm,1)));
axis tight;
%%
cfg = [];
cfg.hilbert = 'angle';

[VC1] = ft_preprocessing(cfg,VC1);

cfg = [];
cfg.hilbert = 'abs';

[VC2] = ft_preprocessing(cfg,VC2);
for it = 1:length(VC2.trial)
    VC2.trial{it} = VC2.trial{it}.^2;
end;
%%
concat1 = zeros(length(VC1.label),length(VC1.time{1})*length(VC1.trial));
concat2 = zeros(length(VC2.label),length(VC2.time{1})*length(VC2.trial));

idx = 1:length(VC1.time{1});
for it = 1:length(VC1.trial)
    concat1(:,idx) = VC1.trial{it};
    concat2(:,idx) = VC2.trial{it};
    idx = idx + length(VC1.time{1});
end;
%%
pbins = -pi:pi/4:pi;

[PAH] = zeros(size(concat1,1),length(pbins));

for it = 1:size(concat1,1)
    
    phi = concat1(it,:);
    [phi,s_idx] = sort(phi);
    
    a = concat2(it,s_idx);
    
    for kt = 1:length(pbins)-1
        
        [idx] = find(phi  >= pbins(kt) & phi < pbins(kt+1));        
        PAH(it,kt) = mean(a(idx));
        
    end;
    PAH(it,kt+1) = PAH(it,1);
end;
%%
PAH = PAH./repmat(sum(PAH,2),[1 size(PAH,2)]);
H = -sum(log(PAH).*PAH,2);
n = length(pbins);
MI = (log(n)-H)./log(n);
%%
cfg = [];
cfg.output = 'fourier';
cfg.method = 'mtmconvol';
cfg.foi = 8:1:12;
cfg.toi =  VC1.time{1};
cfg.t_ftimwin = 6./cfg.foi;
cfg.pad = 'maxperlen';
cfg.taper = 'hanning';
cfg.tapsmofrq = 1;
cfg.keeptrials = 'yes';

phi = ft_freqanalysis(cfg,VC1);
phi.fourierspctrm = angle(phi.fourierspctrm);


cfg = [];
cfg.output = 'pow';
cfg.method = 'mtmconvol';
cfg.foi = 60:80;
cfg.toi = VC2.time{1};
cfg.t_ftimwin = 6./cfg.foi;
cfg.pad = 'maxperlen';
cfg.taper = 'dpss';
cfg.tapsmofrq = 15;
cfg.keeptrials = 'yes';

pow = ft_freqanalysis(cfg,VC2);
%%
concat1 = zeros(length(phi.label),length(phi.freq),length(phi.time)*size(phi.fourierspctrm,1));
idx = 1:length(phi.time);
for it = 1:size(phi.fourierspctrm,1)
    concat1(:,:,idx) = phi.fourierspctrm(it,:,:,:);
    idx = idx + length(phi.time);
end;
phi = rmfield(phi,'fourierspctrm');

concat2 = zeros(length(pow.label),length(pow.freq),length(pow.time)*size(pow.powspctrm,1));
idx = 1:length(pow.time);
for it = 1:size(pow.powspctrm,1)
    concat2(:,:,idx) = pow.powspctrm(it,:,:,:);
    idx = idx + length(pow.time);
end;
pow = rmfield(pow,'powspctrm');
%%
pbins = -pi:pi/4:pi;
PH = zeros(length(pow.label),length(phi.freq),length(pow.freq),length(pbins));
k = 0;
n_it = length(phi.freq)*length(pow.freq);
for it = 1:length(phi.freq)    
    for jt = 1:length(pow.freq)
        k = k+1;
        fprintf([num2str(k),'/',num2str(n_it)]);
        for lt = 1:length(phi.label)
            p = concat1(lt,it,:);
            del_idx = find(isnan(p));
            
            for kt = 1:length(pbins)-1
                
                idx = [];
                idx = find(p >= pbins(kt) & p <pbins(kt+1));

                if any(ismember(idx,del_idx))
                    idx(ismember(idx,del_idx)) = [];
                end;
                if isempty(idx)
                    error('phase bins cannot be empty');
                end;
                a = concat2(lt,jt,idx);
                a(isnan(a)) = [];
                PH(lt,it,jt,kt) = mean(a);
                
            end;
            PH(lt,it,jt,end) = PH(lt,it,jt,1);
        end;
        fprintf('\n');        
    end;
end;
clear concat*;
%%
X = PH./repmat(sum(PH,4),[1 1 1 size(PH,4)]);%normalize by sum over bins
H = -sum(X.*log(X),4);%shannon entropy
N = length(pbins);
MI = (log(N)-H)./log(N);% Kullback-Leibler Divergence
%%
MIsource = source1{1};
MIsource.avg.MI = NaN(size(MIsource.avg.pow));
idx = find(MIsource.inside==1);
MIsource.avg.MI(idx) =MI;% squeeze(mean(mean(MI,3),2));
%%
cfg = [];
cfg.parameter = 'MI';

[int3] = ft_sourceinterpolate(cfg,MIsource,norm);
%%
cfg = [];
cfg.method = 'slice';
cfg.funparameter = 'MI';
cfg.maskparameter = 'MI';
cfg.opacitymap = 'rampup';
%cfg.opacitylim = [max(MIsource.avg.MI)*0.9 max(MIsource.avg.MI)];%'maxabs';
cfg.funcolormap = 'jet';
 
ft_sourceplot(cfg,int3);
