%%
if matlabpool('size') ==0
    matlabpool(24);
end;
%%
restoredefaultpath;
addpath('~froux/froux/fieldtrip-20151020/');
ft_defaults;
savepath = '/bcbl/home/home_a-f/froux/dipole_simulation_thalamicPAC/matFiles/';
%%
tl = 1:30;

Fs = 300:100:1200;
pfoi = 4:2:20;
afoi = 45:5:140;

dfp = mean(diff(pfoi));
dfa = mean(diff(afoi));

pbins = -pi:pi/8:pi;
%%  
MI = zeros(length(Fs),length(tl),length(pfoi),length(afoi));
SNR = zeros(length(Fs),length(tl));
parfor it = 1:length(Fs)
    snr = zeros(length(tl),1);
   [mi3] = zeros(length(tl),length(pfoi),length(afoi));
   for jt = 1:length(tl) 
    
        t = 0:1/Fs(it):tl(jt);
        
        asig = sin(2*pi*10.*t);%
        gsig = (sin(2*pi*10.*t)+1).*sin(2*pi*70.*t); 
        
        sig = asig+gsig+1;%
        noise = (mean(sig)+2*std(sig).*randn(1,length(t)));
        
        snr(jt) = 10*log10(mean(sig)/mean(noise));
        
        sig = sig + noise;
        
        dum = [];
        dum.label = {'dummy_channel'};
        dum.cfg = [];
        dum.trial{1} = sig;
        dum.time{1} = t;
        
        [mi2] = zeros(length(pfoi),length(afoi));
        for kt = 1:length(pfoi)
            
            cfg = [];
            cfg.hilbert = 'angle';
            cfg.bpfilter = 'yes';
            cfg.bpfreq = [pfoi(kt)-dfp pfoi(kt)+dfp];
            cfg.bpfilttype = 'fir';
            
            [phi] = ft_preprocessing(cfg,dum);
            p = phi.trial{1};
            
            [mi] = zeros(length(afoi),1);
            for mt = 1:length(afoi)
                
                cfg = [];
                cfg.hilbert = 'abs';
                cfg.bpfilter = 'yes';
                cfg.bpfreq = [afoi(mt)-dfa afoi(mt)+dfa];
                cfg.bpfilttype = 'fir';
            
                [amp] = ft_preprocessing(cfg,dum);
                a = amp.trial{1}.^2;
                
                PAH = zeros(1,length(pbins));
                for zt = 1:length(pbins)-1
                    a2 = a;
                    pbins2 = pbins;
                    pidx = find( p >= pbins2(zt) & p < pbins2(zt+1));
                    
                    PAH(zt) = mean(a2(pidx));
                    
                end;
                PAH(end) = PAH(1);
                
                n = length(pbins);
                PAH = PAH./sum(PAH,2);
                H = -sum(log(PAH).*PAH,2);
                mi(mt) = (log(n)-H)./log(n);
                
            end;
            mi2(kt,:) = mi;
        end;
        mi3(jt,:,:) = mi2;
        
        
   end;
   MI(it,:,:,:) = mi3;
   SNR(it,:) = snr;
end;
matlabpool close;
%%
save([savepath,'parameters4MI_alphaGammaCoupling.mat'],'MI','SNR','pbins','afoi','pfoi','Fs','tl');