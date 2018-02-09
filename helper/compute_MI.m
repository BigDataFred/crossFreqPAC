function [MI,PAH] = compute_MI(sig,pfoi,afoi,pbins)
%%
dum = [];
dum.label = {'dummy_channel'};
dum.cfg = [];
dum.trial{1} = sig.sig;
dum.time{1} = sig.t;

[MI] = zeros(length(pfoi),length(afoi));
[X] = zeros(length(pfoi),length(afoi),length(pbins));

a = zeros(length(afoi),length(dum.trial{1}));
for mt = 1:length(afoi)
    
    cfg = [];
    cfg.hilbert = 'abs';
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [afoi(mt)-5 afoi(mt)+5];
    cfg.bpfilttype = 'fir';
    
    [amp] = ft_preprocessing(cfg,dum);
    
    a(mt,:) = amp.trial{1};%.^2;
end;

p = zeros(length(pfoi),length(dum.trial{1}));
for kt = 1:length( pfoi )
    cfg = [];
    cfg.hilbert = 'angle';
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [pfoi(kt)-1 pfoi(kt)+1];
    cfg.bpfilttype = 'fir';
    
    [phi] = ft_preprocessing(cfg,dum);
    p(kt,:) = phi.trial{1};
end;

parfor mt = 1:size(a,1)
    %fprintf([num2str(mt),'/',num2str(size(a,1))]);
    [pah] = zeros(length(pfoi),length(pbins));
    [mi] = zeros(length(pfoi),1);
    
    for kt = 1:length(pfoi)
        
        p2 = p(kt,:);
        
        PAH = zeros(1,length(pbins));
        for zt = 1:length(pbins)-1
            a2 = a(mt,:);
            pbins2 = pbins;            
            PAH(zt) = mean(a2(p2 >= pbins2(zt) & p2 < pbins2(zt+1)));
            
        end;
        PAH(end) = PAH(1);
        
        n = length(pbins);
        PAH = PAH./sum(PAH,2);
        H = -sum(log(PAH).*PAH,2);
        mi(kt,1) = (log(n)-H)./log(n);
        pah(kt,:) = PAH;
    end;
    MI(:,mt) = mi;
    X(:,mt,:) = pah;
    %fprintf('\n');
end;
PAH = X;


