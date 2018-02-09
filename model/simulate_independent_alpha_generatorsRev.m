function simulate_independent_alpha_generatorsRev(scf,slab)

if nargin <2
    slab = [];
end;

%%
rand('state',sum(100*clock));

%%
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));
addpath(genpath('/home/rouxf/prj/TC/mcode/'));

[savepath] = '/home/rouxf/prj/TC/matFilesRev/';

%%
if isempty(gcp('Nocreate'))
   parpool(36,'SpmdEnabled',false); 
end;

%%
Fs = 1.2e3;
ntrl = 600;
t = 1/Fs:1/Fs:1;% should be at least 1 min

f0 = [9.7:.1:11.3];% frequency for alpha
f1 = [64:3:73];% frequency for gamma

x1 = [];
x2 = [];

phi1 = ([-105:1:-75]./180)*pi;% TC phase coupling
%phi1 = ([-15:1:15]./180)*pi;

phi2 = -pi:pi/10000:pi;%random phase noise

[dat] = zeros(6,ntrl,length(t));

%%
for it = 1:ntrl% loop over trials
    
    fprintf([num2str(it),'/',num2str(ntrl)]);
    
    % initialize alpha signal
    a = zeros(size(dat,1),length(t));
    
    % draw random phase jitter
    o1 = phi2(randperm(length(phi2)));% large
    o2 = phi1(randperm(length(phi1)));% small
    o3 = phi2(randperm(length(phi2)));% large
        
    tmp1 = [];
    tmp2 = [];
    c1 = 0;
    c2 = 0;
    
    fa = f0(randperm(length(f0))); % draw random alpha frequency

    for jt = 1:size(dat,1)        
        if ismember(jt,[1 2])
            a(jt,:) = 2*sin(2*pi*fa(1).*t+o1(1));% this decouples or couples L&R thalamus                         
            c1 = c1+1;
            tmp1(c1) = o1(1);
            o1(1) = [];% comment this to sync
        elseif ismember(jt,[3 4])
            a(jt,:) = 2*sin(2*pi*fa(1).*t+tmp1(jt-2)+o2(1));% this couples BA7 to thalamus             
            o2(1) = [];% comment this to sync
        else
            c2 = c2+1;
            a(jt,:) = 6*sin(2*pi*fa(2).*t+o3(1));% this couples or decouples L&R BA 18   
            tmp2(c2) = o3(1);
            o3(1) = [];%comment this to sync
        end;                        
    end;        
    
    fg = f1(randperm(length(f1))); % draw random gamma frequency
    
    %set o to zero for perfect sync
    o = phi2(randperm(length(phi2)));%zeros(1,length(phi2));%        
    
    g = zeros(size(dat,1),length(t));        
    g(3,:) = .25*sin(2*pi*fg(1).*t+o(1)).*(1*(a(3,:)+1));
    g(4,:) = .25*sin(2*pi*fg(1).*t+o(2)).*(1*(a(4,:)+1));% set o(2) to o(1) to sync gamma
    
    g(5,:) = 2*sin(2*pi*fg(2).*t+o(3)).*(1*(sin(2*pi*fa(2).*t-pi/2+tmp2(1))+1));
    g(6,:) = 2*sin(2*pi*fg(2).*t+o(4)).*(1*(sin(2*pi*fa(2).*t-pi/2+tmp2(2))+1));% set o(4) to o(3) to sync gamma
    
    x = a+g;
    %add gaussian white noise to each time series
    for jt = 1:size(x,1)
        n = (2*std(x(jt,:)))*randn(1,size(x,2));
        x(jt,:) = x(jt,:)+n;
    end;
    
    dat(:,it,:) = x;
    
    fprintf('\n');
    
end;

%normalize amplitude range
for it = 1:size(dat,1)    
    for jt = 1:size(dat,2)
        x = squeeze( dat(it,jt,:) );
        dat(it,jt,:) = x./max(abs(x));%normalize from -1 to 1
    end;
end;

% apply a scaling factor if selected
dat([1:2],:,:) = dat([1:2],:,:).*scf(1);
dat([3:4],:,:) = dat([3:4],:,:).*scf(2);
dat([5:6],:,:) = dat([5:6],:,:).*scf(3);

% %% concatenate over time
ix = 1:5;
x2 = [];
for it = 1:size(dat,2)/5
    x1 = [];
    for jt = 1:length(ix)
        x1 = [x1 squeeze(dat(:,ix(jt),:))];
    end;
    ix = ix+5;
    x2(:,it,:) = x1;
end;
dat = [];
dat = x2;

%%
params = [];
params.Fs = Fs;

params.f0 = f0;
params.f1 = f1;

params.phi1 = phi1;
params.phi2 = phi2;

t = 1/params.Fs:1/params.Fs:size(dat,3)/params.Fs;
params.t = t;
params.ntrl = ntrl/(length(t)/Fs);

%%
r = zeros(size(dat,2),size(dat,1),size(dat,1));
parfor kt = 1:size(dat,2)
    fprintf([num2str(kt),'/',num2str(size(dat,2))]);
    dum = zeros(size(dat,1),size(dat,1));
    for it = 1:size(dat,1)
        
        [x] = squeeze(dat(it,kt,:));
        
        for jt = 1:size(dat,1)
            
            [y] = squeeze(dat(jt,kt,:));
            dum(it,jt) = corr(x,y,'type','Spearman');
            
        end;
        
    end;
    r(kt,:,:) = dum;
    fprintf('\n');
end;
r = squeeze(mean(r,1));
savename = ['independent_alphaGenerators_timeSeries',num2str(slab),'.mat'];
save([savepath,savename],'dat','params','r');

savename = ['independent_alphaGenerators_amplitudeCorrelationMatrix',num2str(slab),'.mat'];
save([savepath,savename],'r');

%%
TW = 3;
ntapers = 2*TW-1;
params                  =[];
params.Fs               = Fs;
params.pad              = -1;
params.tapers           = [TW,ntapers];
params.fpass            = [0 100];
params.trialave         =1;
params.err              = 0;

nfreq = length(params.fpass(1):1/(size(dat,3)/Fs):params.fpass(2));

C = zeros(nfreq,size(dat,1),size(dat,1));
S1 = zeros(nfreq,size(dat,1),size(dat,1));
S2 = zeros(nfreq,size(dat,1),size(dat,1));

for it = 1:size(dat,1)
    for jt = 1:size(dat,1)
        [C(:,it,jt),~,~,S1(:,it,jt),S2(:,it,jt),f]=coherencyc(squeeze(dat(it,:,:))',squeeze(dat(jt,:,:))',params);
        %[c,~,~,s1,s2,f]=coherencyc(squeeze(dat(it,:,:))',squeeze(dat(jt,:,:))',params);
    end;
end;
savename = ['independent_alphaGenerators_coherenceMatrixANDspectra',num2str(slab),'.mat'];

save([savepath,savename],'C','S1','S2','f');

%%
pfoi = 4:2:20;
afoi = 45:5:140;
pbins = -pi:pi/8:pi;

MI = zeros(size(dat,1),size(dat,2),length(pfoi),length(afoi));

PAH = cell(1,size(MI,1));

for it = 1:size(MI,1)
    X1 = zeros(size(dat,2),length(pfoi),length(afoi));
    X2 = zeros(size(dat,2),length(pfoi),length(afoi),length(pbins));
    parfor jt = 1:size(MI,2)
        fprintf([num2str(jt),'/',num2str(size(MI,2))])
        x = [];
        x.t = 1/Fs:1/Fs:size(dat,3)/Fs;        
        x.sig = squeeze(dat(it,jt,:))';
        [X1(jt,:,:),X2(jt,:,:,:)] = compute_MI(x,pfoi,afoi,pbins);
        fprintf('\n');
    end;
    MI(it,:,:,:) = X1;
    PAH{it}      = X2;
end;
delete(gcp);

savename = ['independent_alphaGenerators_modulationIndexANDpah',num2str(slab),'.mat'];

save([savepath,savename],'MI','PAH','pfoi','afoi','pbins');

