function simulate_independent_alpha_generators(scf,slab)

if nargin <2
    slab = [];
end;

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
ntrl = 100;
t = 0:1/Fs:1;

f0 = 10;% frequency for alpha
f1 = [64:3:73];% frequency for gamma

x1 = [];
x2 = [];

phi1 = ([-105:1:-75]./180)*pi;% TC phase coupling
%phi1 = ([-15:1:15]./180)*pi;
phi2 = -pi:pi/100:pi;

[dat] = zeros(6,ntrl,length(t));

%%
for it = 1:ntrl% loop over trials
    
    fprintf([num2str(it),'/',num2str(ntrl)]);
    
    % initialize alpha signal
    a = zeros(size(dat,1),length(t));
    for jt = 1:size(dat,1)
        a(jt,:) = sin(2*pi*f0.*t);
        a(jt,:) = hilbert(a(jt,:));
    end;
    
    % TC coupling
    o1 = phi2(randperm(length(phi2)));
    o1 = o1(1);
    o2 = phi2(randperm(length(phi2)));
    o2 = o2(1);
    
    % TC coupling
    o3 = phi1(randperm(length(phi1)));
    o3 = o3(1);
    o4 = phi1(randperm(length(phi1)));
    o4 = o4(1);
    
    o5 = phi2(randperm(length(phi2)));
    o5 = o5(1);
    o6 = phi2(randperm(length(phi2)));    
    o6 = o6(1);
    
    a(1,:) = a(1,:)*exp(1i*0);% TH L
    a(2,:) = a(2,:)*exp(1i*0);% TH R
    a(3,:) = a(3,:)*exp(1i*o3);%BA7 L    
    %a(4,:) = a(4,:)*exp(1i*o4);
    a(5,:) = a(5,:)*exp(1i*o5);% BA18L
    %a(6,:) = a(6,:)*exp(1i*o6);
    a(4,:) = a(4,:)*exp(1i*o3);% BA7 R
    a(6,:) = a(6,:)*exp(1i*o5);% BA18 R
    
    %a(1,:) = a(1,:)*exp(1i*o1);
    %a(3,:) = a(3,:)*exp(1i*o1);
    %a(2,:) = a(2,:)*exp(1i*o2);
    %a(4,:) = a(4,:)*exp(1i*o2);
    
    a(1,:) = real( a(1,:) );
    a(2,:) = real( a(2,:) );
    a(3,:) = real( a(3,:) );    
    a(4,:) = real( a(4,:) );
    a(5,:) = real( a(5,:) );    
    a(6,:) = real( a(6,:) );
    
    fg = f1(randperm(length(f1)));
    
    g = zeros(size(dat,1),length(t));        
    g(3,:) = sin(2*pi*fg(1).*t).*(a(3,:)+1);
    g(4,:) = sin(2*pi*fg(2).*t).*(a(4,:)+1);
    
    g(5,:) = sin(2*pi*fg(3).*t).*(real(hilbert(a(5,:))*exp(1i*-pi/2))+1);
    g(6,:) = sin(2*pi*fg(4).*t).*(real(hilbert(a(6,:))*exp(1i*-pi/2))+1);
    
    o1 = phi2(randperm(length(phi2)));
    o1 = 0;%o1(1);
    o2 = phi2(randperm(length(phi2)));
    o2 = 0;%o1(1);
    o3 = phi2(randperm(length(phi2)));
    o4 = phi2(randperm(length(phi2)));
    o3 = 0;%o3(1);
    o4 = 0;%o4(1);
    
    g(3,:) = real( g(3,:)*exp(1i*o1) );
    g(4,:) = real( g(4,:)*exp(1i*o2) );
    g(5,:) = real( g(5,:)*exp(1i*o3) );
    g(6,:) = real( g(6,:)*exp(1i*o4) );
        
    x = a+g;    
    for jt = 1:size(x,1)
        n = mean(x(jt,:))+(std(x(jt,:))*2.5)*randn(1,length(t));
        x(jt,:) = x(jt,:)+n;
    end;
    
    dat(:,it,:) = x;
    
    fprintf('\n');
    
end;
% scf2(1) = scf(1)/max(max(max(dat(1:2,:,:))));
% scf2(2) = scf(2)/max(max(max(dat(3:4,:,:))));
% scf2(3) = scf(3)/max(max(max(dat(5:6,:,:))));
% 
% chck = [max(max(max(dat(1:2,:).*scf2(1)))) max(max(max(dat(3:4,:).*scf2(2)))) max(max(max(dat(5:6,:).*scf2(3))))];
% 
% if chck ~= scf
%     error('scaling factor is out of range');
% end;

scf2 = scf;

for it = 1:size(dat,1)    
    for jt = 1:size(dat,2)
        x = squeeze( dat(it,jt,:) );
        dat(it,jt,:) = x./max(abs(x));%(x -min(x))./(max(x)-min(x));
    end;
end;

dat([1:2],:,:) = dat([1:2],:,:).*scf2(1);
dat([3:4],:,:) = dat([3:4],:,:).*scf2(2);
dat([5:6],:,:) = dat([5:6],:,:).*scf2(3);

%dat([3:4],:,:) = dat([3:4],:,:)./100;
%dat([1:2 5:6],:,:) = dat([1:2 5:6],:,:)./10;

params = [];
params.Fs = Fs;
params.ntrl = ntrl;
params.t = t;

params.f0 = f0;
params.f1 = f1;

params.phi1 = phi1;
params.phi2 = phi2;

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

nfreq = length(params.fpass(1):params.fpass(2));

C = zeros(nfreq,size(dat,1),size(dat,1));
S1 = zeros(nfreq,size(dat,1),size(dat,1));
S2 = zeros(nfreq,size(dat,1),size(dat,1));

for it = 1:size(dat,1)
    for jt = 1:size(dat,1)
        [C(:,it,jt),~,~,S1(:,it,jt),S2(:,it,jt),f]=coherencyc(squeeze(dat(it,:,:))',squeeze(dat(jt,:,:))',params);
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
        x.t = t;        
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

