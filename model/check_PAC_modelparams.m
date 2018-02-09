%%
restoredefaultpath;
addpath('~froux/froux/fieldtrip-20151020/');
ft_defaults;
%%
Fs = 1.2e3;
t = 0:1/Fs:7;
fm = 10;
fc = 70;
B = 1.5;
A = 2.5;
[FM] = sin(2*pi*fm.*t);

[FC] = sin(2*pi*fc.*t);
%%
noise = 0+1.*randn(1,length(t));

Y(1,:) = B*FM+FC;
Y(2,:) = [1+B*FM].*FC;
Y(3,:) = A*FM+([1+B*FM].*FC);
%%
snr = zeros(size(Y,1),1);
for it = 1:size(Y,1)
    b = sqrt(mean(Y(it,:).^2));
    a = sqrt(mean(noise.^2));
    snr(it) = 10*log10(b/a);
    Y(it,:) = Y(it,:) +noise;
end;
%%
pfoi = 4:2:20;
afoi = 45:5:140;
pbins = -pi:pi/8:pi;

parpool(length(pfoi));

MI = zeros(size(Y,1),length(pfoi),length(afoi));
for it = 1:size(Y,1)    
    x.sig = Y(it,:);
    x.t = t;
    [MI(it,:,:)] = compute_MI(x,pfoi,afoi,pbins);
end;

delete(gcp);
%%
l = {'B*fm+fc','(1+B*fm)*fc','Afm+((1+B*fm)*fc)'};
%%
nfft = 2^nextpow2(size(Y,2));
%%
[pow] = fft(Y',nfft)./Fs;

[pow] = pow';
[pow] = pow.*conj(pow);

[pow] = pow(:,1:nfft/2+1);

f = Fs/2*linspace(0,1,nfft/2+1);
%%
a = zeros(size(Y,1),1);
b =zeros(size(Y,1),1);
c =zeros(size(Y,1),1);

figure;
for it = 1:size(Y,1)
    subplot(3,size(Y,1),it);
    a(it) = gca;
    plot(t,Y(it,:));
    title(l{it});
    yl = get(a(it),'YLim');
    text(.3,min(min(yl))+1,['SNR:',num2str(round(snr(it)*10)/10)]);

end;

yl2 = zeros(length(a),2);
for it = 1:size(pow,1)
    subplot(3,size(Y,1),size(pow,1)+it);
    b(it) = gca;
    plot(f,pow(it,:),'r');%./max(pow(it,:))
    yl2(it,:) = get(b(it),'YLim');      

    xlim([0 100]);
end;

ca = zeros(size(MI,1),2);
for it = 1:size(MI,1)
    subplot(3,size(Y,1),size(pow,1)*2+it);
    c(it) = gca;
    pcolor(pfoi,afoi,squeeze(MI(it,:,:))');
    shading interp;lighting phong;
    ca(it,:) = caxis;
end;

set(gcf,'Color','w');
set([a b],'LineWidth',1);
set([a b],'box','off');

yl = zeros(length(a),2);
for it = 1:length(a)
    xlabel(a(it),'Time [s]');
    ylabel(a(it),'Amp. [a.u.]');
    yl(it,:) = get(a(it),'YLim');      
    caxis(c(it),[min(min(ca)) max(max(ca))]);
end;

set(a,'YLim',[min(min(yl)) max(max(yl))]);
set(b,'YLim',[min(min(yl2)) max(max(yl2))]);

set(a,'XLim',[0 1]);
set(a,'XTick',[0 1]);

for it = 1:length(b)
    xlabel(b(it),'Frequency [Hz]');
    ylabel(b(it),'Power [a.u.]');
end;

set(b,'XTick',[fm fc]);
set(b,'XTickLabel',{'fm' 'fc'});


for it = 1:length(c)
    xlabel(c(it),'Frequency [Hz]');
    ylabel(c(it),'Frequency [Hz]');
end;
colormap jet;