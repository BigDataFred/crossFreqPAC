%%
Fs = 10e3;
t = 0:1/Fs:1;

phi = [-pi:pi/24:pi];


sig = sin(2*pi*2.*t);
sig2 = zeros(length(phi),length(sig));
r = zeros(length(phi),1);

figure;
for xt = 1:length(phi)
    
    sig2(xt,:) = real(hilbert(sig).*exp(1i*phi(xt)));
    r(xt) = corr(sig2(xt,:)',sig');
    
%     hold on;axis tight;
%     plot(t,sig)
%     plot(t,sig2(xt,:),'r');
%     title(['Phi:',num2str(phi(xt)*(180/pi)),' ,r:',num2str(round(r(xt)*100)/100)]);
%     pause(.45);
%     clf;
    
end;
%%
figure;
subplot(221);
hold on;
[~,idx] = min(abs(r));
idx(2) = length(r)-idx(1)+1;
plot(phi,r,'b-','LineWidth',3);
plot([phi(idx(1)) phi(idx(1))],[min(r) r(idx(1))],'k--');
plot([phi(idx(2)) phi(idx(2))],[min(r) r(idx(2))],'k--');
plot([min(phi) phi(idx(1))],[r(idx(1)) r(idx(1))],'k--');
plot([phi(idx(2)) max(phi)],[r(idx(2)) r(idx(2))],'k--');
set(gca,'XTick',[-pi:pi/2:pi]);
set(gca,'XTickLabel',[-pi:pi/2:pi].*180/pi);
axis tight;
xlabel('\alpha-phase lag [deg]','Fontsize',14);
ylabel('Correlation coefficient [au]','Fontsize',14);
set(gca,'LineWidth',2);
box on;
set(gca,'Fontsize',12);
set(gca,'LineWidth',3);
set(gca,'FontName','Arial');
set(gca,'ticklength',get(gca,'ticklength')*2);

subplot(222);
as = (.1/(2*pi));
X = phi.*(180/pi);
Y = phi.*as;
hold on;
plot(X,Y.*1000,'b-','LineWidth',3);
plot([X(idx(1)) X(idx(1))],[min(Y) Y(idx(1))].*1000,'k--');
plot([X(idx(2)) X(idx(2))],[min(Y) Y(idx(2))].*1000,'k--');
plot([min(X) X(idx(1))],[Y(idx(1)) Y(idx(1))].*1000,'k--');
plot([X(idx(2)) max(X)],[Y(idx(2)) Y(idx(2))].*1000,'k--');
set(gca,'XTick',[-pi:pi/2:pi].*180/pi);
set(gca,'YTick',[Y(idx(1)) 0 Y(idx(2))].*1e3);
%set(gca,'YTick',round(min(Y)*1000)/001000:as:round(max(Y)*1000)/1000);
axis tight;
xlabel('\alpha-phase lag [deg]','Fontsize',14);
ylabel('Time lag [ms]','Fontsize',14);
set(gca,'LineWidth',2);
box on;
set(gca,'Fontsize',12);
set(gca,'LineWidth',3);
set(gca,'FontName','Arial');
set(gca,'ticklength',get(gca,'ticklength')*2);

subplot(223);
hold on;
plot(Y.*1000,r,'b-','LineWidth',3);
plot([Y(idx(1)) Y(idx(1))].*1000,[min(r) r(idx(1))],'k--');
plot([Y(idx(2)) Y(idx(2))].*1000,[min(r) r(idx(2))],'k--');
plot([min(Y) Y(idx(1))].*1000,[r(idx(1)) r(idx(1))],'k--');
plot([Y(idx(2)) max(Y)].*1000,[r(idx(2)) r(idx(2))],'k--');
axis tight;
set(gca,'XTick',[Y(idx(1)) 0 Y(idx(2))].*1e3);
%set(gca,'XTick',round(min(Y)*10)/10:.125:round(max(Y)*10)/10);
xlabel('Time lag [ms]','Fontsize',14);
ylabel('Correlation coefficient [a.u.]','Fontsize',14);
set(gca,'LineWidth',2);
box on;
set(gca,'Fontsize',12);
set(gca,'LineWidth',3);
set(gca,'FontName','Arial');
set(gca,'ticklength',get(gca,'ticklength')*2);

subplot(224);
x = as+.025.*randn(1,1000);
[w,f] = hist(x);
bar(f.*1000,w);
axis tight;
xlabel('TC delay [ms]','Fontsize',14);
ylabel('Count','Fontsize',14);
set(gca,'LineWidth',2);
box on;
set(gca,'Fontsize',12);
set(gca,'LineWidth',3);
set(gca,'FontName','Arial');
set(gca,'ticklength',get(gca,'ticklength')*2);

set(gcf,'PaperPositionMode','auto');
set(gcf,'Color','w');

%%
savepath = '/home/rouxf/prj/TC/figures/';
savename = 'phaseOffsetANDampCorrelation';

for it = 1       
    pos = get(gcf,'PaperPosition');
    dim = [diff([pos(4) pos(2)]) diff([pos(1) pos(3)])];
    set(gcf,'PaperPositionMode','auto');
    if dim(2)>dim(1)
        set(gcf,'PaperOrientation','landscape');
    end;
    print(it,'-dsvg','-opengl','-r600',[savepath,savename,'.svg']);
    print(it,'-dpdf','-opengl','-r600',[savepath,savename,'.pdf']);
end;