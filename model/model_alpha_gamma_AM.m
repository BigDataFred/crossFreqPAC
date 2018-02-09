%%
[phi,psi,x] = wavefun('coif5',10);

t = 1:length(psi);
t = t.*1/1000;
%%
sig = sin(2*pi*2.*t);
sig2 = -real(hilbert(sin(2*pi*0.2.*t)).*exp(1i*pi/3));
sig3 =sig.*-sig2;
sig4 =real(hilbert(sig3).*exp(1i*pi/0.85));
sig = phi.*sig;
%%
sig = sig(find(t>=9 & t<=14));
sig2 = sig2(find(t>=9 & t<=14));
sig3 = sig3(find(t>=9 & t<=14));
sig4 = sig4(find(t>=9 & t<=14));
t = t(find(t>=9&t<=14));

figure;
subplot(411);
plot([sig sig sig],'r','LineWidth',3);
axis off;axis tight;ylim([-1.5 1.5]);
subplot(412);
plot([sig2 sig2 sig2],'k','LineWidth',3);
axis off;axis tight;ylim([-1.5 1.5]);
subplot(413);
hold on;
plot([sig3(1:fix(length(sig4)/2)) sig3(fix(length(sig4)/2)+1:length(sig4)).*0.6...
   sig3(1:fix(length(sig4)/2)) sig3(fix(length(sig4)/2)+1:length(sig4)).*0.6...
   sig3(1:fix(length(sig4)/2)) sig3(fix(length(sig4)/2)+1:length(sig4)).*0.6],'r','LineWidth',3);
axis off;axis tight;ylim([-1.5 1.5]);
subplot(414);
hold on;
plot([sig2 sig2 sig2].*0.5,'k','LineWidth',3);
axis off;axis tight;ylim([-1.5 1.5]);

figure;
subplot(413);
hold on;
plot([sig3 sig3 sig3],'r','LineWidth',3);
axis off;axis tight;ylim([-1.5 1.5]);
subplot(414);
plot([real(hilbert(sig4(1:fix(length(sig4)/2))).*exp(1i*-pi/0.85)) real(hilbert(sig4(fix(length(sig4)/2+1):length(sig4))).*exp(1i*2*pi)) ...
    real(hilbert(sig4(1:fix(length(sig4)/2))).*exp(1i*-pi/0.85)) real(hilbert(sig4(fix(length(sig4)/2+1):length(sig4))).*exp(1i*0.5*pi)) ...
    real(hilbert(sig4(1:fix(length(sig4)/2))).*exp(1i*-pi/0.85)) real(hilbert(sig4(fix(length(sig4)/2+1):length(sig4))).*exp(1i*2*pi))],'k','LineWidth',3);
axis off;axis tight;ylim([-1.5 1.5]);