%%
rand('state',sum(100*clock))

Fs = 1.2e3;
t = 1/Fs:1/Fs:5; % 1 min

mX = [];
MI = {};

[amp1] = 6;%;5:15;
[amp2] = 1.5;%:.1:1;%1:.5:2;
[amp3] = 3;%0.5:.5:5;
[amp4] = 1;%0.5:.5:5;

pfoi = 4:2:20;
afoi = 40:5:100;
pbins = -pi:pi/8:pi;
for it = 1:length(amp1)
    for jt = 1:length(amp2)
        for kt = 1:length(amp3)
            for lt = 1:length(amp4)
                
                a = ( amp1(it)*sin(2*pi*10.*t) );
                %a = a + std(a)*3*randn(1,length(t));
                g = ( amp2(jt)*sin(2*pi*60.*t)).*(sin(2*pi*10.*t)+1);
                %g = g + std(g)*3*randn(1,length(t));
                
                % figure;
                % subplot(121);
                % plot(t,a);
                % axis tight;
                % subplot(122);
                % plot(t,g);
                % axis tight;
                
                x = [];
                x.t = t;
                x.sig = (a+g)+(3*std(a+g)*randn(1,length(t)));
                %x.sig = x.sig -mean(x.sig);
                %x.sig = gradient(x.sig);
                
                [MI{it,jt,kt},~] = compute_MI(x,pfoi,afoi,pbins);
                mX(it,jt,kt) = max(max( MI{it,jt,kt} ));
            end;
        end;
    end;
end;
%figure;plot(mX(:),'bs-');
figure;pcolor(pfoi,afoi,squeeze(MI{end,end,end})');axis xy;colorbar;shading interp;lighting phong;