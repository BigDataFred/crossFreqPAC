%%
MNIcoord =[-14 -24 12;... % left and right thalamus MNI coordinates
    10 -24 12;...
    -18 -38 46;... % left and right parietal cortex MNI coordinates
    28 -28 46;...
    -19 -79 20;... % left and right occipital cortex MNI coordinates
    21 -83 19];

d = zeros(size(MNIcoord,1),size(MNIcoord,1)-1);
for it = 1:size(MNIcoord,1)
    x = MNIcoord(it,:);
    for jt = 1:size(MNIcoord,1)
        y = MNIcoord(jt,:);
        d(it,jt) = sqrt(sum((x-y).^2));
    end;
end;
d = d./10,; % convert mm to cm        

figure;
imagesc(d);
set(gca,'XTick',1:6);
set(gca,'YTick',1:6);
set(gca,'XTickLabel',{'Th-l','Th-r','BA7-l','BA7-r','BA18-l','BA18-r'});
set(gca,'YTickLabel',{'Th-l','Th-r','BA7-l','BA7-r','BA18-l','BA18-r'});
cb = colorbar;
zlab = get(cb,'YLabel');
set(zlab,'String','d [cm]');