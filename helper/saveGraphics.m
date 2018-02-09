%%
savepath = '/home/rouxf/prj/TC/figures/';
cnt = 3;
for it =1:2
    cnt = cnt+1;
    pos = get(gcf,'PaperPosition');
    dim = [diff([pos(4) pos(2)]) diff([pos(1) pos(3)])];
    set(gcf,'PaperPositionMode','auto');
    it
    print(it,'-dsvg','-opengl','-r800',[savepath,'sourcePACTC',num2str(cnt),'.svg']);
    print(it,'-dpdf','-opengl','-r800',[savepath,'sourcePACTC',num2str(cnt),'.pdf']);
    print(it,'-dtiff','-opengl','-r800',[savepath,'sourcePACTC',num2str(cnt),'.tiff']);
end;