
%% load invidual grid
p2d = '/home/rouxf/prj/TC/matFiles/';
dat = load([p2d,'individual_MNIwarped_grid.mat']);
Sgrid = dat.grid;
tm = dat.tm;

%% load invidual headmodel
load([p2d,'individual_headmodel.mat']);

%%
MNIcoord =[-14 -24 12;... % left and right thalamus MNI coordinates
            10 -24 12;...
           -18 -38 46;... % left and right parietal cortex MNI coordinates
            28 -28 46;...
           -19 -79 20;... % left and right occipital cortex MNI coordinates
            21 -83 19];
        
[CTFcoord,Gridcoord,VCidx] = convert_MNI2Source_gridCoord(tm,MNIcoord,Sgrid);

idx = find(Sgrid.inside ==1);
if Sgrid.pos(idx(VCidx),:)-Gridcoord ~= 0
    error('VC index does not lign up with output grid coordinates');
end;

%%
figure;
hold on;
ft_plot_vol(hdm,'facealpha',.5,'surfaceonly',1);
ft_plot_mesh(Sgrid.pos(Sgrid.inside,:),'vertexcolor','r','surfaceonly',1);
view(-44,30);

figure;
hold on;
ft_plot_vol(hdm,'facealpha',.5,'surfaceonly',1);
for it = 1:size(Gridcoord,1)
    ft_plot_dipole(Gridcoord(it,:), [-1 0 0],'unit',hdm.unit,'color','r');
end;
view(-44,30);
set(gcf,'Color','w');

figure;
hold on;
ft_plot_vol(hdm,'facealpha',.5,'surfaceonly',1);
for it = 1:size(Gridcoord,1)
    ft_plot_dipole(Gridcoord(it,:), [-1 0 0],'unit',hdm.unit,'color','r');
end;
view(88,40);
set(gcf,'Color','w');

figure;
hold on;
ft_plot_vol(hdm,'facealpha',.5,'surfaceonly',1);
for it = 1:size(Gridcoord,1)
    ft_plot_dipole(Gridcoord(it,:), [-1 0 0],'unit',hdm.unit,'color','r');
end;
view(149,20);
set(gcf,'Color','w');