%%
restoredefaultpath;
addpath('~froux/froux/fieldtrip-20151020/');
ft_defaults;
addpath(genpath('~froux/froux/project_reaction_times/'));   
%%
SIunit = 'mm';

cfg = [];
cfg.dataset = '/bcbl/home/home_a-f/froux/project_reaction_times/THA07_DEVELOPMENT_20090924_12.ds';
cfg.continous = 'yes';
cfg.trl = [1 1*1200 0];
cfg.demean = 'yes';
cfg.detrend = 'yes';
cfg.channel = {'MEG' '-MLP12' '-MLT41' '-MRC12'  '-MRC14' '-MRC25' '-MRP56' '-MRT12' '-MRT21' '-MRT23' '-MLT57' '-MLT52' '-MRF22' '-MRF13' '-MRF24' '-MRF43' '-MLO22' '-MLF25' '-MRO31' '-MRO21'};

[meg_data] = ft_preprocessing(cfg);

meg_data.grad = ft_convert_units(meg_data.grad,SIunit);
%%
p2f = '~froux/froux/project_reaction_times/matFiles/';
load([p2f,'individual_headmodel.mat']);
load([p2f,'individual_MNIwarped_grid.mat']);
load([p2f,'individual_MRIdata.mat']);
load([p2f,'template_volumes.mat'],'template_mri');
%%
MNIcoord = [-30 -66 49;-9 -26 9];% left thalamus and left parietal cortex coordinates
[CTFcoord,Gridcoord,VCidx] = convert_MNI2Source_gridCoord(tm,MNIcoord,grid);

idx = find(grid.inside ==1);
if grid.pos(idx(VCidx),:)-Gridcoord ~= 0
    error('VC index does not lign up with output grid coordinates');
end;
%%
figure;
hold on;
ft_plot_vol(hdm,'facecolor','brain');
ft_plot_sens(meg_data.grad,'coil','yes','coildiameter',8);
%%
ft_determine_coordsys(mri,'interactive','false');hold on;
ft_plot_vol(hdm,'facealpha',.1);
ft_plot_mesh(grid.pos(grid.inside,:),'vertexcolor','g');
%%
ft_determine_coordsys(mri,'interactive','false');hold on;
ft_plot_vol(hdm,'facealpha',.1);
ft_plot_sens(meg_data.grad,'coil','yes','coildiameter',8);
%%
figure;
hold on;
%ft_determine_coordsys(norm,'interactive','false');
ft_plot_vol(hdm,'facecolor','brain','facealpha',.1);
ft_plot_sens(meg_data.grad,'coil','yes','coildiameter',8);

ft_plot_mesh(grid.pos(grid.inside,:),'vertexcolor','g');
%%
figure;
ft_plot_vol(hdm,'edgecolor','none','facealpha',.2);view([-60 30]);
hold on;
plot3(Gridcoord(1,1),Gridcoord(1,2),Gridcoord(1,3),'go','MarkerSize',12,'LineWidth',3);
if size(Gridcoord,1) > 1
    for nt = 1:size(Gridcoord,1)
    plot3(Gridcoord(nt,1),Gridcoord(nt,2),Gridcoord(nt,3),'go','MarkerSize',12,'LineWidth',3);
    end;
end;
ft_plot_mesh(grid.pos(grid.inside,:),'vertexcolor','k');
view([90 0])
%%
for nt = 1:size(MNIcoord,1)
    cfg = [];
    cfg.method = 'ortho';
    %cfg.location = MNIcoord(nt,:);
    cfg.interactive = 'no';
    
    ft_sourceplot(cfg,norm);
end;
%%
cfg = [];
cfg.parameter = 'anatomy';

[int] = ft_sourceinterpolate(cfg,norm,template_mri);
%% create figure handles
cfg = [];
cfg.location = MNIcoord;

loc = norm.transform\[cfg.location(:); 1];
loc = round(loc(1:3));

figure;
ft_plot_ortho(int.anatomy,'style','subplot','location',loc);
subplot(2,2,1);
hold on;
plot3(40,40,40,'g.','MarkerSize',18);
plot3(54,40,40,'g.','MarkerSize',18);
camlight;camlight;camlight;camlight;

subplot(2,2,2);
hold on;
plot3(40,51,40,'g.','MarkerSize',18);
plot3(54,51,40,'g.','MarkerSize',18);
camlight;camlight;camlight;camlight;

subplot(2,2,4);
hold on;
plot3(40,50,40,'g.','MarkerSize',18);
plot3(54,50,40,'g.','MarkerSize',18);
camlight;camlight;camlight;camlight;

set(gcf,'Color','white');
%%
cfg                     = [];
cfg.roi                 = [10 -24 12];
cfg.sphere              = 0.1;
cfg.round2nearestvoxel  = 'yes';

[mask] = ft_volumelookup(cfg, template_mri);
template_mri.mask = mask;
%%
cfg = [];
cfg.maskparameter = 'mask';

ft_sourceplot(cfg,template_mri);