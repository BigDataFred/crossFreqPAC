%% visual check of geometrical objects
figure;
hold on;
ft_plot_vol(hdm,'facealpha',.7,'surfaceonly',0);
%ft_plot_mesh(Sgrid.pos(Sgrid.inside,:),'facecolor','cortex','surfaceonly',0);
gI1 = gC{1};
gI2 = gC{2};
gI3 = gC{3};
gI4 = gC{4};
gI5 = gC{5};
gI6 = gC{6};

for it = 1:size(gI1,1)
    plot3(gI1(it,1),gI1(it,2),gI1(it,3),'ro','MarkerFaceColor','r');
end;
for it = 1:size(gI2,1)
    gI2;
    plot3(gI2(it,1),gI2(it,2),gI2(it,3),'bo','MarkerFaceColor','b');
end;
for it = 1:size(gI3,1)
    plot3(gI3(it,1),gI3(it,2),gI3(it,3),'go','MarkerFaceColor','g');
end;
for it = 1:size(gI4,1)
    plot3(gI4(it,1),gI4(it,2),gI4(it,3),'co','MarkerFaceColor','c');
end;
for it = 1:size(gI5,1)
    plot3(gI5(it,1),gI5(it,2),gI5(it,3),'ko','MarkerFaceColor','k');
end;
for it = 1:size(gI6,1)
    plot3(gI6(it,1),gI6(it,2),gI6(it,3),'yo','MarkerFaceColor','y');
end;

view(-86,29);

%%
dum = [];
dum.dim = template_grid.dim;
dum.inside = template_grid.inside;
dum.pos = template_grid.pos;
dum.method = 'average';
dum.avg.pow = NaN(length(dum.pos),1);
dum.avg.pow(gridIdx) = 0;
dum.avg.pow(gridIdx(VCidx)) = 1;


cfg                     = [];
cfg.filename            = 'ROIsanityCheckBA18';
cfg.parameter           = 'pow';
cfg.filetype            = 'nifti';

ft_volumewrite( cfg , dum );

%% normalize individual T1 to template brain
ft = dir('/home/rouxf/tbx/fieldtrip-********');
[spm_template] = ['/home/rouxf/tbx/',ft.name,'/external/spm8/templates/T1.nii'];

cfg = [];
cfg.spmversion  = 'spm8';
cfg.template = spm_template;
cfg.coordinates = 'ctf';
cfg.nonlinear = 'no';

[norm] = ft_volumenormalise(cfg,mri);

[norm] = ft_convert_units(norm,'mm');

%%
cfg = [];
cfg.parameter = 'avg.pow';

[int] = ft_sourceinterpolate(cfg,dum,norm);

%%
cfg = [];
cfg.funparameter = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.opacitymap = 'rampup';
cfg.method = 'slice';
cfg.nslices         = 40;

ft_sourceplot( cfg , int );