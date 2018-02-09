%%
ft = dir('/home/rouxf/tbx/fieldtrip-********');
addpath(['/home/rouxf/tbx/',ft.name]);
ft_defaults;

%%
p2f = '~rouxf/prj/TC/matFiles/';
load([p2f,'template_volumes.mat']);

%%
[template_grid] = ft_convert_units(template_grid,'mm');
[template_hdm] = ft_convert_units(template_hdm,'mm');
[template_mri] = ft_convert_units(template_mri,'mm');

%%
pInf = '/home/rouxf/prj/TC/WFU_PickAtlas_3.0.5b/wfu_pickatlas/MNI_atlas_templates/';
fn = 'atlas116.nii';

[atlasWFU] = ft_read_atlas([pInf,fn]);

%%
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';

maskedGrid = ft_sourceinterpolate( cfg , atlasWFU, template_grid );

%%
[ROIs] = {'Thalamus_L','Parietal_Sup_L','Frontal_Mid_L'};

[gridIdx] = find(template_grid.inside);

selIx = {};
selIx{1} = find(strcmp(maskedGrid.tissuelabel,ROIs{1}));
selIx{2} = find(strcmp(maskedGrid.tissuelabel,ROIs{2}));
selIx{3} = find(strcmp(maskedGrid.tissuelabel,ROIs{3}));

atlasIdx = cell(1,length(selIx));
GridCoord = cell(1,length(selIx));
for it = 1:length( selIx )
    
    atlasIdx{it} = find(maskedGrid.tissue == selIx{it});
    atlasIdx{it}(template_grid.inside(atlasIdx{it})==0) = [];
    atlasIdx{it} = find(ismember(gridIdx,atlasIdx{it}));
    
    GridCoord{it} = template_grid.pos(gridIdx(atlasIdx{it}),:);
    
end;

%%
dum = [];
dum.dim = template_grid.dim;
dum.inside = template_grid.inside;
dum.pos = template_grid.pos;
dum.method = 'average';
dum.avg.pow = NaN(length(dum.pos),1);
dum.avg.pow(gridIdx) = 0;
for it = 1:length(atlasIdx)
    dum.avg.pow(gridIdx([atlasIdx{it}])) = 1e6;
end;

%%
cfg = [];
cfg.parameter = 'avg.pow';

[int] = ft_sourceinterpolate(cfg,dum,template_mri);
        
%%
figure;
hold on;
ft_plot_vol(template_hdm,'facealpha',.7,'surfaceonly',0);
%ft_plot_mesh(template_grid.pos(template_grid.inside,:),'facecolor','cortex','surfaceonly',0);

for it = 1:size(GridCoord{1},1)
    plot3(GridCoord{1}(it,1),GridCoord{1}(it,2),GridCoord{1}(it,3),'ro','MarkerFaceColor','r');
end;
for it = 1:size(GridCoord{2},1)
    plot3(GridCoord{2}(it,1),GridCoord{2}(it,2),GridCoord{2}(it,3),'bo','MarkerFaceColor','b');
end;
for it = 1:size(GridCoord{3},1)
    plot3(GridCoord{3}(it,1),GridCoord{3}(it,2),GridCoord{3}(it,3),'go','MarkerFaceColor','g');
end;

view(-86,29);

%%
cfg = [];
cfg.funparameter = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.opacitymap = 'rampup';
cfg.method = 'slice';

ft_sourceplot( cfg , int );