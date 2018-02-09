restoredefaultpath;
ft = dir('/home/rouxf/tbx/fieldtrip-*');
addpath(['/home/rouxf/tbx/',ft.name]);           
ft_defaults;
savepath = '~rouxf/prj/TC/matFiles/'; 

%%
outname = ([savepath,'MNI_MRI_template.mat']);

template = ['/home/rouxf/tbx/',ft.name,'/external/spm8/templates/T1.nii'];
[template_mri] = ft_read_mri(template);
template_mri.coordsys = 'spm';
save(outname,'template_mri');

load(outname,'template_mri');

[template_mri] = ft_convert_units(template_mri,'cm');

%% segment the template brain and construct a volume conduction model (headmodel)
outname =[savepath,'MNI_MRI_template_segmented.mat'];

cfg = [];
cfg.spmversion = 'spm8';
cfg.units = 'mm';
cfg.coordsys = 'spm';
cfg.template = template; 
[template_seg] = ft_volumesegment(cfg, template_mri);
save(outname,'template_seg');

load([savepath,'MNI_MRI_template_segmented.mat'],'template_seg');
template_seg.anatomy = template_mri.anatomy;

[template_seg] = ft_convert_units(template_seg,'cm');

%% build headmodel (volume conduction model)
cfg = [];
cfg.method = 'singleshell';

[template_hdm] = ft_prepare_headmodel(cfg, template_seg);
template_hdm = ft_convert_units(template_hdm,'cm');

%% this dummy gradiometer definition is required for prepare_dipole_grid and headmodelplot
cfg = [];
cfg.layout = 'CTF275.lay';

[lay] = ft_prepare_layout(cfg);

template_grad = [];
template_grad.pnt = [];
template_grad.ori = [];
template_grad.tra = [];
template_grad.label = {};
template_grad.unit = 'cm';

%% construct the dipole grid in the template brain coordinates
% % the source units are in cm
% % the negative inwardshift means an outward shift of the brain surface for
% % inside/outside detection

gridres = 0.1;% resolution in cm
inward_shift =0;

cfg = [];
cfg.grid.xgrid = -20:gridres:20;
cfg.grid.ygrid = -20:gridres:20;
cfg.grid.zgrid = -10:gridres:20;
cfg.grid.tight = 'yes';
cfg.inwardshift = inward_shift;
cfg.headmodel = template_hdm;
cfg.grad = template_grad;

[template_grid] = ft_prepare_sourcemodel(cfg);

[template_grid] = ft_convert_units(template_grid,'cm');

template_grid.res = cfg.grid.xgrid(2)-cfg.grid.xgrid(1); 
template_grid.inwardshift = inward_shift;

%%
%make a figure of the single subject headmodel, sensor positions and
%grid positions
ft_determine_coordsys(template_mri,'interactive','false');hold on;
ft_plot_vol(template_hdm);

ft_determine_coordsys(template_mri,'interactive','false');hold on;
ft_plot_mesh(template_grid.pos(template_grid.inside,:),'vertexcolor','r');

figure;
hold on;
ft_plot_vol(template_hdm,'facecolor','brain','facealpha',.1);
ft_plot_mesh(template_grid.pos(template_grid.inside,:),'vertexcolor','r');

%%
save([savepath,'template_volumesII.mat'],'template_hdm','template_grid','template_seg','template_mri');
