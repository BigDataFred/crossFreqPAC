%%
restoredefaultpath;
ft = dir(['~rouxf/tbx/fieldtrip-*']);
addpath(['~rouxf/tbx/',ft.name]);
ft_defaults;

savepath = '~rouxf/prj/TC/matFiles/'; 
p2d      = '~rouxf/prj/TC/matFiles/';

%%
load([savepath,'individual_MRIdata.mat'],'seg','mri');
load([savepath,'individual_headmodel.mat'],'hdm');

%%
SIunit = 'mm';
[spm_template] = ['/home/rouxf/tbx/',ft.name,'/external/spm8/templates/T1.nii'];

%%
%[mri] = ft_read_mri('~rouxf/prj/TC/THA07_V2.mri');

load([p2d,'template_volumesII.mat']);
%load([p2d,'THA07_V2.mri_MRI4MNI.mat']);

%% make sure all geo-objects have same units
% all template objects are in MNI-space
[template_mri] = ft_convert_units(template_mri,SIunit);
[template_seg] = ft_convert_units(template_seg,SIunit);
[template_hdm] = ft_convert_units(template_hdm,SIunit);
[template_grid] = ft_convert_units(template_grid,SIunit);

%[mri] = ft_convert_units(mri,SIunit); % individual  MRI (CTF head-space)

%% segment individual mri (head-space)
% cfg = [];
% cfg.spmversion = 'spm8';
% cfg.units = 'mm';
% cfg.coordsys = 'spm';
% cfg.template = spm_template; 
% 
% [seg] = ft_volumesegment(cfg, mri);
% 
% [seg] = ft_convert_units(seg,SIunit);

% %% make a headmodel (head-space)
% cfg = [];
% cfg.method = 'singleshell';
% cfg.tissue = 'brain';
% 
% [hdm] = ft_prepare_headmodel(cfg,seg);

%% align and warp individual MRI onto MNI template (MNI-space)
cfg = [];      
cfg.template = spm_template;
%cfg.coordinates = 'ctf';
cfg.nonlinear = 'no';

[norm] = ft_volumenormalise(cfg, mri);% this returns the warped invidiual MRI

%% mare an individual grid in head-space that is warped onto the template grid in MNI space

grid = template_grid;

tm = norm.cfg.final;% transformation matrix - from headspace to MNI space
tm = inv(tm);% reverse order of transform, ie from MNI space to ind. headspace

% reverse-normalise the template grid to individual subject
grid.pos = warp_apply(tm, template_grid.pos, 'homogenous'); %xyz-pos of grid are expressed in individual headspace, but each grid point
%has a corresponding xyz-point in the MNI grid 


grid = ft_convert_units(grid,SIunit);
norm = ft_convert_units(norm,SIunit);

%%
save([savepath,'individual_MNIwarped_gridII.mat'],'grid','norm','template_grid','tm');
%save([savepath,'individual_MRIdata.mat'],'seg','mri');
%save([savepath,'individual_headmodel.mat'],'hdm');






















