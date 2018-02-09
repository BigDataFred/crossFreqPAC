%%
restoredefaultpath;
addpath('~froux/froux/fieldtrip-20151020/');
ft_defaults;
%%
p2df = '/bcbl/home/home_a-f/froux/project_reaction_times/matFiles/';
load([p2df,'THA07_V2.mri_MRI4MNI.mat']);
%%
cfg = [];
cfg.spmversion  = 'spm8';
cfg.template = '/bcbl/home/home_a-f/froux/spm8 2/templates/T1.nii';
cfg.coordinates = 'ctf';
cfg.nonlinear = 'no';

[norm] = ft_volumenormalise(cfg,mri);
%%
cfg = [];
cfg.viewmode = 'ortho';
cfg.location = [-44 8 43];%[-8 -60 31];%[10 -24 12];%[-14 -24 12];

ft_sourceplot(cfg,norm);
%%
[14 -27 12]
[-18 -27 12]
[-44 8 43]