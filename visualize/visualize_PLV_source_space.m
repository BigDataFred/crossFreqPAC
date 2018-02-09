%%
restoredefaultpath;
addpath('~froux/froux/fieldtrip-20151020/');
addpath(genpath('~froux/froux/project_reaction_times/mcode/'));
ft_defaults;
%%
p2df = '/bcbl/home/home_a-f/froux/project_reaction_times/matFiles/';
%% load anatomical volume data
load([p2df,'individual_MRIdata.mat']);

%%
files = dir([p2df,'whole_brain_alpha_PLV_seed_based_ortho_virtual_channels_1:150Hz_1x_L_thalamic_and_2x_L_cortical_dipoles_alphagammaPAC_SNR1_sanity.mat']);
%%
for it = 1:length(files)
    load([p2df,'lcmv_spatial_filter_1:150Hz_1x_L_thalamic_and_2x_L_cortical_dipoles_alphagammaPAC_SNR1_sanity.mat']);
    
    [lcmv] = ft_convert_units(lcmv,'mm');
    
    load([p2df,files(it).name]);
    %% normalize individual T1 to template brain
    [spm_template] = '/bcbl/home/home_a-f/froux/fieldtrip-20151020/external/spm8/templates/T1.nii';
    
    cfg = [];
    cfg.spmversion  = 'spm8';
    cfg.template = spm_template;
    cfg.coordinates = 'ctf';
    cfg.nonlinear = 'no';
    
    [norm] = ft_volumenormalise(cfg,mri);
    
    [norm] = ft_convert_units(norm,'mm');
     
    %%
    lcmv.avg.pow = zeros(size(lcmv.avg.pow));
    lcmv.avg.pow(lcmv.inside) = PLV;
    %
    cfg = [];
    cfg.parameter = 'avg.pow';
    
    [int] = ft_sourceinterpolate(cfg,lcmv,norm);
    %%
    cfg                 = [];
    cfg.parameter       = 'pow';
    cffg.filetype       = 'nifti';
    cfg.filename        = '/bcbl/home/home_a-f/froux/project_reaction_times/matFiles/alpha_PLV_seed_region';
    
    ft_volumewrite(cfg,int);
    %%
    cfg = [];
    cfg.funparameter = 'pow';
    cfg.method = 'slice';
    cfg.nslices = 40;
    %cfg.funcolorlim = [.95 1];
    
    ft_sourceplot(cfg,int);
    
    set(gcf,'Color','w');    
    
end;
