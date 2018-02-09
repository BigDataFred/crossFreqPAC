function contrast_PACvsNoPAC_DICS(file1,file2)

%%
p2df = '~rouxf/prj/TC/matFiles/';
p2df2 = '~rouxf/prj/TC/matFilesRev/';

%%

load([p2df2,file1]);
load([p2df2,file2]);

raw1 = save_data1;
raw2 = save_data2;

%%
source1 = cell(length(raw1),1);
source2 = cell(length(raw2),1);

source3 = cell(length(raw1),1);
source4 = cell(length(raw2),1);

%% compute beamformer and source maps
for xt = 1:length(raw1)
    
    cfg = [];
    cfg.channel = {'MEG'};
    
    [dum1] = ft_selectdata(cfg,raw1{xt});
    [dum2] = ft_selectdata(cfg,raw2{xt});
    clear meg_data mri norm mni* loc r;
    
    %%
    grad1 = dum1.grad;
    grad2 = dum2.grad;
    
    %%
    dat = load([p2df,'individual_MNIwarped_grid.mat'],'grid','template_grid');
    source_grid = dat.grid;
    template_grid = dat.template_grid;
    
    load([p2df,'individual_headmodel.mat'],'hdm');
    
    Nx = length(template_grid.xgrid);
    Ny = length(template_grid.ygrid);
    Nz = length(template_grid.zgrid);
    
    %         %%
    %
    %         [common_dat] = ft_appenddata([],dum1,dum2);
    %
    %         cfg                     = [];
    %         cfg.method              = 'mtmfft';
    %         cfg.output              = 'powandcsd';
    %         cfg.foilim              = [10 10];
    %         cfg.taper               = 'dpss';
    %         cfg.tapsmofrq           = 1;
    %         cfg.pad                 = 'maxperlen';
    %         cfg.trials              = 'all';
    %
    %         [common_csd1] = ft_freqanalysis(cfg,common_dat);
    %
    %         cfg                     = [];
    %         cfg.method              = 'mtmfft';
    %         cfg.output              = 'powandcsd';
    %         cfg.foilim              = [70 70];
    %         cfg.taper               = 'dpss';
    %         cfg.tapsmofrq           = 10;
    %         cfg.pad                 = 'maxperlen';
    %         cfg.trials              = 'all';
    %
    %         [common_csd2] = ft_freqanalysis(cfg,common_dat);
    %         clear common_dat;
    %
    %         %%
    %         cfg                         = [];
    %         cfg.channel                 = {'MEG'};
    %         cfg.grid                    = source_grid;
    %         cfg.headmodel               = hdm;
    %         cfg.grad                    = grad1;
    %         cfg.method                  ='dics';
    %         cfg.grid.dim                =[Nx Ny Nz];
    %         %cfg.normalize              = 'yes';
    %         cfg.keepleadfield           = 'no';
    %
    %         cfg.dics.fixedori           = 'yes';
    %         %cfg.dics.lambda            = '1%'; %
    %         cfg.dics.powmethod          = 'trace';
    %         cfg.dics.projectnoise       = 'yes';
    %         cfg.dics.keepfilter         ='yes';% these filters are for computing virtual electrodes later
    %         cfg.dics.projectmom         = 'no';
    %         cfg.dics.keepmom            = 'no';
    %         cfg.dics.realfilter         = 'yes';
    %         cfg.dics.keepleadfield      = 'no';
    %
    %         cfg.frequency               = common_csd1.freq;
    %         [common_filt1] = ft_sourceanalysis(cfg,common_csd1);
    %         clear common_csd1;
    %
    %         cfg.frequency               = common_csd2.freq;
    %         [common_filt2] = ft_sourceanalysis(cfg,common_csd2);
    %         clear common_csd2;
    
    %%
    cfg                     = [];
    cfg.method              = 'mtmfft';
    cfg.output              = 'powandcsd';
    cfg.foilim              = [10 10];
    cfg.taper               = 'dpss';
    cfg.tapsmofrq           = 2;
    cfg.pad                 = 'maxperlen';
    cfg.trials              = 'all';
    
    [csd1] = ft_freqanalysis(cfg,dum1);
    [csd2] = ft_freqanalysis(cfg,dum2);
    
    cfg                     = [];
    cfg.method              = 'mtmfft';
    cfg.output              = 'powandcsd';
    cfg.foilim              = [70 70];
    cfg.taper               = 'dpss';
    cfg.tapsmofrq           = 15;
    cfg.pad                 = 'maxperlen';
    cfg.trials              = 'all';
    
    [csd3] = ft_freqanalysis(cfg,dum1);
    [csd4] = ft_freqanalysis(cfg,dum2);
    
    clear dum*;
    
    %%
    cfg                             = [];
    cfg.channel                     = {'MEG'};
    cfg.grid                        = source_grid;
    %cfg.grid.filter                 = common_filt1.avg.filter;  clear common_filt1;
    cfg.headmodel                   = hdm;
    cfg.method                      ='dics';
    cfg.grid.dim                    =[Nx Ny Nz];
    %cfg.normalize                  = 'yes';
    cfg.keepleadfield               = 'no';
    
    cfg.dics.fixedori               = 'yes';
    %cfg.dics.lambda                = '1%'; %
    cfg.dics.powmethod              = 'trace';
    cfg.dics.projectnoise           = 'yes';
    cfg.dics.keepfilter             ='yes';% these filters are for computing virtual electrodes later
    cfg.dics.projectmom             = 'no';
    cfg.dics.keepmom                = 'no';
    cfg.dics.realfilter             = 'yes';
    cfg.dics.keepleadfield          = 'no';
    cfg.grad                        = grad1;
    cfg.frequency                   = csd1.freq;
    
    [source1{xt}] = ft_sourceanalysis(cfg,csd1);
    source1{xt} = rmfield(source1{xt},'cfg');
    
    cfg.grad                        = grad2;
    cfg.frequency                   = csd2.freq;
    [source2{xt}] = ft_sourceanalysis(cfg,csd2);
    source2{xt} = rmfield(source2{xt},'cfg');
    
    %align the beamformer with the grid positions
    source1{xt}.pos = template_grid.pos;
    source1{xt}.dim = template_grid.dim;
    
    %align the beamformer with the grid positions
    source2{xt}.pos = template_grid.pos;
    source2{xt}.dim = template_grid.dim;
    
    cfg                             = [];
    cfg.channel                     = {'MEG'};
    cfg.grid                        = source_grid;
    %cfg.grid.filter                 = common_filt2.avg.filter;  clear common_filt1;
    cfg.headmodel                   = hdm;
    cfg.method                      ='dics';
    cfg.grid.dim                    =[Nx Ny Nz];
    %cfg.normalize                  = 'yes';
    cfg.keepleadfield               = 'no';
    
    cfg.dics.fixedori               = 'yes';
    %cfg.dics.lambda                = '1%'; %
    cfg.dics.powmethod              = 'trace';
    cfg.dics.projectnoise           = 'yes';
    cfg.dics.keepfilter             ='yes';% these filters are for computing virtual electrodes later
    cfg.dics.projectmom             = 'no';
    cfg.dics.keepmom                = 'no';
    cfg.dics.realfilter             = 'yes';
    cfg.dics.keepleadfield          = 'no';
    
    cfg.grad                        = grad1;
    cfg.frequency                   = csd3.freq;
    [source3{xt}] = ft_sourceanalysis(cfg,csd3);
    source3{xt} = rmfield(source3{xt},'cfg');
    
    cfg.grad                        = grad2;
    cfg.frequency                   = csd4.freq;
    [source4{xt}] = ft_sourceanalysis(cfg,csd4);
    source4{xt} = rmfield(source4{xt},'cfg');
    
    %align the beamformer with the grid positions
    source3{xt}.pos = template_grid.pos;
    source3{xt}.dim = template_grid.dim;
    
    %align the beamformer with the grid positions
    source4{xt}.pos = template_grid.pos;
    source4{xt}.dim = template_grid.dim;
    clear grad*;
    
end;

%% save output
savepath = p2df2;

savename1 = ['dics_PAC_sourceMap_alpha_',file1];
savename2 = ['dics_noPAC_sourceMap_alpha_',file2];
savename3 = ['dics_PAC_sourceMap_gamma_',file1];
savename4 = ['dics_noPAC_sourceMap_gamma_',file2];

save([savepath,savename1],'source1');
save([savepath,savename2],'source2');
save([savepath,savename3],'source3');
save([savepath,savename4],'source4');
