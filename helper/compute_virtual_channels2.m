function compute_virtual_channels2(bpfreq,fn)

%% default params
if nargin ==0
    bpfreq = [60 80];
    fn = '1x_L_thalamic_and_2x_L_cortical_dipoles_alphagammaPAC_SNR1_nRandSources10_simMEG.mat';
end;

%% path for read and write 
[p2df] = '~rouxf/prj/TC/matFiles/';
[p2df2] = '~rouxf/prj/TC/matFilesRev/';
[savepath] = p2df2;

%% load anatomical volume data
% template grid (MNI template)
load([p2df,'template_volumes.mat'],'template_grid');

%load the forward model: inward shift:0, gird res: 1cm
dat = load([p2df,'individual_MNIwarped_grid.mat']);
sgrid = dat.grid;
tm = dat.tm;

dat = load([p2df,'individual_headmodel.mat']);
hdm = dat.hdm;
clear dat;

%% get the individual grid dimensions
%to align grid with anatomical template
Nx = length(template_grid.xgrid);
Ny = length(template_grid.ygrid);
Nz = length(template_grid.zgrid);

%% load the MEG data 
load([p2df2,fn]);

%%
if exist('save_data1','var')
    raw1 = save_data1; clear save_data1;
elseif exist('dc2','var')
    raw1{1} = dc2; clear dc2;
end;

%%
for nt = 1:length(raw1)
    
    %% apply band-pass filter to simulated MEG-signals
    if ~isempty(bpfreq)
        cfg = [];
        cfg.channel = ft_channelselection('MEG',raw1{1}.grad.label);
        cfg.bpfilter = 'yes';
        cfg.bpfreq = bpfreq;% frequency range of interest
        cfg.bpfilttype = 'fir';
        cfg.padtype = 'mirror';
        cfg.padding = 2*max(raw1{nt}.time{1});% length of padding = 2x total length of recording
        
        [raw1{nt}] = ft_preprocessing(cfg,raw1{nt});
        
        % stopre the pre-proessing info for later usage
        raw1{nt}.bpfilter = cfg.bpfilter;
        raw1{nt}.bpfreq = cfg.bpfreq;
        raw1{nt}.bpfilttype = cfg.bpfilttype;
        raw1{nt}.padtype = cfg.padtype;
        raw1{nt}.padding = cfg.padding;
    end;
        
    %% compute the broad-band covariance matrix
    
    % concatenate the single trial data into one long trial
    concat = zeros(length(raw1{nt}.label),length(raw1{nt}.trial)*length(raw1{nt}.time{1}));
    idx = 1:length(raw1{nt}.trial{1});
    for it = 1:length(raw1{nt}.trial)
        concat(:,idx) = raw1{nt}.trial{it};
        idx = idx+length(raw1{nt}.trial{it});
    end;
    
    % create a dummy structure
    dum = struct;
    dum.label = raw1{nt}.label;% MEG-chan labels
    dum.fsample = raw1{nt}.fsample;
    dum.time{1} = 0:1/dum.fsample:(size(concat,2)-1)/dum.fsample;% fake time axis
    dum.trial{1} = concat;% concatenated signals
    
    % compute the MEG-chan covariance matrix
    cfg = [];
    cfg.channel = {'MEG'};
    cfg.covariance = 'yes';
    cfg.keeptrials = 'no';%average over trials
    
    [cov] = ft_timelockanalysis(cfg,dum);
    clear dum;
    
    %% build the broadband filter using LCMV beamformer
    
    % first filter is built from averaged covariance matrix
    cfg = [];
    cfg.grad = raw1{nt}.grad;
    cfg.grid = sgrid;
    cfg.headmodel = hdm;
    cfg.senstype = 'meg';
    cfg.method = 'lcmv';
    cfg.normalize = 'yes';
    
    %cfg.grid.pos = Gridcoord;
    %cfg.grid.mom = [-1 0 0];
    
    cfg.lcmv.fixedori = 'yes';% let dipole rotate freely = no (or not = yes)
    cfg.lcmv.keepfilter = 'yes';
    cfg.lcmv.projectmom = 'no';
    cfg.lcmv.projectnoise = 'yes';
    %cfg.lcmv.lambda = '10%';
    cfg.lcmv.lambda = '1%';
    
    cfg.grid.units = 'mm';
    cfg.grid.dim=[Nx Ny Nz];
    
    [lcmv] = ft_sourceanalysis(cfg,cov);
    clear cov;
    lcmv = rmfield(lcmv,'cfg');% get rid of large cfg-info
    
    % align beamformer with grid-space
    lcmv.pos = template_grid.pos;
    lcmv.dim = template_grid.dim;
    lcmv.nai = log(lcmv.avg.pow./lcmv.avg.noise);
    
    %% clear some data
    lcmv = rmfield(lcmv,'time');
    lcmv.avg = rmfield(lcmv.avg,'noise');
    if isfield(lcmv.avg,'mom')
        lcmv.avg = rmfield(lcmv.avg,'mom');
    end;       
    
    %% reconstruct the virtual channels by using the broadband filters
    
    % find those grid indices that are within the grid-space
    idx = find(lcmv.inside==1);
    
    % dummy structure for virtual channels
    VC= raw1{nt};
    VC.trial = {};
    VC.label = {};
    
    % get the number of trials
    N = length(raw1{nt}.trial);
    
    %preallocate and open parpool
    trial = zeros(length(idx),N,length(raw1{nt}.time{1}));
    
    [filt] = lcmv.avg.filter(idx);
    
    %% save data
    %save the lcmv-power-map
    if ~isempty(bpfreq)
        save([savepath,'lcmv_spatial_filter_',num2str(bpfreq(1)),':',num2str(bpfreq(2)),'Hz_',fn],'lcmv');
        %save([savepath,'lcmv_spatial_filter_',num2str(bpfreq(1)),':',num2str(bpfreq(2)),'Hz_fixedPosAndMom.mat'],'lcmv');
    else
        save([savepath,'lcmv_spatial_filter_broadband_',fn],'lcmv');
        %save([savepath,'lcmv_spatial_filter_broadband_fixedPosAndMom.mat'],'lcmv');
    end;
    clear lcmv;
    
    %%
    % build labels of virtual channels
    VC_label = cell(length(idx),1);
    for jt = 1:length(idx)
        VC_label(jt) = {['virtual_channel',num2str(idx(jt))]};
    end
    
    tN = length(raw1{nt}.time{1});
    Z = raw1{nt}.trial;
    clear raw1;
    
     %%
    if isempty(gcp('nocreate'))
        parpool(36,'SpmdEnabled', false);
    end;
    
    %%
    for it = 1:N% loop over trials
        
        fprintf([num2str(it),'/',num2str(N)]);
                
        % get the broadband MEG data
        sig = Z{it};
        
        %pre-allocate
        t = zeros(length(idx),tN);
        parfor jt = 1:length(idx)% loop over grid points
            
            
            [X] = (filt{jt}*sig)';% compute broadband virtual channel for each grid point
            
            [Y] = sorted_EVD(X);
            
            t(jt,:) = Y(:,1);
        end;
        trial(:,it,:) = t;
        clear t;
        
        fprintf('\n');
    end;
    clear Z;
    
    %%
    delete(gcp);
    
    %%
    % save the VC-labels
    VC.label = VC_label;
    clear VC_label;
    
    %%
    %re-organize VC data into ft-compatible trial structure
    for it = 1:size(trial,2)
        VC.trial{it} = squeeze(trial(:,it,:));
    end;
    clear trial;
    
    %% clear some data
    VC= rmfield(VC,'grad');
    VC= rmfield(VC,'cfg');    
    
    %% save the virtual channel data
    if ~isempty(bpfreq)
        save([savepath,'virtual_channels_',num2str(bpfreq(1)),':',num2str(bpfreq(2)),'Hz_',fn],'VC','-v7.3');
        %save([savepath,'virtual_channels_',num2str(bpfreq(1)),':',num2str(bpfreq(2)),'Hz_fixedPosAndMom.mat'],'VC','-v7.3');
    else
        save([savepath,'virtual_channels_broadband_',fn],'VC','-v7.3');
        %save([savepath,'virtual_channels_broadband_fixedPosAndMom.mat'],'VC','-v7.3');
    end;    
    
end;

%%
clear raw1;

%%
%exit;

