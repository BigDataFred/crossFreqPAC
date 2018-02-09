%% set path defs
restoredefaultpath;
addpath('~froux/froux/fieldtrip-20151020/');
addpath('~froux/froux/project_reaction_times/mcode/helper/');    

ft_defaults;
%% load VC data
path2files = '/bcbl/home/home_a-f/froux/project_reaction_times/matFiles/';
[files] = dir([path2files,'virtual_channels_1:150Hz_1x_L_thalamic_and_2x_L_cortical_dipoles_alphagammaPAC_SNR1_sanity.mat']);
%%
ref = 3;

%%
for zt = 1:length(files)
    VC = load([path2files,files(zt).name]);
    % VC = load([path2files,'virtual_channels_1:150Hz_fixedPosAndMom.mat']);
    % for it = 1:length(VC.VC.trial)
    %   VC.VC.trial{it} = VC.VC.trial{it}';
    % end;
    %%
    [alpha] = VC.VC;
    
    clear VC;
    
    load('/bcbl/home/home_a-f/froux/project_reaction_times/matFiles/1x_L_thalamic_and_2x_L_cortical_dipoles_alphagammaPAC_SNR1_sanity.mat','VCidx');
    
    ref = VCidx(ref);
    %% hilbert transform the virtual channel data
    cfg = [];
    cfg.channel = alpha.label;%extract alpha phase
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [9 11];
    cfg.bpfilttype = 'fir';
    cfg.padding = max(alpha.time{1});%
    cfg.padtype = 'mirror';
    cfg.continuous = 'no';
    
    [alpha] = ft_preprocessing(cfg,alpha);% single trial alpha activity
    
        %%
    idx = setdiff(1:length(alpha.label),ref);%do orthogonalization for all except ref channel
    
    for nt = 1:length(alpha.trial) 
        
        signal = alpha.trial{nt};
        X = signal(ref,:);
        Y = signal(idx,:);
        
        parfor xt = 1:length(idx)       
            
            [Yo] = orthogonalize_time_domain( X , Y(xt,:) );
            
            if any(isnan(Yo))
                error('NaNs detected');
            end;
            Y(xt,:) = Yo;
        end;      
        
        alpha.trial{nt}(ref,:) = X;
        alpha.trial{nt}(idx,:) = Y;% replace the original with the orthogonalized part
        
    end;
    
    %%
    cfg = [];
    cfg.hilbert = 'angle';% phase in rad
    
    [alpha] = ft_preprocessing( cfg, alpha);
    
        %% concatenate
    [concat1] = zeros(length(alpha.label),length(alpha.time{1})*length(alpha.trial));
    
    idx = 1:length(alpha.time{1});
    for it = 1:length(alpha.trial)
        concat1(:,idx) = alpha.trial{it};
        idx = idx + length(alpha.time{it});
    end;
    clear alpha;
    
    %% start parallel pool
    if matlabpool('size') == 0
        matlabpool 64;%open;
    end;
    
    %% phase bins
    [pbins] = -pi:pi/8:pi;
    %preallocate
    [PLV] = zeros(size(concat1,1),1);
    
    %% loop over virtual channels
    [phi] = concat1;clear concat1;% extract the alpha phase from seed region
    
    phi1 = phi(ref,:);
    n = length(phi1);
    
    parfor it = 1:size(phi,1)%loop over virtual channels
        fprintf([num2str(it),'/',num2str(size(phi,1))]);
        
        phi2 = phi(it,:);
        
        PLV(it) = 1/n*abs(sum(exp(1i.*(phi1-phi2))));
        
        fprintf('\n');
    end;
    
    %% close parallel pool
    clear phi*;
    matlabpool close;  
    
    %% save data
    savepath = '/bcbl/home/home_a-f/froux/project_reaction_times/matFiles/';
    save([savepath,'whole_brain_alpha_PLV_seed_based_ortho_',files(zt).name],'PLV');
    
end;
%% clear and exit
clear;
exit;