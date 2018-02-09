function simulate_dipole_ROItest(param1,param2,param3)

if nargin ==0
    param1 = -pi/2;%-pi:pi/4:pi;%phase lag param
    param2 = 1;%0.2:.3:2.1;% SNR param
    param3 = 0;%;10e1; % number of random sources
end;

%%
p2d = '~rouxf/prj/TC/matFiles/';
[savepath] = p2d;

%% set the unit for geometrical data
SIunit = 'mm';

%% load single subject MEG data
load([savepath,'resting_state_recording_THA07_ICAcleanedMEGdataEO.mat']);
load([savepath,'resting_state_recording_THA07_ICAcleanedMEGdataEC.mat']);

load([savepath,'independent_alphaGenerators_timeSeries.mat']);
SigDat = dat;
for it = 1:size(r,1)
    r(it,it) = NaN;
end;

%%
ntrl = params.ntrl;%length(dc1.trial);% number of trials
trlt = max(params.t);%max(dc1.time{1});% duration of each trial in s
Fs = params.Fs;% sampling freq

%%
cfg = [];
cfg.layout = 'CTF275.lay';

[lay] = ft_prepare_layout(cfg);

meg_data1 = dc1; clear dc1;
meg_data2 = dc2; clear dc2;

meg_data1.grad = ft_convert_units(meg_data1.grad,SIunit);
meg_data2.grad = ft_convert_units(meg_data2.grad,SIunit);

meg_data1.fsample = Fs;
meg_data2.fsample = Fs;

%% load invidual grid
dat = load([p2d,'individual_MNIwarped_grid.mat']);
Sgrid = dat.grid;
tm = dat.tm;
load([p2d,'individual_headmodel.mat']);
load([p2d,'template_volumes.mat']);
load([p2d,'individual_MRIdata.mat']);

% dat = load([p2d,'template_volumes.mat']);
% Sgrid = dat.template_grid;
% tm = template_mri.transform;

%%
hdm = ft_convert_units(hdm,'mm');
Sgrid = ft_convert_units(Sgrid,'mm');
template_grid = ft_convert_units(template_grid,'mm');
mri = ft_convert_units(mri,'mm');
template_mri = ft_convert_units(template_mri,'mm');

%%
%pInf = '/home/rouxf/prj/TC/WFU_PickAtlas_3.0.5b/wfu_pickatlas/MNI_atlas_templates/';
%fn = 'atlas116.nii';

%pInf = '/home/rouxf/tbx/fieldtrip-20170618/template/atlas/afni/';
%fn = 'TTatlas+tlrc.HEAD';

%pInf = '/home/rouxf/tbx/fieldtrip-20170618/template/atlas/spm_anatomy/';
%fn = 'AllAreas_v17.hdr';

pInf = '/home/rouxf/tbx/fieldtrip-20170618/template/atlas/aal/';
fn = 'ROI_MNI_V4.nii';

%pInf = '/home/rouxf/tbx/fieldtrip-20170618/template/atlas/brainweb/';
%fn = 'brainweb_discrete.mat';

[atlas] = ft_read_atlas([pInf,fn]);

%%
cfg                     = [];
cfg.interpmethod        = 'nearest';
cfg.parameter           = 'tissue';

[maskedGrid] = ft_sourceinterpolate( cfg , atlas, template_grid );

%%
[ROIs] = {'Thalamus_L','Thalamus_R','Precuneus_L','Precuneus_R','Occipital_Inf_L','Occipital_Inf_R'};

[gridIdx] = find(Sgrid.inside);

selIx = {};
for it = 1:length(ROIs)
    selIx{it} = find(strcmp(atlas.tissuelabel,ROIs{it}));
end;

atlasIdx = cell(1,length(selIx));
GridCoord = cell(1,length(selIx));
VCidx = cell(1,length(selIx));
for it = 1:length( selIx )
    fprintf([num2str(it),'/',num2str(length( selIx ))]);
    atlasIdx{it} = find(maskedGrid.tissue == selIx{it});
    
    VCidx{it} = find(ismember(gridIdx,atlasIdx{it}));
    GridCoord{it} = Sgrid.pos(gridIdx(VCidx{it}),:);
    
    %GridCoord{it} = Sgrid.pos(atlasIdx{it},:);
    fprintf('\n');
end;

cnt = 0;sel = [];
for it = 1:length(GridCoord)
    
    if ~isempty(GridCoord{it})
        cnt = cnt+1;
        sel(cnt) = it;
    end;
end;

[GridCoord] = GridCoord(sel);
[atlasIdx]  = atlasIdx(sel);
[VCidx]     = VCidx(sel);
gC=GridCoord;
gClab = [];
idx = 1:length(gC{1});
for it = 1:length(gC)
    gClab(idx) = it*ones(1,length(idx));
    if it < length(gC)
        idx = idx(end)+1:idx(end)+length(gC{it+1});
    end;
end;
gClab = gClab';

a1 = [];a2 = [];a3 = [];
for it = 1:length( VCidx)
    a1 = [a1;GridCoord{it}]; 
    a2 = [a2;atlasIdx{it}]; 
    a3 = [a3;VCidx{it}];    
end;
[GridCoord] = a1;
[atlasIdx]  = a2;
[VCidx]     = a3;


%%
savePath = '/home/rouxf/prj/TC/matFiles/';
saveName = 'sourcePositions_ATLASbasedROI_';
for it = 1:length(ROIs)
    saveName = [saveName,ROIs{it},'_'];
end;

%% compute the coordinates in CTF sensor space for dipoles in the ROIs
CTFcoord     = zeros(length(VCidx),3);
GridCoord2   = zeros(length(VCidx),3);
VCidx2       = zeros(length(VCidx),1);



idx = find(Sgrid.inside ==1);
MNIcoord = template_grid.pos(idx(VCidx),:);
for it = 1:length(VCidx)
    fprintf([num2str(it),'/',num2str(length(VCidx))]);
           
        [CTFcoord(it,:),GridCoord2(it,:),VCidx2(it)] = convert_MNI2Source_gridCoord(tm,MNIcoord(it,:),Sgrid);
    
    fprintf('\n');
end;

for it = 1:length( VCidx )
    if any(Sgrid.pos(idx(VCidx(it)),:)-GridCoord2(it,:) ~= 0)
        error('VC index does not lign up with output grid coordinates');
    end;
end;

if ~isequal(GridCoord,GridCoord2)
    error('grid coordinates must match');
end;

%%
sanityCheck_dipolSimRoi;

%%
raw1 = cell(length(param1),length(param2));
raw2 = cell(length(param1),length(param2));

%%
for xt = 1:length(param1)
    for yt = 1:length(param2)
        
        %Orientation of the dipole [qx qy qz] in MNI corrdinates: ARS
        %note: x-axis goes towards R, y-axis goes towards anterior, z-axis
        %goes towards superior
        %[1 0 0] points to Nasion
        %[0 1 0] points to LPA
        %[0 0 1] points to vertex
        %[-1 0 0] points to Inion
        %[0 -1 0] points to RPA
        %[0 0 -1] points away from vertex towards the feet
        
        cfg = [];
        cfg.dipoleunit = 'nA*m';
        cfg.channelunit = 'T/m';
        cfg.headmodel = hdm;%forward model
        cfg.grad = meg_data2.grad;% gradiometer positions                        
        cfg.ntrials = ntrl;
        cfg.triallength = trlt;
        cfg.fsample = Fs;%meg_data.fsample;
        cfg.relnoise    = 0;
        
        cfg.dip.signal = cell(1,cfg.ntrials);
        
        %snr1 = zeros(cfg.ntrials,1);
        for it = 1:cfg.ntrials
            
            %%
            rand_mom = [];
            rand_pos = [];
                       
            idx = find(Sgrid.inside ==1);
            if Sgrid.pos(idx(VCidx),:)-GridCoord ~= 0
                error('VC index does not lign up with output grid coordinates');
            end;
            idx(VCidx) = [];
            
            mom = [-1 0 1];
            for zt = 1:param3
                idx = idx( randperm(length(idx)) );
                ix = idx(1);
                idx(1) = [];
                rand_pos = [rand_pos;Sgrid.pos(ix,:)];
                fl = zeros(1,3);
                for vt = 1:3
                    x = randperm(3);
                    fl(vt) = x(1);
                end;
                rand_mom = [rand_mom;mom(fl)'];
            end;
            
            %Position of dipole in cartesian space (x,y,z)
            cfg.dip.pos = [GridCoord;rand_pos];
            cfg.dip.mom = [];
            for gt = 1:size(GridCoord,1)
                if ismember(gClab(gt),[1 2])
                    for vt = 1:3
                        x = randperm(3);
                        fl(vt) = x(1);
                    end;
                    cfg.dip.mom = [cfg.dip.mom mom(fl)];
                else
                    cfg.dip.mom = [cfg.dip.mom -1 0 0];
                end;
            end;
            cfg.dip.mom = [cfg.dip.mom zeros(1,param3*3)];
            cfg.dip.mom = cfg.dip.mom';                       
            cfg.dip.mom(size(GridCoord,1)*3+1:end) = rand_mom;
            
            %[SNR] = compute_SNR(gsig,asig);
            
            if max(max(abs(r))) <.6
                f = 1;
                %snr1(it) = 10*log10(sqrt(mean(abs(sig2).^2))./sqrt(mean(abs(sig1).^2)));
            else
                error('Signal correlations must remain < 0.6');
            end;
            
            for gt = 1:length(gC)
                cfg.dip.signal{it} = [cfg.dip.signal{it};repmat(squeeze(SigDat(gt,it,:)),[1 length(gC{gt})])'];%  2 x thalamus + 6 x cortex
            end;
            
            for jt = 1:size(rand_pos,1)
                rand_sig = mean(mean(cfg.dip.signal{it}(3:4,:)))+1*mean(std(cfg.dip.signal{it}(3:4,:),0,2)).*randn(1,cfg.triallength*cfg.fsample+1);
                rand_sig = rand_sig./10;
                cfg.dip.signal{it} = [cfg.dip.signal{it};rand_sig];
            end;
            
        end;
        
        [raw1{xt,yt}] = ft_dipolesimulation(cfg);
        
        cfg = [];
        cfg.channel = {'MEG'};
        
        [raw1{xt,yt}] = ft_selectdata(cfg,raw1{xt,yt});
        
        cfg = [];
        cfg.channel = setdiff(meg_data2.label,setdiff(meg_data2.label,raw1{xt,yt}.label))
        
        [meg_data2] = ft_selectdata(cfg,meg_data2);
        
        %SNR1 = zeros(length(raw1{xt,yt}.trial),length(raw1{xt,yt}.label));
        for it = 1:length(raw1{xt,yt}.trial)
            
            m = mean(raw1{xt,yt}.trial{it},2)*ones(1,length(raw1{xt,yt}.time{it}));
            sd = std(raw1{xt,yt}.trial{it},0,2)*ones(1,length(raw1{xt,yt}.time{it}));
            chan_noise = 2*m+2*sd.*randn(length(raw1{xt,yt}.label),length(raw1{xt,yt}.time{1}));
            
            %[SNR1(it,:)] = compute_SNR(raw1{xt,yt}.trial{it},chan_noise);
            
            raw1{xt,yt}.trial{it} = raw1{xt,yt}.trial{it} + chan_noise;%
        end;
        
        %%
        cfg = [];
        cfg.dipoleunit = 'nA*m';
        cfg.channelunit = 'T/m';
        cfg.headmodel = hdm;%forward model
        cfg.grad = meg_data2.grad;% gradiometer positions                        
        cfg.ntrials = ntrl;
        cfg.triallength = trlt;
        cfg.fsample = Fs;%meg_data.fsample;
        %cfg.relnoise    = 0.15;
        cfg.relnoise    = 1.75;
        
        A = 1;
        cfg.dip.signal = cell(1,cfg.ntrials);
        %snr2 = zeros(cfg.ntrials,1);
        for it = 1:cfg.ntrials                        
            
            %%
            rand_mom = [];
            rand_pos = [];
            
            idx = find(Sgrid.inside ==1);
            if Sgrid.pos(idx(VCidx),:)-GridCoord ~= 0
                error('VC index does not lign up with output grid coordinates');
            end;
            idx(VCidx) = [];
            
            mom = [-1 0 1];
            for zt = 1:param3
                idx = idx( randperm(length(idx)) );
                ix = idx(1);
                idx(1) = [];
                rand_pos = [rand_pos;Sgrid.pos(ix,:)];
                fl = zeros(1,3);
                for vt = 1:3
                    x = randperm(3);
                    fl(vt) = x(1);
                end;
                rand_mom = [rand_mom;mom(fl)'];
            end;
            
            sig = zeros(size(GridCoord,1),cfg.triallength*cfg.fsample+1);
            for jt = 1:size(GridCoord,1)
                [sig(jt,:)] = A*(.5+(1.*randn(1,cfg.triallength*cfg.fsample+1)));
                if ismember(gClab(jt),[1 2 5 6])
                    sig(jt,:) = sig(jt,:)./10;
                else
                    sig(jt,:) = sig(jt,:)./100;
                end;
            end;            
            
            if isempty(gcp('NoCreate'))
                parpool(36,'SpmdEnabled',false);
            end;
            
%             r = zeros(size(sig,1),size(sig,1));
%             for jt = 1:size(sig,1)
%                 fprintf([num2str(jt),'/',num2str(size(sig,1))]);
%                 x1 = sig(jt,:)';
%                 dum = zeros(1,size(sig,1));
%                 parfor kt = 1:size(sig,1)
%                     dum(kt) = corr(x1,sig(kt,:)');
%                 end;
%                 r(jt,:) = dum;
%                 r(jt,jt) = 0;
%                 fprintf('\n');
%             end;
%             
%             if max(max(abs(r))) <.6
%                 %snr2(it) = 10*log10(sqrt(mean(abs(sig2).^2))./sqrt(mean(abs(sig1).^2)));
%             else
%                 error('Signal correlations must remain <0.6');
%             end;
            
            %Position of dipole in cartesian space (x,y,z)
            cfg.dip.pos = [GridCoord;rand_pos];
            cfg.dip.mom = [];
            for gt = 1:size(GridCoord,1)
                if ismember(gClab(gt),[1 2])
%                     for vt = 1:3
%                         x = randperm(3);
%                         fl(vt) = x(1);
%                     end;
%                     cfg.dip.mom = [cfg.dip.mom mom(fl)];
                    cfg.dip.mom = [cfg.dip.mom -1 0 0];
                else
                    cfg.dip.mom = [cfg.dip.mom -1 0 0];
                end;
            end;
            cfg.dip.mom = [cfg.dip.mom zeros(1,param3*3)];
            cfg.dip.mom = cfg.dip.mom';                       
            cfg.dip.mom(size(GridCoord,1)*3+1:end) = rand_mom;
                    
            cfg.dip.signal{it} = sig;
                        
            for jt = 1:size(rand_pos,1)
                rand_sig = mean(mean(sig(3:4,:)))+1*(mean(std(sig(3:4,:),0,2)).*randn(1,cfg.triallength*cfg.fsample+1));
                rand_sig = rand_sig./10;
                cfg.dip.signal{it} = [cfg.dip.signal{it};rand_sig];
            end;
            
        end;
        
        [raw2{xt,yt}] = ft_dipolesimulation(cfg);
        
        cfg = [];
        cfg.channel = meg_data2.label;
        
        [raw2{xt,yt}] = ft_selectdata(cfg,raw2{xt,yt});
        
        %SNR2 = zeros(length(raw2{xt,yt}.trial),length(raw2{xt,yt}.label));
        for it = 1:length(raw2{xt,yt}.trial)
            
            m = mean(raw2{xt,yt}.trial{it},2)*ones(1,length(raw2{xt,yt}.time{it}));
            sd = std(raw2{xt,yt}.trial{it},0,2)*ones(1,length(raw2{xt,yt}.time{it}));
            chan_noise = 2*m+2*sd.*randn(length(raw2{xt,yt}.label),length(raw2{xt,yt}.time{1}));
            
            %[SNR2(it,:)] = compute_SNR(raw2{xt,yt}.trial{it},chan_noise);
            
            raw2{xt,yt}.trial{it} = raw2{xt,yt}.trial{it} + chan_noise;%
        end;
    end;
end;

%%
clear sig* t;
for it = 1:length(param2)
    save_data1 = raw1(:,it);
    save_data2 = raw2(:,it);
    
    savename5 = ['2x_L_thalamic_and_4x_L_cortical_dipoles_alphagammaPAC_SNR',num2str(param2(it)),'_nRandSources',num2str(param3),'_simMEG_ROI.mat'];
    savename6 = ['2x_L_thalamic_and_4x_L_cortical_dipoles_nonrhythmic_SNR',num2str(param2(it)),'_nRandSources',num2str(param3),'_simMEG_ROI.mat'];
    
    save([savepath,savename5],'save_data1','CTFcoord','GridCoord','MNIcoord','VCidx','-v7.3');
    save([savepath,savename6],'save_data2','CTFcoord','GridCoord','MNIcoord','VCidx','-v7.3');
end;