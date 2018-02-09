function simulate_dipole_ROI(param1,param2,param3,scf)

if nargin ==0
    param1 = -pi/2;%-pi:pi/4:pi;%phase lag param
    param2 = 1;%0.2:.3:2.1;% SNR param
    param3 = 1e1;%;10e1; % number of random sources
    scf = [30 3 30];
end;

%%
% MNIcoord =[-14 -24 12;... % left and right thalamus MNI coordinates
%     10 -24 12;...
%     -18 -38 46;... % left and right parietal cortex MNI coordinates
%     28 -28 46;...
%     -19 -79 20;... % left and right occipital cortex MNI coordinates
%     21 -83 19];

% MNIcoord =[-10 -20 10;... % left and right thalamus MNI coordinates
%     10 -20 10;...
%     -24 -50 48;... % left and right parietal cortex MNI coordinates
%     18 -46 50;...
%     -24 -78 -4;... % left and right occipital cortex MNI coordinates
%     16 -76 -4];

MNIcoord =[-14 -24 12;... % left and right thalamus MNI coordinates
    10 -24 12;...
    -18 -38 46;... % left and right parietal cortex MNI coordinates
    28 -28 46;...
    -19 -79 20;... % left and right occipital cortex MNI coordinates
    21 -83 19];

%%
p2d = '~rouxf/prj/TC/matFiles/';
p2d2 = '~rouxf/prj/TC/matFilesRev/';
[savepath] = p2d2;

%% set the unit for geometrical data
SIunit = 'mm';

%% load single subject real MEG data
load([p2d,'resting_state_recording_THA07_ICAcleanedMEGdataEO.mat']);%dc1
load([p2d,'resting_state_recording_THA07_ICAcleanedMEGdataEC.mat']);%dc2

%load([p2d2,'resting_state_recording_THA07_MEGdataEOnoCleaning.mat']);%dc1
%load([p2d2,'resting_state_recording_THA07_MEGdataECnoCleaning.mat']);%dc2

chck = dir([p2d2,'independent_alphaGenerators_timeSeries*.mat']);
chck = chck(end).name;
slab = chck(regexp(chck,'.mat')-1);

load([savepath,chck]);

if max(max(max(dat))) > 1
    error('amplitude factor is badly scaled');
end;

SigDat = dat;
SigDat(1:2,:,:) = SigDat(1:2,:,:).*scf(1);
SigDat(3:4,:,:) = SigDat(3:4,:,:).*scf(2);
SigDat(5:6,:,:) = SigDat(5:6,:,:).*scf(3);

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

pInf = '/home/rouxf/tbx/fieldtrip-20170618/template/atlas/afni/';
fn = 'TTatlas+tlrc.HEAD';

%pInf = '/home/rouxf/tbx/fieldtrip-20170618/template/atlas/spm_anatomy/';
%fn = 'AllAreas_v17.hdr';

%pInf = '/home/rouxf/tbx/fieldtrip-20170618/template/atlas/aal/';
%fn = 'ROI_MNI_V4.nii';

%pInf = '/home/rouxf/tbx/fieldtrip-20170618/template/atlas/brainweb/';
%fn = 'brainweb_discrete.mat';

[atlas] = ft_read_atlas([pInf,fn]);

%%
cfg                     = [];
cfg.interpmethod        = 'nearest';
%cfg.parameter           = 'tissue';
cfg.parameter           = 'brick1';

[maskedGrid] = ft_sourceinterpolate( cfg , atlas, template_grid );

%%
%[ROIs] = {'Thalamus_L','Thalamus_R','Precuneus_L','Precuneus_R','Occipital_Inf_L','Occipital_Inf_R'};
[ROIs] = {'Pulvinar','Brodmann area 7','Brodmann area 18'};
%[ROIs] = {'Lateral Geniculum Body','Brodmann area 7','Brodmann area 18'};

[gridIdx] = find(Sgrid.inside);

selIx = {};
for it = 1:length(ROIs)
    %selIx{it} = find(strcmp(atlas.tissuelabel,ROIs{it}));
    selIx{it} = find(strcmp(atlas.brick1label,ROIs{it}));
end;

atlasIdx = cell(1,length(selIx));
GridCoord = cell(1,length(selIx));
VCidx = cell(1,length(selIx));
for it = 1:length( selIx )
    fprintf([num2str(it),'/',num2str(length( selIx ))]);
    %atlasIdx{it} = find(maskedGrid.tissue == selIx{it});
    atlasIdx{it} = find(maskedGrid.brick1 == selIx{it});
    
    VCidx{it} = find(ismember(gridIdx,atlasIdx{it}));
%     if it == 2
%         VCidx{it} = VCidx{it}(find(template_grid.pos(gridIdx(VCidx{it}),3)>=43 & template_grid.pos(gridIdx(VCidx{it}),3)<=49));
%     elseif it ==3
%        VCidx{it} = VCidx{it}(find(template_grid.pos(gridIdx(VCidx{it}),3)>=17 & template_grid.pos(gridIdx(VCidx{it}),3)<=23));
%     end;        
    GridCoord{it} = Sgrid.pos(gridIdx(VCidx{it}),:);    
    fprintf('\n');
end;

cnt = 0;
VCidx2      = cell(1,length(VCidx)*2); 
GridCoord2  = cell(1,length(VCidx)*2); 
atlasIdx2   = cell(1,length(VCidx)*2); 
for it = 1:length(VCidx)    
    
    [mniC] = template_grid.pos(gridIdx(VCidx{it}),:);
    Lix = find(sign(mniC(:,1))==-1);
    Rix = find(sign(mniC(:,1))==1);

    cnt = cnt+1;
    VCidx2{cnt} = VCidx{it}(Lix);

    [~,mix]=min(sum(abs(template_grid.pos(gridIdx(VCidx2{cnt}),:)-repmat(MNIcoord(cnt,:),[length(VCidx2{cnt}) 1])),2).^2);    
    %VCidx2{cnt} = VCidx2{cnt}(mix);
    %Lix = Lix(mix);    
    GridCoord2{cnt} = GridCoord{it}(Lix,:);
    atlasIdx2{cnt} =  atlasIdx{it}(Lix);
    
    cnt = cnt+1;
    VCidx2{cnt} = VCidx{it}(Rix);
    [~,mix]=min(sum(abs(template_grid.pos(gridIdx(VCidx2{cnt}),:)-repmat(MNIcoord(cnt,:),[length(VCidx2{cnt}) 1])),2).^2);    
    %VCidx2{cnt} = VCidx2{cnt}(mix);
    %Rix = Rix(mix);    
    GridCoord2{cnt} = GridCoord{it}(Rix,:);
    atlasIdx2{cnt} =  atlasIdx{it}(Rix);    
    
end;
VCidx = VCidx2;
GridCoord = GridCoord2;
atlasIdx = atlasIdx2;


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
idx = 1:size(gC{1},1);
for it = 1:length(gC)
    gClab(idx) = it*ones(1,size(idx,1));
    if it < length(gC)
        idx = idx(end)+1:idx(end)+size(gC{it+1},1);
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

%% select x % of sources in ROI
[selIx] = cell(1,length( unique( gClab) ) );
for it = 1:length( unique( gClab ) )
    n = round(length(find(gClab==it))*0.3);
    selIx{it} = find(gClab == it);
    rIx = randperm(length(selIx{it}));
    selIx{it} = selIx{it}(rIx(1:n));
    gC{it} = gC{it}(rIx(1:n),:);
end;
selIx = [selIx{1};selIx{2};selIx{3};selIx{4};selIx{5};selIx{6}];

gClab = gClab( selIx );
GridCoord = GridCoord(selIx,:);
atlasIdx = atlasIdx(selIx,:);
VCidx = VCidx(selIx);

%%
savePath = '/home/rouxf/prj/TC/matFiles/';
saveName = 'sourcePositions_ATLASbasedROI_';
for it = 1:length(ROIs)
    saveName = [saveName,ROIs{it},'_'];
end;
saveName(regexp(saveName,' ')) = [];

% %% section added for revision to use Center of Gravity of ROIs
% % to go back to original line 260-261 must be uncommented!
% idx = find(Sgrid.inside ==1);
% MNIcoord = template_grid.pos(idx(VCidx),:);
% 
% cx = [];
% tmp = [];
% for it = 1:6
%     
%     A = MNIcoord(gClab == it,:);
%     dum = [mean(A(:,1)) mean(A(:,2)) mean(A(:,3))]
%     d = sqrt(sum((ones(size(template_grid.pos(idx(VCidx),:),1),1)*dum - template_grid.pos(idx(VCidx),:)).^2,2));
%    [~,mIx] = min(d);
%    tmp(it) = VCidx(mIx);
%    cx(it,:) = template_grid.pos(idx(VCidx(mIx)),:)
% end;
% MNIcoord = cx;
% VCidx = tmp;
% GridCoord = Sgrid.pos(idx(VCidx),:);

%% compute the coordinates in CTF sensor space for dipoles in the ROIs
CTFcoord     = zeros(length(VCidx),3);
GridCoord2   = zeros(length(VCidx),3);
VCidx2       = zeros(length(VCidx),1);

%uncomment this part if COG of ROIS is commented !!!
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
%sanityCheck_dipolSimRoi;

%%
raw1 = cell(length(param1),length(param2(1)));
raw2 = cell(length(param1),length(param2(1)));

%%
for xt = 1:length(param1)
    for yt = 1:length(param2(1))
        
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
        cfg.chanunit = {};
        for it = 1:length(meg_data2.grad.chanunit)
           meg_data2.grad.chanunit(it) = {'T'};
           cfg.chanunit(it) = {'T'};%'T';
        end;
        
        cfg.headmodel = hdm;%forward model
        cfg.grad = meg_data2.grad;% gradiometer positions     
        cfg.headmodel = ft_convert_units(cfg.headmodel,'m');
        cfg.grad = ft_convert_units(cfg.grad,'m');
        cfg.ntrials = ntrl;
        cfg.triallength = trlt;
        cfg.fsample = Fs;%meg_data.fsample;
        cfg.relnoise    = 0;
        
        cfg.dip.signal = cell(1,cfg.ntrials);
        cfg.dip.pos = cell(1,ntrl);
        cfg.dip.mom = cell(1,ntrl);
        
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
            cfg.dip.pos{it} = [GridCoord;rand_pos].*1e-3;% convert from mm to m
            cfg.dip.mom{it} = [];
            for gt = 1:size(GridCoord,1)
                if ismember(gClab(gt),[1 2])
                    if param2(1) == 2
                        fl = zeros(1,3);
                        for vt = 1:3
                            x = randperm(3);
                            fl(vt) = x(1);
                        end;
                        cfg.dip.mom{it} = [cfg.dip.mom{it} mom(fl)];
                    else
                        cfg.dip.mom{it} = [cfg.dip.mom{it} -1 0 0];
                    end;
                else
                    cfg.dip.mom{it} = [cfg.dip.mom{it} -1 0 0];
                end;
            end;
            cfg.dip.mom{it} = [cfg.dip.mom{it} zeros(1,param3*3)];
            cfg.dip.mom{it} = cfg.dip.mom{it}';                       
            cfg.dip.mom{it}(size(GridCoord,1)*3+1:end) = rand_mom;
            
            %[SNR] = compute_SNR(gsig,asig);
            
            if max(max(abs(r))) <.6
                f = 1;
                %snr1(it) = 10*log10(sqrt(mean(abs(sig2).^2))./sqrt(mean(abs(sig1).^2)));
            else
                error('Signal correlations must remain < 0.6');
            end;
            
            sig = zeros(length(VCidx),size(SigDat,3));
            for gt = 1:length(gC)
                fprintf([num2str(gt),'/',num2str( length(gC) )]);
                sig(gClab==gt,:) = repmat(squeeze(SigDat(gt,it,:)),[1 size(gC{gt},1)])';%  2 x thalamus + 6 x cortex
                fprintf('\n');
            end;
            cfg.dip.signal{it} = sig;
            
            M = mean(mean(cfg.dip.signal{it}(ismember(gClab,[3 4]),:)));
            SD = mean(std(cfg.dip.signal{it}(ismember(gClab,[3 4]),:),0,2));
            nsamp = cfg.triallength;
            Fs = cfg.fsample;
            sig = zeros( size(rand_pos,1),nsamp *Fs );
            parfor jt = 1:size(rand_pos,1)
                fprintf([num2str(jt),'/',num2str(size(rand_pos,1) )]);
                rand_sig =M+SD.*randn(1,nsamp*Fs);
                rand_sig = rand_sig./10;
                sig(jt,:)= rand_sig;
                fprintf('\n');
            end;
            cfg.dip.signal{it} = [ cfg.dip.signal{it};sig];
                        
        end;
                        
        [raw1{xt,yt}] = ft_dipolesimulation(cfg);
        
        cfg = [];
        cfg.channel = {'MEG'};
        
        [raw1{xt,yt}] = ft_selectdata(cfg,raw1{xt,yt});
        
        cfg = [];
        cfg.channel = setdiff(meg_data1.label,setdiff(meg_data1.label,raw1{xt,yt}.label))
        
        [meg_data1] = ft_selectdata(cfg,meg_data1);
        
        cfg                     = [];
        cfg.length              = 1.0008;
        
        [dum] = ft_redefinetrial( cfg, meg_data1 );
        
        for it = 1:length(raw1{xt,yt}.trial)            
            [raw1{xt,yt}.trial{it}] = raw1{xt,yt}.trial{it};% + dum.trial{it};                        
        end;
                
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
        cfg.chanunit = {};
        for it = 1:length(meg_data2.grad.chanunit)
           meg_data2.grad.chanunit(it) = {'T'};
           cfg.chanunit(it) = {'T'};%'T';
        end;
        
        cfg.headmodel = hdm;%forward model
        cfg.grad = meg_data2.grad;% gradiometer positions     
        cfg.headmodel = ft_convert_units(cfg.headmodel,'m');
        cfg.grad = ft_convert_units(cfg.grad,'m');
        cfg.ntrials = ntrl;
        cfg.triallength = trlt;
        cfg.fsample = Fs;%meg_data.fsample;
        %cfg.relnoise    = 0.15;
        cfg.relnoise    = 1.75;                
        
        cfg.dip.signal = cell(1,cfg.ntrials);
        cfg.dip.pos = cell(1,ntrl);
        cfg.dip.mom = cell(1,ntrl);
        
        A = repmat(scf,[2 1]);
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
            
            sig = zeros(size(GridCoord,1),cfg.triallength*cfg.fsample);
            for jt = 1:size(GridCoord,1)
                [sig(jt,:)] = A(gClab(jt))*(0+(1.*randn(1,cfg.triallength*cfg.fsample)));
            end;            
            
            if isempty(gcp('NoCreate'))
                parpool(36,'SpmdEnabled',false);
            end;            
            
            %Position of dipole in cartesian space (x,y,z)
            cfg.dip.pos{it} = [GridCoord;rand_pos].*1e-3;
            cfg.dip.mom{it} = [];
            for gt = 1:size(GridCoord,1)
                if ismember(gClab(gt),[1 2])
                    if param2(1) == 2
                        fl = zeros(1,3);
                        for vt = 1:3
                            x = randperm(3);
                            fl(vt) = x(1);
                        end;
                        cfg.dip.mom{it} = [cfg.dip.mom{it} mom(fl)];
                    else
                        cfg.dip.mom{it} = [cfg.dip.mom{it} -1 0 0];
                    end;
                else
                    cfg.dip.mom{it} = [cfg.dip.mom{it} -1 0 0];
                end;
            end;
            cfg.dip.mom{it} = [cfg.dip.mom{it} zeros(1,param3*3)];
            cfg.dip.mom{it} = cfg.dip.mom{it}';                       
            cfg.dip.mom{it}(size(GridCoord,1)*3+1:end) = rand_mom;
                    
            cfg.dip.signal{it} = sig;
                        
            M = mean(mean(cfg.dip.signal{it}(ismember(gClab,[3 4]),:)));
            SD = mean(std(cfg.dip.signal{it}(ismember(gClab,[3 4]),:),0,2));
            nsamp = cfg.triallength;
            Fs = cfg.fsample;
            sig = zeros( size(rand_pos,1),nsamp *Fs );
            parfor jt = 1:size(rand_pos,1)
                fprintf([num2str(jt),'/',num2str(size(rand_pos,1) )]);
                rand_sig =M+SD.*randn(1,nsamp*Fs);
                rand_sig = rand_sig./10;
                sig(jt,:)= rand_sig;
                fprintf('\n');
            end;
            cfg.dip.signal{it} = [ cfg.dip.signal{it};sig];
                        
        end;
        
        [raw2{xt,yt}] = ft_dipolesimulation(cfg);
        
        cfg = [];
        cfg.channel = meg_data1.label;
        
        [raw2{xt,yt}] = ft_selectdata(cfg,raw2{xt,yt});
        
        cfg                     = [];
        cfg.length              = 1.0008;
        
        [dum] = ft_redefinetrial( cfg, meg_data1 );
        
        for it = 1:length(raw1{xt,yt}.trial)
            [raw2{xt,yt}.trial{it}] = raw2{xt,yt}.trial{it};% + dum.trial{it};
        end;
        
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
for it = 1:length(param2(1))
    save_data1 = raw1(:,it);
    save_data2 = raw2(:,it);
    
    savename5 = [saveName,'dipoles_alphagammaPAC_SNR',num2str(param1),'_nRandSources',num2str(param3),'_simMEG_ROI_fieldMode',num2str(param2(1)),'.mat'];
    savename6 = [saveName,'dipoles_nonrhythmic_SNR',num2str(param1),'_nRandSources',num2str(param3),'_simMEG_ROI_fieldMode',num2str(param2(1)),'.mat'];
    
    save([savepath,savename5],'save_data1','CTFcoord','GridCoord','MNIcoord','VCidx','gClab','-v7.3');
    save([savepath,savename6],'save_data2','CTFcoord','GridCoord','MNIcoord','VCidx','gClab','-v7.3');
end;