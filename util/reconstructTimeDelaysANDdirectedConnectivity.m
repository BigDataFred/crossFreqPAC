%% set path defs
restoredefaultpath;

%fieldtrip
ft = dir('~rouxf/tbx/fieldtrip-*');
addpath(['~rouxf/tbx/',ft.name]);
ft_defaults;

%local scripts/functions
addpath(genpath('/home/rouxf/prj/TC/mcode/'));

p2dat = '/home/rouxf/prj/TC/matFiles/';
load([p2dat,'template_volumes.mat']);

%% normalize individual T1 to template brain
%load([p2df,'THA07_V2.mri_MRI4MNI.mat']);
[mri] = ft_read_mri(['~rouxf/prj/TC/THA07_V2.mri']);
ft = dir('/home/rouxf/tbx/fieldtrip-*');
[spm_template] = ['/home/rouxf/tbx/',ft.name,'/external/spm8/templates/T1.nii'];

cfg = [];      
cfg.template = spm_template;
%cfg.coordinates = 'ctf';
cfg.nonlinear = 'no';

[norm] = ft_volumenormalise(cfg,mri);

%%
for jt = 1
    
    %%    
    file1 = [ 'sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources1000_simMEG_ROI',num2str(jt),'_fieldMode1.mat' ];
    %load([p2dat,file1],'VCidx','gClab');

    %%
    [file2] = [ 'virtual_channels_1:100Hz_sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources1000_simMEG_ROI',num2str(jt),'_fieldMode1.mat' ];
    file2 = dir([p2dat,file2]);
    
    %load([p2dat,file2.name]);
    
    %%
    file3 = [ 'dics_PAC_sourceMap_alpha_sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources1000_simMEG_ROI',num2str(jt),'_fieldMode1.mat' ];
    file3 = dir([p2dat,file3]);
    %load([p2dat,file3.name]);
    file4 = [ 'dics_PAC_sourceMap_gamma_sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources1000_simMEG_ROI',num2str(jt),'_fieldMode1.mat' ];
    file4 = dir([p2dat,file4]);
    %load([p2dat,file4.name]);
    
    %%
    gIx = find(template_grid.inside==1);
    
    ref1 = source1{1}.avg.filter{gIx(VCidx(1))}';
    ref2 = source1{1}.avg.filter{gIx(VCidx(2))}';
    
    ix1 = setdiff(gIx,gIx(VCidx(1)));
    ix2 = setdiff(gIx,gIx(VCidx(2)));
    
    rxy1 = zeros(1,length(ix1));
    rxy2 = zeros(1,length(ix1));
    for it = 1:length(ix1)
        rxy1(it) = corr(ref1,source1{1}.avg.filter{ix1(it)}','Type','Spearman');
        rxy2(it) = corr(ref2,source1{1}.avg.filter{ix2(it)}','Type','Spearman');
    end;       
    
    c1 = template_grid.pos(gIx(VCidx(1)),:);
    d1 = zeros(1,length(ix1));
    for it = 1:length(ix1)
        c2 = template_grid.pos(ix1(it),:);
        d1(it) = sqrt(sum((c1-c2).^2));
    end;
    d1 = d1.*10;
    
    c1 = template_grid.pos(gIx(VCidx(2)),:);
    d2 = zeros(1,length(ix2));
    for it = 1:length(ix2)
        c2 = template_grid.pos(ix2(it),:);
        d2(it) = sqrt(sum((c1-c2).^2));
    end;
    d2 = d2.*10;
    
     %% autoregressive model
    cfg                     = [];
    cfg.channel             = VC.label(VCidx);
    
    [dum] = ft_selectdata( cfg , VC );
    
    %%
    C = zeros(2,2401);
    selIx = [1 5; 2 6];
    mIx = [];
    for kt = 1:size(selIx,1)
        for it = 1:length(dum.trial)
            sig1 = dum.trial{it}(selIx(kt,1),:);
            sig2 = dum.trial{it}(selIx(kt,2),:);                        
%             [dum1] = orthogonalize_time_domain(sig2 , sig1);
%             [dum2] = orthogonalize_time_domain(sig1 , sig2);  
%             dum.trial{it}(selIx(kt,1),:) = dum1;
%             dum.trial{it}(selIx(kt,2),:) = dum2;
            [c,lags] = xcorr(sig1,sig2,'coeff');
            C(kt,:) =C(kt,:)+c;
            [~,mIx(kt,it)] = max(c);
        end;
        C(kt,:) = C(kt,:)./it;
    end;
    
    %%
    cfg                     = [];
    cfg.resamplefs          = 200;
    cfg.detrend             = 'no';
    cfg.channel             = VC.label(VCidx([1 3]));
    
    [dumds] = ft_resampledata(cfg,dum);
    
    granger = cell(1,length(dumds.trial));
    for nt = 1:length(dumds.trial)
        
        cfg                     = [];
        cfg.trials              = nt;
        
        [dum2] = ft_selectdata(cfg,dumds);
        
        cfg                     = [];
        cfg.order               = 10;
        cfg.toolbox             = 'bsmart';
        cfg.zscore              = 'yes';
        
        [mdata] = ft_mvaranalysis(cfg,dum2);
        
        %% transfer function
        cfg                     = [];
        cfg.method              = 'mvar';
        
        [mfreq] = ft_freqanalysis(cfg,mdata);
        
        %%
        cfg                     = [];
        cfg.method              = 'granger';
        %cfg.granger.sfmethod   = 'bivariate';
        
        [granger{nt}] = ft_connectivityanalysis(cfg,mfreq);
    end;
end;

%%
v1 = []; v2 = []; v3 = []; v4 = [];
Gm1 = zeros(size(granger{1}.freq)); Gm2 = zeros(size(granger{1}.freq)); 
Gm3 = zeros(size(granger{1}.freq)); Gm4 = zeros(size(granger{1}.freq));
selIx = find(granger{1}.freq >= 9 & granger{1}.freq <=11);
for it = 1:length(granger)
    
    [v1(it)] = squeeze( mean(granger{it}.grangerspctrm(1,5,selIx),3));
    [v2(it)] = squeeze( mean(granger{it}.grangerspctrm(5,1,selIx),3));
    [v3(it)] = squeeze( mean(granger{it}.grangerspctrm(2,6,selIx),3));
    [v4(it)] = squeeze( mean(granger{it}.grangerspctrm(6,2,selIx),3));
    
    Gm1 = Gm1 + squeeze(granger{it}.grangerspctrm(1,5,:))';
    Gm2 = Gm2 + squeeze(granger{it}.grangerspctrm(5,1,:))';
    Gm3 = Gm3 + squeeze(granger{it}.grangerspctrm(2,6,:))';
    Gm4 = Gm4 + squeeze(granger{it}.grangerspctrm(6,2,:))';
    
end;
Gm1 = Gm1./it;
Gm2 = Gm2./it;
Gm3 = Gm3./it;
Gm4 = Gm4./it;

%%

figure
subplot(221);
a = gca;
hold on;
plot(granger{1}.freq,(Gm1+Gm3)./2,'LineWidth',3);
plot(granger{1}.freq,(Gm2+Gm4)./2,'r','LineWidth',3);
xlabel('Frequency [Hz]','Fontsize',14);
ylabel('Granger causality [a.u.]','Fontsize',14);
set(gca,'XTick',[10 50 90]);

subplot(222);
a = [a gca];
hold on;
for it = 1:length(v1);
    plot([1 2],[v1(it) v2(it)],'k');
end;
plot(ones(1,length(v1)),v1,'ko','MarkerFaceColor','b');
plot(2*ones(1,length(v1)),v2,'ko','MarkerFaceColor','r');

for it = 1:length(v1);
    plot([3 4],[v1(it) v2(it)],'k');
end;
plot(3*ones(1,length(v1)),v1,'ko','MarkerFaceColor','b');
plot(4*ones(1,length(v1)),v2,'ko','MarkerFaceColor','r');
axis tight;
xlim([0 5]);
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',{'Th' 'Ctx' 'Th' 'Ctx'});
ylabel('Granger causality [a.u.]','Fontsize',14);

subplot(223);
a = [a gca];
hold on;
[n1,x1] = hist(lags(mIx(:)),[-500:25:500]);
h = bar(x1,n1);
plot(median(lags(mIx(:)))*ones(1,2),[0 max(n1)],'r--','LineWidth',3);
set(h,'FaceColor',[.75 .75 .75],'EdgeColor','w');
axis tight;
xlabel('Time-lag [ms]','Fontsize',14);
ylabel('Count','Fontsize',14);

subplot(224);
a = [a gca];
[~,sIx1] = sort(d1);
[~,sIx2] = sort(d2);
rxy = (rxy1(sIx1)+rxy2(sIx2))./2;
d = (d1(sIx1)+d2(sIx2))./2;

hold on;
plot(d,rxy,'b.');
plot([d(1) d(end)],[.35 .35],'k--','LineWidth',3);
plot([d(1) d(end)],[-.35 -.35],'k--','LineWidth',3);
axis tight;ylim([-1 1]);

xlabel('Distance [mm]','Fontsize',14);
ylabel('Correlation Coefficient [a.u.]','Fontsize',14);

set(gcf,'Color','w');
set(a,'Fontsize',12);
set(a,'LineWidth',3);
set(a,'ticklength',get(gca,'ticklength')*2);
set(a,'FontName','Arial');

%%
savepath = '/home/rouxf/prj/TC/figuresrm/';
pos = get(gcf,'PaperPosition');
dim = [diff([pos(4) pos(2)]) diff([pos(1) pos(3)])];
set(gcf,'PaperPositionMode','auto');
if dim(2)>dim(1)
    set(gcf,'PaperOrientation','landscape');
end;

print(gcf,'-dsvg','-opengl','-r1000',[savepath,'granger_Causality',num2str(it),'.svg']);
print(gcf,'-dpdf','-opengl','-r1000',[savepath,'granger_Causality',num2str(it),'.pdf']);

 %%
 
 mask1 = rxy1 >  0.45;
 mask2 = rxy1 < -0.45;
 mask3 = rxy2 >  0.45;
 mask4 = rxy2 < -0.45;
 
 dum1 = template_grid;
 dum1.avg.pow = NaN(1,length(template_grid.inside));
 dum1.avg.pow(ix1) = abs((rxy1.*(mask1+mask2) + rxy2.*(mask3+mask4))./2);
 
 cfg = [];
 cfg.parameter = 'all';
 
 [int1] = ft_sourceinterpolate(cfg,dum1,norm);
 
%  
%  cfg                 = [];
%  cfg.flipdim         = 'yes';
%  cfg.parameter       = 'pow';
%  cfg.filetype       = 'nifti';
%  cfg.coordsys        = 'spm';
%  cfg.datatype        = 'double';
%  cfg.crosshair       = 'no';
%  
%  cfg.filename        = [ p2dat,'spatial_leak_LH'];
%  ft_volumewrite(cfg,int1);
%  cfg.filename        = [ p2dat,'spatial_leak_RH'];
%  ft_volumewrite(cfg,int2);
%  
 
 cfg = [];
 cfg.method          = 'ortho';
 cfg.downsample      = 2;
 cfg.funparameter    = 'pow';
 cfg.maskparameter   = cfg.funparameter;
 cfg.opacitymap      = 'rampup';
 cfg.opacitylim      = [0 1];
 cfg.funcolormap     = 'jet';
 cfg.colorbar        = 'no';
  cfg.crosshair          = 'no';
%  cfg.nslices = 52;
%  cfg.slicerange    = [65 141];
 cfg.location           = [0 -32 32];
 ft_sourceplot(cfg,int1);
 
 %%
figure;
colorbar;colormap jet;
caxis([0 1]);ca = caxis;
axis off;
cb= colorbar;
zlab = get(cb,'Ylabel');
set(zlab,'String','r_{xy} [a.u.]','Fontsize',14);
set(cb,'YTick',round([ca(1) ca(2)].*100)/100);
set(cb,'YTickLabel',round([ca(1) ca(2)].*100)/100,'Fontsize',14);
set(cb,'YLim',round([ca(1) ca(2)].*100)/100);
set(cb,'Fontsize',12);
set(cb,'LineWidth',3);
set(cb,'FontName','Arial');
set(cb,'ticklength',get(gca,'ticklength')*2);
set(gcf,'Color','w');

%%
savepath = '/home/rouxf/prj/TC/figures/';
pos = get(gcf,'PaperPosition');
dim = [diff([pos(4) pos(2)]) diff([pos(1) pos(3)])];
set(gcf,'PaperPositionMode','auto');
if dim(2)>dim(1)
    set(gcf,'PaperOrientation','landscape');
end;

print(gcf,'-dtiff','-opengl','-r600',[savepath,'spatial_leakage_thalamus_colorbar',num2str(1),'.tif']);
print(gcf,'-dsvg','-opengl','-r600',[savepath,'spatial_leakage_thalamus_colorbar',num2str(1),'.svg']);
print(gcf,'-dpdf','-opengl','-r600',[savepath,'spatial_leakage_thalamus_colorbar',num2str(1),'.pdf']);
