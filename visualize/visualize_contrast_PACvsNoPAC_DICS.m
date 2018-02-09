function visualize_contrast_PACvsNoPAC_DICS(files1,files2,file3)

[p2df] = '~rouxf/prj/TC/matFiles/';
[p2df2] = '~rouxf/prj/TC/matFilesRev/';

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
simMode = 'simMEG';
SNR = 1;%[1 2];%sort([1 [0.2:0.3:2] ]);
nSources = 0;%10e1;
pbins = -pi/2;%:%-pi:pi/4:pi;
pbins = pbins*(180/pi);

%%
y  = zeros(length(SNR),length(pbins),2);
for nt = 1:length(SNR)%[1 length(SNR)]
    
    snr = num2str(SNR(nt));
    
    %%
    for it = 1:length(files1)
        load([p2df2,files1{it}]);
    end;
    
    for it = 1:length(files2)
        load([p2df2,files2{it}]);
    end;
    
    %%
    load([p2df2,file3],'VCidx');
    
    %%
    for it = 1:length(source1)
        dum1 = source1{it};
        dum1.avg.pow = (source1{it}.avg.pow-source2{it}.avg.pow)./(source2{it}.avg.pow);
        dum1.avg.pow = dum1.avg.pow.*(dum1.avg.pow > max(dum1.avg.pow)*0.1);
        idx = find(dum1.inside ==1);
        
        y(nt,it,1) = dum1.avg.pow(idx(VCidx(1)));
        y(nt,it,2) = dum1.avg.pow(idx(VCidx(2)));
        
        dum2 = source3{it};
        dum2.avg.pow =(source3{it}.avg.pow-source4{it}.avg.pow)./(source4{it}.avg.pow);
        dum2.avg.pow = dum2.avg.pow.*(dum2.avg.pow > max(dum2.avg.pow)*0.1);
        
        cfg = [];
        cfg.parameter = 'avg.pow';
        
        [int1] = ft_sourceinterpolate(cfg,dum1,norm);
        [int2] = ft_sourceinterpolate(cfg,dum2,norm);
        
        cfg = [];
        cfg.method          = 'ortho';
        cfg.funparameter    = 'pow';
        cfg.maskparameter   = cfg.funparameter;
        cfg.opacitymap      = 'rampup';
        cfg.funcolormap     = 'jet';
        cfg.colorbar        = 'no';
        cfg.crosshair          = 'no';
        %cfg.nslices         = 40;
        %cfg.funcolorlim     = [-5 5];
        cfg.downsample = 2;
        cfg.opacitylim = [max(dum1.avg.pow)*0.01 max(dum1.avg.pow)];
        ft_sourceplot(cfg,int1);
        
                
        %cfg.funcolorlim     = [-2.5 2.5];
        cfg.opacitylim = [max(dum2.avg.pow)*0.01 max(dum2.avg.pow)];
        ft_sourceplot(cfg,int2);        
        
%         %%
%         cfg                 = [];
%         cfg.parameter       = 'pow';
%         cffg.filetype       = 'nifti';
%         cfg.filename        = [p2df,'alpha_beamformer_surf'];
%         cfg.coordsys      = 'spm';
%         cfg.scaling       ='no';
%         
%         ft_volumewrite(cfg,int1);
%         
%         %%
%         cfg                 = [];
%         cfg.parameter       = 'pow';
%         cffg.filetype       = 'nifti';
%         cfg.filename        = [p2df,'gamma_beamformer_surf'];
%         cfg.coordsys      = 'spm';
%         cfg.scaling       ='no';
%         
%         ft_volumewrite(cfg,int2);
        
    end;
end;
%%

% figure;
% h = [];
% for it = 1:size(y,1)
%     subplot(3,3,it);
%     hold on;
%     h(1) = plot(pbins,log10(squeeze(y(it,:,1))),'bs-','LineWidth',3);
%     h(2) = plot(pbins,log10(squeeze(y(it,:,2))),'rs-','LineWidth',3);
%     ylim([0 2]);
%     set(gca,'XTick',[-180:90:180]);
%     title(['SNR:',num2str(SNR(it))],'FontWeight','bold');
%     xlabel('Phase [deg]');
%     ylabel('Alpha power [dB]');
% end;
% legend(h,'Left parietal cortex','Left thalamus');
% set(gcf,'Color','w');

% %%
% cfg = [];
% cfg.location = mni2ctf3;
%
% loc = int2.transform\[cfg.location(:); 1];
% loc = round(loc(1:3));
% %%
% ana = norm.anatomy;
% amin = min(ana(:));
% amax = max(ana(:));
% ana = (ana-amin)./(amax-amin);
%
% ori = eye(3);
% ori = ori(3,:);
%
% figure;
%
% k = 0;
% for it = 2:-1:0
%     k = k+1;
%     subplot(8,2,k);
%     a =  gca;
%     set(gcf,'currentaxes',a);
%     set(a,'DataAspectRatio',[1 1 1]);
%     hold on;
%     ft_plot_slice(ana,'location',[loc(1:2);loc(3)+it],'orientation',ori,'doscale',false,'transform',eye(4));%
%     ft_plot_slice(int2.pow,'location',[loc(1:2);loc(3)+it],'datmask',int2.pow,'opacitylim',[max(max(max(int2.pow)))/100*50 max(max(max(int2.pow)))],'style','flat','colormap','hot','colorlim',[3 4],'orientation',ori,'transform',eye(4));
%     axis tight;axis off;
%     set(a, 'view', [0 90]);
%     title(['+',num2str(it),'mm']);
%     camlight;camlight;camlight;camlight;
%
% end;
% for it = 1:13
%     k = k+1;
%     subplot(8,2,k);
%     a =  gca;
%     set(gcf,'currentaxes',a);
%     set(a,'DataAspectRatio',[1 1 1]);
%     hold on;
%     ft_plot_slice(ana,'location',[loc(1:2);loc(3)-it],'orientation',ori,'doscale',false,'transform',eye(4));%
%     ft_plot_slice(int2.pow,'location',[loc(1:2);loc(3)-it],'datmask',int2.pow,'opacitylim',[3 4],'style','flat','colormap','hot','colorlim',[3 4],'orientation',ori,'transform',eye(4));
%     axis tight;axis off;
%     set(a, 'view', [0 90]);
%     title(['-',num2str(it),'mm']);
%     camlight;camlight;camlight;camlight;
%
% end;
%
% set(gcf,'Color','w');