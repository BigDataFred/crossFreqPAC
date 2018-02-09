function visualize_PAC_source_space(file1,file2)

%%
[ p2df ] = '~rouxf/prj/TC/matFiles/';
[ p2df2 ] = '~rouxf/prj/TC/matFilesRev/';

%%
if nargin ==0
    file1 = 'whole_brain_alphaGammaPAC_reversed_seed_based_virtual_channels_8:12Hz_1x_L_thalamic_and_2x_L_cortical_dipoles_alphagammaPAC_SNR1_nRandSources0_simMEG.mat';
    file2 = 'lcmv_spatial_filter_1:100Hz_1x_L_thalamic_and_2x_L_cortical_dipoles_alphagammaPAC_SNR1_nRandSources0_simMEG.mat';
end;

%% load anatomical volume data
load([p2df,'individual_MRIdata.mat']);

load([p2df,'template_volumes.mat'],'template_grid');

tmp = 'whole_brain_alphaGammaPAC_local_virtual_channels_\d{1,2}:\d{1,2}Hz_';
dum = file1;
% dum([regexp(file1,tmp):regexp(file1,'Hz_')+2]) = [];
% if length(dum)==length(file1)
%     tmp = 'whole_brain_alphaGammaPAC_seed_based_\d{1,6}_virtual_channels_\d{1,2}:\d{1,2}Hz_';
%     dum = file1;
% 	dum([regexp(file1,tmp):regexp(file1,'Hz_')+2]) = [];
% end;
% if length(dum)==length(file1)
%     tmp = 'whole_brain_alphaGammaPAC_reversed_seed_based_\d{1,6}_virtual_channels_\d{1,2}:\d{1,2}Hz_';
%     dum = file1;
% 	dum([regexp(file1,tmp):regexp(file1,'Hz_')+2]) = [];
% end;
% if length(dum)==length(file1)
%     tmp = 'whole_brain_alphaGammaPAC_reversed_seed_based_ortho_\d{1,6}_virtual_channels_\d{1,2}:\d{1,2}Hz_';
%     dum = file1;
% 	dum([regexp(file1,tmp):regexp(file1,'Hz_')+2]) = [];
% end;
% 
% if ~isempty(dum(regexp(dum,'alphagammaPAC'):regexp(dum,'alphagammaPAC')+12))
%     dum = [dum(1:regexp(dum,'alphagammaPAC')-1),'nonrhythmic' dum(regexp(dum,'alphagammaPAC')+13:length(dum))];
% end;

load([p2df2,dum],'VCidx');

%%
[files] = dir([p2df2,file1]);

%Y = zeros(length(files),length(VCidx));
%ph = zeros(length(files),1);

%%
for it = 1:size(files,2)
    
    load([p2df2,file2]);
    [lcmv] = ft_convert_units(lcmv,'mm');
    
    load([p2df2,files(it).name]);
    
    %% normalize individual T1 to template brain
    ft = dir('/home/rouxf/tbx/fieldtrip-*');
    [spm_template] = ['/home/rouxf/tbx/',ft.name,'/external/spm8/templates/T1.nii'];
    
    cfg = [];
    cfg.spmversion  = 'spm8';
    cfg.template = spm_template;
    cfg.coordinates = 'ctf';
    cfg.nonlinear = 'no';
    
    [norm] = ft_volumenormalise(cfg,mri);
    
    [norm] = ft_convert_units(norm,'mm');
    
    %%
    lcmv.avg.pow = zeros(size(lcmv.avg.pow));
    n = min(size(MI));
    
    for zt = 1:n
        
        if sign(size(MI,2)-1)==1
            lcmv.avg.pow(lcmv.inside) = MI(zt,:);
        else
            lcmv.avg.pow(lcmv.inside) = MI;
             %lcmv.avg.pow =  lcmv.avg.pow.*( lcmv.avg.pow > max( lcmv.avg.pow)*0.2);
        end;
        
        %%
%         for kt = 1:length(VCidx)
%             Y(it,kt) = MI(VCidx(kt));
%         end;
        
%         %%
%         lcmv.avg.pow(selIx) = 0;
%         lcmv.avg.pow(selIx(VCidx(3))) = 100;               
        
        %%
        cfg = [];
        cfg.parameter = 'avg.pow';
        
        [int] = ft_sourceinterpolate(cfg,lcmv,norm);
        
%         %%
%         cfg                 = [];
%         cfg.flipdim         = 'yes';
%         cfg.parameter       = 'pow';
%         cfg.filetype       = 'nifti';
%         %cfg.filename        = [p2df,'agc_surf_whole_brain';
%         %cfg.filename        = [p2df,'agc_surf_seed_region';
%         cfg.filename        = [p2df,'agc_surf_whole_brain'];
%         %cfg.scaling           = 'no';
%         cfg.coordsys        = 'spm';
%         cfg.datatype        = 'double';
%         cfg.crosshair       = 'no';
%         ft_volumewrite(cfg,int);
        
        %%
        cfg = [];
        cfg.method          = 'slice';
        cfg.funparameter    = 'pow';
        cfg.maskparameter   = cfg.funparameter;
        cfg.opacitymap      = 'rampup';
        cfg.funcolormap     = 'jet';
        cfg.colorbar        = 'no';
        cfg.nslices = 52;       
        cfg.slicerange    = [65 141];
        
        %cfg.funcolorlim = [max(lcmv.avg.pow)*0.5 max(lcmv.avg.pow)];
        
        ft_sourceplot(cfg,int);
        %%
%         %for bt = 1:length(selIx)
%         figure;
%         hold on;
%         x = mean(PAH(find(MI > max(MI)*0.9),:),1);
%         x = normalize_range(x);
%         bar(pbins,x);
%         y = sin(pbins);
%         y = normalize_range(y);
%         plot(pbins,y,'r');
%         y = cos(pbins);
%         y = normalize_range(y);
%         plot(pbins,y,'b');
%         axis tight;box on;
%         xlabel('\alpha-phase [rad]','Fontsize',14);
%         ylabel('\gamma-power [a.u.]','Fontsize',14);
%         set(gca,'XTick',[-pi -pi/2 0 pi/2 pi]);
%         set(gca,'XTickLabel',{'\pi' '-\pi/2' '0' '\pi/2' '\pi'});
%         set(gcf,'Color','w');
%         set(gca,'Fontsize',12);
%         set(gca,'LineWidth',3);
%         set(gca,'ticklength',get(gca,'ticklength')*2);
%         set(gca,'FontName','Arial');

        %end;
        %%
%         %%
%         [v,s_idx] = sort(lcmv.avg.pow(find(lcmv.inside)));
%         sel_idx = find(lcmv.inside==1);
%         
%         %MNIcoords = ft_warp_apply(inv(norm.transform),[lcmv.pos(sel_idx(s_idx(end-1)),:)],'homogenous')
%         
%         n = length(VCidx);
%         k = n/2;
%         
%         figure;                   
%         dim = size(PAH);
%         
%         if length(dim)==2
%             for kt = 1:length(VCidx)
%                 subplot(k,2,kt);
%                 hold on;
%                 x = squeeze(PAH(VCidx(kt),:));
%                 x = normalize_range(x);
%                 bar(pbins,x);
%                 y = sin(pbins);
%                 y = normalize_range(y);
%                 plot(pbins,y,'r');
%                 y = cos(pbins);
%                 y = normalize_range(y);
%                 plot(pbins,y,'b');
%                 axis tight;box on;
%                 xlabel('Alpha source phase [rad]');
%                 ylabel('Gammma source power [a.u.]');
%             end;
%             %[~,ix] = max(squeeze(PAH(VCidx(1),:)));
%             %[~,ix2] = max(squeeze(PAH(VCidx(2),:)));
%         else
%             plot(pbins,squeeze(PAH(zt,VCidx(1),:)),'r--','LineWidth',3);
%             plot(pbins,squeeze(PAH(zt,VCidx(2),:)),'b-o','LineWidth',3);
%             plot(pbins,squeeze(PAH(zt,VCidx(3),:)),'k-o','LineWidth',3);
%             [~,ix] = max(squeeze(PAH(zt,VCidx(1),:)));
%             [~,ix2] = max(squeeze(PAH(zt,VCidx(2),:)));
%         end;
%         
%         set(gcf,'Color','w');                
        
        %ph(it) = angle(exp(1i*(pbins(ix)-pbins(ix2))))*(180/pi);
        
%         %%
%         [i1,i2,i3] = ind2sub(size(int.pow),find(int.pow == max(max(max(int.pow)))));
%         MNIcoords = warp_apply(int.transform,[mean(i1) mean(i2) mean(i3)],'homogenous')
%         
%         idx = find(lcmv.inside==1);
%         idxL = find(sign(lcmv.pos(idx,1))==-1);
%         idxR = find(sign(lcmv.pos(idx,1))==1);
%         
%         [~,midx1] = max(lcmv.avg.pow(idx(idxL)));
%         Gridcoords1 = lcmv.pos(idx(idxL(midx1)),:)
%         
%         [~,midx2] = max(lcmv.avg.pow(idx(idxR)));
%         Gridcoords2 = lcmv.pos(idx(idxR(midx2)),:)
%         
%         figure;
%         subplot(222);
%         hold on;
%         if length(dim) == 2
%             plot(pbins,squeeze(PAH(idxL(midx1),:)),'b-o','LineWidth',3);
%             plot(pbins,squeeze(PAH(idxR(midx2),:)),'r--','LineWidth',3);
%         else
%             plot(pbins,squeeze(PAH(zt,idxL(midx1),:)),'b-o','LineWidth',3);
%             plot(pbins,squeeze(PAH(zt,idxR(midx2),:)),'r--','LineWidth',3);
%         end;
%         
%         axis tight;box off;
%         xlabel('Alpha source phase [rad]');
%         ylabel('Gammma source power [a.u.]');
%         set(gcf,'Color','w');
        
    end;
end;

% %%
% figure;
% subplot(221);
% hold on;
% h = [];
% h(1) = bar([1],[Y(1,1)],'FaceColor',[.9 0 0]);
% bar([2],[Y(2,1)],'FaceColor',[.9 0 0]);
% h(2) = bar([3],[Y(1,2)],'FaceColor',[0 0 .9]);
% bar([4],[Y(2,2)],'FaceColor',[0 0 .9]);
% 
% set(gca,'XTick',[1 2 3 4]);
% xlabel('Phase lag [deg]');
% set(gca,'XTickLabel',{'-90' '0' '-90' '0'});
% ylabel('Modulation index [au]');
% legend(h,'Parietal cortex','Thalamus');
% 
% subplot(223);
% hold on;
% bar([1],[Y(1,1)/Y(2,1)],'FaceColor',[.9 0 0]);
% bar([2],[Y(1,2)/Y(2,2)],'FaceColor',[0 0 .9]);
% xlim([-1 4]);
% 
% ylabel('Gain factor [au]');
% set(gca,'XTick',[]);
% set(gcf,'Color','white');
