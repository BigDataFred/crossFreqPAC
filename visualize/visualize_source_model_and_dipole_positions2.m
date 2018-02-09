savepath = '/home/rouxf/prj/TC/figures/figure1/';

%%
MNIcoord =[-10 -20 10;... % left and right thalamus MNI coordinates
    10 -20 10;...
    -24 -50 48;... % left and right parietal cortex MNI coordinates
    18 -46 50;...
    -24 -78 -4;... % left and right occipital cortex MNI coordinates
    16 -76 -4];

%% load invidual grid
p2d = '/home/rouxf/prj/TC/matFiles/';
dat = load([p2d,'individual_MNIwarped_grid.mat']);
Sgrid = dat.grid;
tm = dat.tm;

%% load invidual headmodel
load([p2d,'individual_headmodel.mat']);
load([p2d,'template_volumes.mat']);

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
gC2         = cell(1,length(VCidx)*2);
for it = 1:length(VCidx)    
    
    [mniC] = template_grid.pos(gridIdx(VCidx{it}),:);
    Lix = find(sign(mniC(:,1))==-1);
    Rix = find(sign(mniC(:,1))==1);

    cnt = cnt+1;
    VCidx2{cnt} = VCidx{it}(Lix);

    [~,mix]=min(sum(abs(template_grid.pos(gridIdx(VCidx2{cnt}),:)-repmat(MNIcoord(cnt,:),[length(VCidx2{cnt}) 1])),2).^2);    
    VCidx2{cnt} = VCidx2{cnt}(mix);
    Lix = Lix(mix);    
    GridCoord2{cnt} = GridCoord{it}(Lix,:);
    atlasIdx2{cnt} =  atlasIdx{it}(Lix);
    
    gC2{cnt} = GridCoord{it}(find(sign(mniC(:,1))==-1),:);
    
    cnt = cnt+1;
    VCidx2{cnt} = VCidx{it}(Rix);
    [~,mix]=min(sum(abs(template_grid.pos(gridIdx(VCidx2{cnt}),:)-repmat(MNIcoord(cnt,:),[length(VCidx2{cnt}) 1])),2).^2);    
    VCidx2{cnt} = VCidx2{cnt}(mix);
    Rix = Rix(mix);    
    GridCoord2{cnt} = GridCoord{it}(Rix,:);
    atlasIdx2{cnt} =  atlasIdx{it}(Rix);    
    
    gC2{cnt} = GridCoord{it}(find(sign(mniC(:,1))==1),:);
    
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

idx = find(Sgrid.inside ==1);
MNIcoord = template_grid.pos(idx(VCidx),:);

% %% visual check of geometrical objects
% figure;
% hold on;
% ft_plot_vol(hdm,'facealpha',.7,'surfaceonly',0);
% %ft_plot_mesh(Sgrid.pos(Sgrid.inside,:),'facecolor','cortex','surfaceonly',0);
% gI1 = gC2{1};
% gI2 = gC2{2};
% gI3 = gC2{3};
% gI4 = gC2{4};
% gI5 = gC2{5};
% gI6 = gC2{6};
% 
% for it = 1:size(gI1,1)
%     plot3(gI1(it,1),gI1(it,2),gI1(it,3),'ro','MarkerFaceColor','b');
% end;
% for it = 1:size(gI2,1)
%     gI2;
%     plot3(gI2(it,1),gI2(it,2),gI2(it,3),'o','Color',[.25 .25 .25],'MarkerFaceColor',[.25 .25 .25]);
% end;
% for it = 1:size(gI3,1)
%     plot3(gI3(it,1),gI3(it,2),gI3(it,3),'ro','MarkerFaceColor','r');
% end;
% for it = 1:size(gI4,1)
%     plot3(gI4(it,1),gI4(it,2),gI4(it,3),'ko','MarkerFaceColor','k');
% end;
% for it = 1:size(gI5,1)
%     plot3(gI5(it,1),gI5(it,2),gI5(it,3),'o','Color',[0 .5 .5],'MarkerFaceColor',[0 .5 .5]);
% end;
% for it = 1:size(gI6,1)
%     plot3(gI6(it,1),gI6(it,2),gI6(it,3),'o','Color',[.9 .5 .5],'MarkerFaceColor',[.9 .5 .5]);
% end;
% 
% view(-56,22);

%%
C = {'b' [.25 .25 .25] 'r' 'k'  [0 .5 .5] [.9 .5 .5]};

figure;
hold on;
ft_plot_vol(hdm,'facealpha',.5,'surfaceonly',1);
for it = 1:size(GridCoord,1)
    ft_plot_dipole(GridCoord(it,:), [-1 0 0],'unit',hdm.unit,'Color',C{it});
end;
view(-56,22);
set(gcf,'Color','w');

%%
for it = 8       
    pos = get(gcf,'PaperPosition');
    dim = [diff([pos(4) pos(2)]) diff([pos(1) pos(3)])];
    set(gcf,'PaperPositionMode','auto');
    if dim(2)>dim(1)
        set(gcf,'PaperOrientation','landscape');
    end;
    print(it,'-dsvg','-opengl','-r600',[savepath,'sourceModel_anatomy',num2str(it),'.svg']);
end;
