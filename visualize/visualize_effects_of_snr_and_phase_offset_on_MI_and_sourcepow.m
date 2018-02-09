%%
restoredefaultpath;
addpath('/bcbl/home/home_a-f/froux/fieldtrip-20151020/');
addpath(genpath('~froux/froux/project_reaction_times/'));
ft_defaults;
savepath = '~froux/froux/project_reaction_times/matFiles/';
%%
p2d = '~froux/froux/project_reaction_times/matFiles/';
load([p2d,'individual_MNIwarped_grid.mat']);
%%
p2d = '~froux/froux/project_reaction_times/matFiles/';
load([p2d,'individual_headmodel.mat']);
%% compute the coordinates in CTF sensor space for a dipole in the
% left and right occipital cortex
%MNIcoord =[-30 -66 49];% left parietal cortex MNI coordinates
%MNIcoord =[30 -66 49];% right parietal cortex MNI coordinates
%MNIcoord =[-30 -66 49;31 -66 49];% left and right parietal cortex MNI coordinates
%MNIcoord = [-9 -26 9];% left thalamus coordinates
%MNIcoord = [-9 -26 9;9 -26 9];% left and right thalami coordinates

MNIcoord = [-30 -66 49;-9 -26 9];% left thalamus and left parietal cortex coordinates

%MNIcoord =[-30 -66 49;31 -66 49];% left and right parietal cortex MNI coordinates

[CTFcoord,Gridcoord,VCidx] = convert_MNI2Source_gridCoord(tm,MNIcoord,grid);

idx = find(grid.inside ==1);
if grid.pos(idx(VCidx),:)-Gridcoord ~= 0
    error('VC index does not lign up with output grid coordinates');
end;

[CTFcoord;Gridcoord]
%%
p2df = '/bcbl/home/home_a-f/froux/project_reaction_times/matFiles/';
%% load anatomical volume data
load([p2df,'individual_MRIdata.mat']);

load([p2df,'template_volumes.mat'],'template_grid');
%%
load([p2df,'lcmv_spatial_filter_1:150Hz_left_thalamus.mat']);

[lcmv] = ft_convert_units(lcmv,'mm');

%%
p2d = '/bcbl/home/home_a-f/froux/project_reaction_times/matFiles/';

files1 = dir([p2d,'dics_PAC_sourceMap_alpha_left_thalamic_and_left_parietal_dipoles_SNR*.mat']);
files2 = dir([p2d,'dics_PAC_sourceMap_gamma_left_thalamic_and_left_parietal_dipoles_SNR*.mat']);
%%
pow1 = zeros(length(files1),9,2);
mi1 = zeros(length(files1),2);
snr1 = zeros(length(files1),1);

pow2 = zeros(length(files2),9,2);
mi2 = zeros(length(files2),2);
snr2 = zeros(length(files1),1);
%%

for it =1:length(files1)
    
    dat = load([p2d,files1(it).name]);
    %%
    ix1 = regexp(files1(it).name,'SNR')+3;
    ix2 = min(regexp(files1(it).name(ix1:end),'.mat'))-1;
    ix2 = ix1+(ix2-1);
    snr1(it) = str2double(files1(it).name(ix1:ix2));
    %%        
    for kt = 1:length(dat.source1)
        ix = find(dat.source1{kt}.inside==1);
        pow1(it,kt,1) =  log10(dat.source1{kt}.avg.pow(ix(VCidx(1))));
        pow1(it,kt,2) =  log10(dat.source1{kt}.avg.pow(ix(VCidx(2))));
    end;
end;
%%
sel_idx = find(sign(snr1-2.1)==-1);
snr1 = snr1(sel_idx);
pow1 = pow1(sel_idx,:,:);
mi1 = mi1(sel_idx,:);
%%

for it =1:length(files2)
    
    dat = load([p2d,files2(it).name]);
    %%
    ix1 = regexp(files2(it).name,'SNR')+3;
    ix2 = min(regexp(files2(it).name(ix1:end),'.mat'))-1;
    ix2 = ix1+(ix2-1);
    snr2(it) = str2double(files2(it).name(ix1:ix2));
    %%        
    for kt = 1:length(dat.source3)
        ix = find(dat.source3{kt}.inside==1);
        pow2(it,kt,1) =  log10(dat.source3{kt}.avg.pow(ix(VCidx(1))));
        pow2(it,kt,2) =  log10(dat.source3{kt}.avg.pow(ix(VCidx(2))));
    end;
end;
%%
sel_idx = find(sign(snr2-2.1)==-1);
snr2 = snr2(sel_idx);
pow2 = pow2(sel_idx,:,:);
mi2 = mi2(sel_idx,:);
%%
c{1} = [.9 0 0];
c{2} = [.9 .9 0];
c{3} = [0 .9 0];
c{4} = [0 .9 .9];
c{5} = [.75 .75 .75];
c{6} = [.9 .75 0];
c{7} = [0 0 0];
c{8} = [0 0 .9];


pbins = (-pi:pi/4:pi).*(180/pi);
%%
figure;
subplot(221);
a = gca;
hold on;
for it = 1:size(pow1,1)
    plot(pbins,squeeze(pow1(it,:,1)),'s-','Color',c{it},'LineWidth',3);
end;
set(gca,'XTick',[-180:45:180]);
yl = get(gca,'YLim');

lgd = cell(length(snr1),1);
for it = 1:length(snr1)
    lgd(it) = {num2str(snr1(it))};
end;
legend(lgd);
title('Parietal cortex');


subplot(222);
a = [a gca];
hold on;
for it = 1:size(pow1,1)
    plot(pbins,squeeze(pow1(it,:,2)),'s-','Color',c{it},'LineWidth',3);
end;
set(gca,'XTick',[-180:45:180]);
yl = [yl get(gca,'YLim')];

lgd = cell(length(snr1),1);
for it = 1:length(snr1)
    lgd(it) = {num2str(snr1(it))};
end;
legend(lgd);
title('Thalamus');


%set(a,'YLim',[min(min(yl)) max(max(yl))]);
set(gcf,'Color','w');

for it = 1:length(a)
    xlabel(a(it),'Phase [deg]');
    ylabel(a(it),'Alpha power [log]');
end;

subplot(223);
a = gca;
hold on;
for it = 1:size(pow2,1)
    plot(pbins,squeeze(pow2(it,:,1)),'s-','Color',c{it},'LineWidth',3);
end;
set(gca,'XTick',[-180:45:180]);
yl = get(gca,'YLim');

lgd = cell(length(snr2),1);
for it = 1:length(snr2)
    lgd(it) = {num2str(snr2(it))};
end;
legend(lgd);
title('Parietal cortex');


subplot(224);
a = [a gca];
hold on;
for it = 1:size(pow2,1)
    plot(pbins,squeeze(pow2(it,:,2)),'s-','Color',c{it},'LineWidth',3);
end;
set(gca,'XTick',[-180:45:180]);
yl = [yl get(gca,'YLim')];

lgd = cell(length(snr2),1);
for it = 1:length(snr2)
    lgd(it) = {num2str(snr2(it))};
end;
legend(lgd);
title('Thalamus');


%set(a,'YLim',[min(min(yl)) max(max(yl))]);
set(gcf,'Color','w');

for it = 1:length(a)
    xlabel(a(it),'Phase [deg]');
    ylabel(a(it),'Gamma power [log]');
end;

