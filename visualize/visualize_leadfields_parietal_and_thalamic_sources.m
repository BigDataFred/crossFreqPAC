%%
restoredefaultpath;
addpath('/bcbl/home/home_a-f/froux/fieldtrip-20151020/');
addpath(genpath('~froux/froux/project_reaction_times/mcode/'));
ft_defaults;
savepath = '~froux/froux/project_reaction_times/matFiles/';
%%
SIunit = 'mm';
spm_template = '/bcbl/home/home_a-f/froux/fieldtrip-20151020/external/spm8/templates/T1.nii';
%%
cfg = [];
cfg.dataset = '~froux/froux/project_reaction_times/THA07_DEVELOPMENT_20090924_12.ds';
cfg.continous = 'yes';
cfg.trl = [240*1200+1 480*1200 0];
cfg.demean = 'yes';
cfg.detrend = 'yes';
cfg.channel = {'MEG' '-MLP12' '-MLT41' '-MRC12'  '-MRC14' '-MRC25' '-MRP56' '-MRT12' '-MRT21' '-MRT23' '-MLT57' '-MLT52' '-MRF22' '-MRF13' '-MRF24' '-MRF43' '-MLO22' '-MLF25' '-MRO31' '-MRO21'};

[meg_data] = ft_preprocessing(cfg);

meg_data.grad = ft_convert_units(meg_data.grad,SIunit);
%%
cfg = [];
cfg.channel = {'MEG'};
cfg.length = 7;

[meg_data] = ft_redefinetrial(cfg,meg_data);
%%
p2d = '~froux/froux/project_reaction_times/matFiles/';
load([p2d,'individual_MNIwarped_grid.mat']);
ix = find(grid.inside==1);
%%
p2d = '~froux/froux/project_reaction_times/matFiles/';
load([p2d,'individual_headmodel.mat']);
%%

load('/bcbl/home/home_a-f/froux/project_reaction_times/matFiles/template_lead_field.mat');

figure;

subplot(321);
plot(lead_field{ix(VCidx(1))},'k');
axis tight;
aX = get(gca,'YLim');
aX = aX.*2;
ylim(aX);
title('Left parietal cortex');

subplot(322);
hold on;
plot(lead_field{ix(VCidx(2))},'k');
%plot(lead_field{ix(VCidx(2))}.*0.8,'Color',[.9 .9 0]);
%plot(lead_field{ix(VCidx(2))}.*1.1,'Color',[.9 .75 0]);
plot(lead_field{ix(VCidx(2))}.*15,'Color',[.9 0 0]);
axis tight;
ylim(aX);
title('Left thalamus');



cfg = [];
cfg.layout = 'CTF275.lay';
cfg.parameter = 'avg';
cfg.comment = 'no';
cfg.marker = 'off';

subplot(323);
dum = [];
dum.label = meg_data.label;
dum.avg = lead_field{ix(VCidx(1))};
dum.time = 0;
dum.dimord = 'chan_time';
ft_topoplotER(cfg,dum);
ca = caxis;

subplot(324);
dum = [];
dum.label = meg_data.label;
dum.avg = lead_field{ix(VCidx(2))};
dum.time = 0;
dum.dimord = 'chan_time';
ft_topoplotER(cfg,dum);
caxis(ca);

subplot(326);
dum = [];
dum.label = meg_data.label;
dum.avg = lead_field{ix(VCidx(2))}.*15;
dum.time = 0;
dum.dimord = 'chan_time';
ft_topoplotER(cfg,dum);
caxis(ca);