%% set path defs
restoredefaultpath;

%fieldtrip
ft = dir('~rouxf/tbx/fieldtrip-*');
addpath(['~rouxf/tbx/',ft.name]);
ft_defaults;

%local scripts/functions
addpath(genpath('/home/rouxf/prj/TC/mcode/'));

%%
for jt = 1%:4
    
    [file1] = [ 'sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources1000_simMEG_ROI',num2str(jt),'_fieldMode2.mat' ];
    
    compute_virtual_channels2([1 100],file1);
    
end;