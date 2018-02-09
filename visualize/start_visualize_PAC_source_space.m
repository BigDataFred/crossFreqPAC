%% set path defs
restoredefaultpath;

%fieldtrip
ft = dir('~rouxf/tbx/fieldtrip-*');
addpath(['~rouxf/tbx/',ft.name]);
ft_defaults;

%local scripts/functions
addpath(genpath('/home/rouxf/prj/TC/mcode/'));

cd /home/rouxf/prj/TC/matFilesRev/;

%%
file1 = 'sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources100_simMEG_ROI_fieldMode1.mat';
file2 = 'sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_nonrhythmic_SNR1_nRandSources100_simMEG_ROI_fieldMode1.mat';

visualize_simulated_dipole_on_MEGchannels(file1,file2);

%%
file1 = dir('');

visualize_PAC_sensor_level(file1);

%%
f1 = dir('dics_PAC_sourceMap_*_sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources100_simMEG_ROI_fieldMode1.mat');
f2 = dir('dics_noPAC_sourceMap_*_sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_nonrhythmic_SNR1_nRandSources100_simMEG_ROI_fieldMode1.mat');

files1 ={};
for it = 1:length(f1)
    files1(it) = {f1(it).name};
end;
files2 ={};
for it = 1:length(f2)
    files2(it) = {f2(it).name};
end;

file3 = 'sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources100_simMEG_ROI_fieldMode1.mat';

visualize_contrast_PACvsNoPAC_DICS(files1,files2,file3);

%%
file1 = 'whole_brain_alphaGammaPAC_local_virtual_channels_9:11Hz_sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources100_simMEG_ROI_fieldMode1.mat';
file2 = 'lcmv_spatial_filter_9:11Hz_sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources100_simMEG_ROI_fieldMode1.mat';

visualize_PAC_source_space(file1,file2);

%% L & R Thalamus
file1 = 'whole_brain_alphaGammaPAC_seed_based_5_LH_virtual_channels_9:11Hz_sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources100_simMEG_ROI_fieldMode1.mat';
visualize_PAC_source_space(file1,file2);

file1 = 'whole_brain_alphaGammaPAC_seed_based_5_RH_virtual_channels_9:11Hz_sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources100_simMEG_ROI_fieldMode1.mat';
visualize_PAC_source_space(file1,file2);

%% L cuneus
file1 = 'whole_brain_alphaGammaPAC_reversed_seed_based_23_LH_virtual_channels_9:11Hz_sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources100_simMEG_ROI_fieldMode1.mat';
visualize_PAC_source_space(file1,file2);

file1 = 'whole_brain_alphaGammaPAC_reversed_seed_based_ortho_23_LH_virtual_channels_9:11Hz_sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources100_simMEG_ROI_fieldMode1.mat';
visualize_PAC_source_space(file1,file2);

%% R cuneus
file1 = 'whole_brain_alphaGammaPAC_reversed_seed_based_31_RH_virtual_channels_9:11Hz_sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources100_simMEG_ROI_fieldMode1.mat';
visualize_PAC_source_space(file1,file2);

file1 = 'whole_brain_alphaGammaPAC_reversed_seed_based_ortho_31_RH_virtual_channels_9:11Hz_sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR1_nRandSources100_simMEG_ROI_fieldMode1.mat';
visualize_PAC_source_space(file1,file2);
