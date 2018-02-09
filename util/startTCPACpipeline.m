%% set path defs
restoredefaultpath;

%fieldtrip
ft = dir('~rouxf/tbx/fieldtrip-*');
addpath(['~rouxf/tbx/',ft.name]);
ft_defaults;

%local scripts/functions
addpath(genpath('/home/rouxf/prj/TC/mcode/'));

%% open parpool
open_parpool;

%%
cd('~rouxf/prj/TC/mcode/');

%%
% SCALE FACTOR Matrices: dipole activation strength (nA*m )

% scf = [ 0.2 0.01 0.1;...%1=thalamus,2=parietal cortex,3=occipital cortex
%         0.1 0.01 0.1;...        
%         0.05 0.01 0.1;...                
%         0.01 0.01 0.1];

scf = [ 20   1 10;...
        15   1 10;...
        10 1 10;...%1=thalamus (Pulvinar),2=parietal cortex (BA7),3=occipital cortex (BA18)
         5   1 10;...        
         1 1 10];
    
nDip = [0 100 1e3];% number of random dipoles to simulate 

%parameters to manipulate: 1 coherence, 2 patch size

%% simualte the time series

%simulate_independent_alpha_generatorsRev([1 1 1],1);

%%
for lt = 1%1:2
    for jt = 3%1:size(scf,1)
        for kt = 2%1:length(nDip)
            
            snr = scf(jt,1)/scf(jt,3)% SNR between thalamus and BA 18                        
            
            %simulate_dipole_ROI(snr,lt,nDip(kt),scf(jt,:));
            
            file1 = [ 'sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR',num2str(snr),'_nRandSources',num2str(nDip(kt)),'_simMEG_ROI_fieldMode',num2str(lt),'.mat' ];
            file2 = [ 'sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_nonrhythmic_SNR',num2str(snr),'_nRandSources',num2str(nDip(kt)),'_simMEG_ROI_fieldMode',num2str(lt),'.mat' ];
            
            %contrast_PACvsNoPAC_DICS(file1,file2);
            
            open_parpool;
            bpfreq = [9 11];
            %compute_virtual_channels2(bpfreq,file1);
            
            open_parpool;
            bpfreq = [55 85];
            %compute_virtual_channels2(bpfreq,file1);
            
            path2files = '/home/rouxf/prj/TC/matFilesRev/';           
            [lf_file] = dir([path2files,'virtual_channels_9:11Hz_sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR',num2str(snr),'_nRandSources',num2str(nDip(kt)),'_simMEG_ROI_fieldMode',num2str(lt),'.mat']);
            [hf_file] = dir([path2files,'virtual_channels_55:85Hz_sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_alphagammaPAC_SNR',num2str(snr),'_nRandSources',num2str(nDip(kt)),'_simMEG_ROI_fieldMode',num2str(lt),'.mat']);
            files = [lf_file;hf_file];
            
            load([path2files,file1],'VCidx','gClab');            
            
            if ( length(VCidx) ~= length(gClab) )
                error('index assignment is out of range');
            end;
            
            % local PAC
            open_parpool;
            compute_MI_source_level(1,files);
            
            % thalamic seed region PAC
            open_parpool;
            compute_MI_source_level_seed_region_based(1,files,VCidx(find(gClab== 1)),'L'); %left Thalamus
            open_parpool;
            compute_MI_source_level_seed_region_based(1,files,VCidx(find(gClab== 2)),'R'); %right Thalamus
            
            % reversed-seed region PAC
            c = 0;
            hemi = {'L','R'};
            for it = 3:4
                c = c+1;
                % reversed seed region
                open_parpool;% BA 7 amp -> whole brain phase
                compute_MI_source_level_reversed_seed_region_based(1,files,VCidx(find(gClab== it)),hemi{c});
                % reversed seed region with ROI based orthogonalization
                open_parpool;% BA 7 amp -> whole brain phase
                compute_MI_source_level_reversed_seed_region_based_ortho(1,files,VCidx(find(gClab== it)),hemi{c});
            end;
            
%             open_parpool;
%             compute_MI_sensor_level(file1);
            
        end;
    end;
end;

% exit MATLAB
delete(gcp);
exit;
