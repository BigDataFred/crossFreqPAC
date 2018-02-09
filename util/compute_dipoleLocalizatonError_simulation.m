%%
p2d = '/home/rouxf/prj/TC/matFiles/';

%% calculate DLE
load([p2d,'template_volumes.mat']);
[template_grid] = ft_convert_units(template_grid,'mm');
gIx = find(template_grid.inside);

%%
ixM = [1 3; ...
       2 4;...
       3 1;...
       4 2];
   
   %%
   scf = [ 0.2 0.01 0.1;...%1=thalamus,2=parietal cortex,3=occipital cortex
       0.1 0.01 0.1;...
       0.05 0.01 0.1;...
       0.01 0.01 0.1];
   
   SNR = scf(:,1)./scf(:,3);
 
%%
VCfiles = dir([p2d,'sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_*_SNR1_nRandSources1000_simMEG_ROI*_fieldMode1.mat']);
%VCfiles = VCfiles([1 2 3 4 5 9]);
VCfiles = VCfiles(1:length(VCfiles)/2);
id = {};
for it = 1:length(VCfiles)
    x = regexp(VCfiles(it).name,'ROI');
    id(it) = {VCfiles(it).name(max(x):end-4)};
end;
[~,sIDx] = sort(id);
VCfiles = VCfiles(sIDx);

%%
dleM = [];
dleG = [];
ROC1 = [];
cnt = 0;
c2mem = [];
c1mem = [];
for it = 1:4%:length(VCfiles)
    for kt = 1:size(ixM,1)
        
        cnt = cnt+1;
        load([p2d,VCfiles(it).name],'VCidx');        
        
        %%
        kt
        fn = VCfiles(it).name;
        fn = ['whole_brain_alphaGammaPAC_reversed_seed_based_ortho_',num2str(VCidx(ixM(kt,1))),'_virtual_channels_9:11Hz_',fn];
        fn = dir([p2d,fn]);
        
        if isempty( fn)
            fn = VCfiles(it).name;
            fn = ['whole_brain_alphaGammaPAC_seed_based_',num2str(VCidx(ixM(kt,1))),'_virtual_channels_9:11Hz_',fn];            
            fn = dir([p2d,fn]);
        end;
        
        load([p2d,fn.name],'MI');
        
        %%
        dum = zeros(1,length(template_grid.inside));
        dum(gIx) = MI;
        dum(isnan(dum)) = 0;
        %mask = dum >= max(dum)*0.9;
        %dum = dum.*mask;
        
        %%
        [~,mix] = max(dum);
        
        c1 = template_grid.pos(gIx(VCidx(ixM(kt,2))),:);
        c2 = template_grid.pos(mix,:);
        
        c1mem = [c1mem;c1];
        c2mem = [c2mem;c2];
        
        [dleM(it,kt)] = dipoleLocaliziationError(c1,c2);
        
        %%
        trsh = [0.9:-0.2:0.5];
        for jt = 1:length(trsh)
            mask = dum >= max(dum)*trsh(jt);
            
            dipX = find(mask ==1);
            nd = length(mask);
            td = gIx(VCidx(ixM(kt,2)));
            ROC1(it,kt,jt,:) = [ismember(td,dipX) length(setdiff(dipX,td))/nd];
        end;
        
        %%
         mask = dum >= max(dum)*0.5;
        
%         dum2 = zeros(template_grid.dim);
%         for jt = 1:length(gIx)
%             dum2(gIx(jt)) =mask(gIx(jt));
%         end;
%         gIx2 = find(template_grid.inside==0);
%          for jt = 1:length(gIx2)
%             dum2(gIx2(jt)) =mask(gIx2(jt));
%         end;
        
         dum2 = zeros(template_grid.dim);
        for jt = 1:length(template_grid.inside)
            dum2((jt)) =mask((jt));
        end;
        
        [x,y,z] = COG(dum2);
        
        [indx] = sub2ind(size(dum2),x,y,z);
        [indx] = nearest(1:length(template_grid.inside),indx);
        c1 = template_grid.pos(gIx(VCidx(ixM(kt,2))),:);
        c2 = template_grid.pos(indx,:);
        
        [dleG(it,kt)] = dipoleLocaliziationError(c1,c2);
        
    end;
end;

AUC = [];
for jt = 1:size(ROC1,1)
    for kt = 1:size(ROC1,2)
        x = []; y = [];
        for lt = 1:size(ROC1,3)
            cnt = cnt+1;
            x(lt) = [squeeze(ROC1(jt,kt,lt,2))];
            y(lt) = [squeeze(ROC1(jt,kt,lt,1))];
        end;
        [~,sIx] = sort(x);
        x = x(sIx);
        y = y(sIx);
        if y(end) ==1 && x(end) <1
            y(end+1) = 1;
            x(end+1) = 1;
        end;
        AUC(jt,kt) = trapz(x,y);
%         if AUC(jt,kt) >=.8
%             return;
%         end;
    end;
end;

%%
VCfiles = dir([p2d,'sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_*_SNR1_nRandSources1000_simMEG_ROI*_fieldMode2.mat']);
%VCfiles = VCfiles([1 2 3 4 5 9]);
VCfiles = VCfiles(1:length(VCfiles)/2);
id = {};
for it = 1:length(VCfiles)
    x = regexp(VCfiles(it).name,'ROI');
    id(it) = {VCfiles(it).name(max(x):end-4)};
end;
[~,sIDx] = sort(id);
VCfiles = VCfiles(sIDx);

%%
dleM2 = [];
dleG2 = [];
ROC1 = [];
cnt = 0;
c4mem = [];
c3mem = [];
for it = 1:4%:length(VCfiles)
    for kt = 1:size(ixM,1)
        
        cnt = cnt+1;
        load([p2d,VCfiles(it).name],'VCidx');        
        
        %%
        kt
        fn = VCfiles(it).name;
        fn = ['whole_brain_alphaGammaPAC_reversed_seed_based_ortho_',num2str(VCidx(ixM(kt,1))),'_virtual_channels_9:11Hz_',fn];
        fn = dir([p2d,fn]);
        
        if isempty( fn)
            fn = VCfiles(it).name;
            fn = ['whole_brain_alphaGammaPAC_seed_based_',num2str(VCidx(ixM(kt,1))),'_virtual_channels_9:11Hz_',fn];            
            fn = dir([p2d,fn]);
        end;
        
        load([p2d,fn.name],'MI');
        
        %%
        dum = zeros(1,length(template_grid.inside));
        dum(gIx) = MI;
        dum(isnan(dum)) = 0;
        %mask = dum >= max(dum)*0.9;
        %dum = dum.*mask;
        
        %%
        [~,mix] = max(dum);
        
        c1 = template_grid.pos(gIx(VCidx(ixM(kt,2))),:);
        c2 = template_grid.pos(mix,:);               
        c3mem = [c3mem;c1];
        c4mem = [c4mem;c2];
        
        [dleM2(it,kt)] = dipoleLocaliziationError(c1,c2);
        
        %%
        trsh = [0.9:-0.2:0.5];
        for jt = 1:length(trsh)
            mask = dum >= max(dum)*trsh(jt);
            
            dipX = find(mask ==1);
            nd = length(mask);
            td = gIx(VCidx(ixM(kt,2)));
            ROC1(it,kt,jt,:) = [ismember(td,dipX) length(setdiff(dipX,td))/nd];
        end;
        
        %%
         mask = dum >= max(dum)*0.5;
        
%         dum2 = zeros(template_grid.dim);
%         for jt = 1:length(gIx)
%             dum2(gIx(jt)) =mask(gIx(jt));
%         end;
%         gIx2 = find(template_grid.inside==0);
%          for jt = 1:length(gIx2)
%             dum2(gIx2(jt)) =mask(gIx2(jt));
%         end;
        
         dum2 = zeros(template_grid.dim);
        for jt = 1:length(template_grid.inside)
            dum2((jt)) =mask((jt));
        end;
        
        [x,y,z] = COG(dum2);
        
        [indx] = sub2ind(size(dum2),x,y,z);
        [indx] = nearest(1:length(template_grid.inside),indx);
        c1 = template_grid.pos(gIx(VCidx(ixM(kt,2))),:);
        c2 = template_grid.pos(indx,:);
        
        [dleG2(it,kt)] = dipoleLocaliziationError(c1,c2);
        
    end;
end;

AUC2 = [];
for jt = 1:size(ROC1,1)
    for kt = 1:size(ROC1,2)
        x = []; y = [];
        for lt = 1:size(ROC1,3)
            cnt = cnt+1;
            x(lt) = [squeeze(ROC1(jt,kt,lt,2))];
            y(lt) = [squeeze(ROC1(jt,kt,lt,1))];
        end;
        [~,sIx] = sort(x);
        x = x(sIx);
        y = y(sIx);
        if y(end) ==1 && x(end) <1
            y(end+1) = 1;
            x(end+1) = 1;
        end;
        AUC2(jt,kt) = trapz(x,y);
%         if AUC(jt,kt) >=.8
%             return;
%         end;
    end;
end;

%%
VCfiles = dir([p2d,'sourcePositions_ATLASbasedROI_Pulvinar_Brodmannarea7_Brodmannarea18_dipoles_*_SNR1_nRandSources1000_simMEGII_ROI*_fieldMode1.mat']);
%VCfiles = VCfiles([1 2 3 4 5 9]);
VCfiles = VCfiles(1:length(VCfiles)/2);
id = {};
for it = 1:length(VCfiles)
    x = regexp(VCfiles(it).name,'ROI');
    id(it) = {VCfiles(it).name(max(x):end-4)};
end;
[~,sIDx] = sort(id);
VCfiles = VCfiles(sIDx);

%%
dleM3 = [];
dleG3 = [];
ROC1 = [];
cnt = 0;
c5mem = [];
c6mem = [];
for it = 1:4%:length(VCfiles)
    for kt = 1:size(ixM,1)
        
        cnt = cnt+1;
        load([p2d,VCfiles(it).name],'VCidx');        
        
        %%
        kt
        fn = VCfiles(it).name;
        fn = ['whole_brain_alphaGammaPAC_reversed_seed_based_ortho_',num2str(VCidx(ixM(kt,1))),'_virtual_channels_9:11Hz_',fn];
        fn = dir([p2d,fn]);
        
        if isempty( fn)
            fn = VCfiles(it).name;
            fn = ['whole_brain_alphaGammaPAC_seed_based_',num2str(VCidx(ixM(kt,1))),'_virtual_channels_9:11Hz_',fn];            
            fn = dir([p2d,fn]);
        end;
        
        load([p2d,fn.name],'MI');
        
        %%
        dum = zeros(1,length(template_grid.inside));
        dum(gIx) = MI;
        dum(isnan(dum)) = 0;
        %mask = dum >= max(dum)*0.9;
        %dum = dum.*mask;
        
        %%
        [~,mix] = max(dum);
        
        c1 = template_grid.pos(gIx(VCidx(ixM(kt,2))),:);
        c2 = template_grid.pos(mix,:);
        c5mem = [c5mem;c1];
        c6mem = [c6mem;c2];
        
        [dleM3(it,kt)] = dipoleLocaliziationError(c1,c2);
        
        %%
        trsh = [0.9:-0.2:0.5];
        for jt = 1:length(trsh)
            mask = dum >= max(dum)*trsh(jt);
            
            dipX = find(mask ==1);
            nd = length(mask);
            td = gIx(VCidx(ixM(kt,2)));
            ROC1(it,kt,jt,:) = [ismember(td,dipX) length(setdiff(dipX,td))/nd];
        end;
        
        %%
         mask = dum >= max(dum)*0.5;
        
%         dum2 = zeros(template_grid.dim);
%         for jt = 1:length(gIx)
%             dum2(gIx(jt)) =mask(gIx(jt));
%         end;
%         gIx2 = find(template_grid.inside==0);
%          for jt = 1:length(gIx2)
%             dum2(gIx2(jt)) =mask(gIx2(jt));
%         end;
        
         dum2 = zeros(template_grid.dim);
        for jt = 1:length(template_grid.inside)
            dum2((jt)) =mask((jt));
        end;
        
        [x,y,z] = COG(dum2);
        
        [indx] = sub2ind(size(dum2),x,y,z);
        [indx] = nearest(1:length(template_grid.inside),indx);
        c1 = template_grid.pos(gIx(VCidx(ixM(kt,2))),:);
        c2 = template_grid.pos(indx,:);
        
        [dleG3(it,kt)] = dipoleLocaliziationError(c1,c2);
        
    end;
end;

AUC3 = [];
for jt = 1:size(ROC1,1)
    for kt = 1:size(ROC1,2)
        x = []; y = [];
        for lt = 1:size(ROC1,3)
            cnt = cnt+1;
            x(lt) = [squeeze(ROC1(jt,kt,lt,2))];
            y(lt) = [squeeze(ROC1(jt,kt,lt,1))];
        end;
        [~,sIx] = sort(x);
        x = x(sIx);
        y = y(sIx);
        if y(end) ==1 && x(end) <1
            y(end+1) = 1;
            x(end+1) = 1;
        end;
        AUC3(jt,kt) = trapz(x,y);
%         if AUC(jt,kt) >=.8
%             return;
%         end;
    end;
end;