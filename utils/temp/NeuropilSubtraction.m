%%

%load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Fhalo_all_2.mat');
load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Stanage_20210924_F_beh.mat');
load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Stanage_20210924_Fneu_beh.mat');


%%

idx = 1;
wdw = 1:10000;

left = 1;
right = 10000;

%%

this_F = F_beh(idx,wdw);
this_Fneu = Fneu_beh(idx,wdw);
this_FminFneu = this_F - this_Fneu;

figure;

subplot(3,1,1); hold on;
plot(this_F)
for i=1:10;
    xline(paq_beh.sync(i));
end
xlim([left,right])
[rho,pval] = corr(this_F',this_Fneu','Type','Pearson','Rows','Complete')
title(['F, corr(F,Fneu)=',num2str(rho,1)])

subplot(3,1,2)
plot(this_Fneu)
for i=1:10;
    xline(paq_beh.sync(i));
end
xlim([left,right])
title('Fneu')

subplot(3,1,3)
plot(this_FminFneu)
[rho,pval] = corr(this_FminFneu',this_Fneu','Type','Pearson','Rows','Complete')
for i=1:10;
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['F, corr(F-Fneu,Fneu)=',num2str(rho,1)])






%% Adapting SimplePipeline

% get centroids from suite2p output
centroids = {};
centroids.x = nan(length(iscell),1);
centroids.y = nan(length(iscell),1);
for i=1:length(iscell)
    centroids.x(i) = round(s2p_meta.stat{i}.med(2));
    centroids.y(i) = round(s2p_meta.stat{i}.med(1));
end

% 1

% generate halos
roiRadius = 5; % pixels
watershedWidth = 2; % pixels
haloWidth = 2; % multipler of roi radius for alpha of gaussian
removeOverlap = false; % remove overlapping ROI pixels
refineRois = false;
[rois, halos] = dilateCentroids_mod(centroids, roiRadius, watershedWidth, haloWidth, removeOverlap, false);

% subtract cell masks from halos
halos_final_2 = halos;
for j=1:length(halos)
    [~,these_pixels]=setdiff(halos{j}.coords,sub2ind([512,512],s2p_meta.stat{j}.ypix,s2p_meta.stat{j}.xpix));
    halos_final_2{j}.coords = halos{j}.coords(these_pixels);
    halos_final_2{j}.weights = halos{j}.weights(these_pixels);  
    for i=1:length(s2p_meta.stat{j}.xpix)
        halos_final_2{j}.image(s2p_meta.stat{j}.ypix(i),s2p_meta.stat{j}.xpix(i)) = 0;
    end
    temp = floor(length(halos_final_2{j}.coords)/4);
    halos_final_2_1{j}.coords = halos_final_2{j}.coords(1:temp);
    halos_final_2_1{j}.weights = halos_final_2{j}.weights(1:temp);  
    halos_final_2_2{j}.coords = halos_final_2{j}.coords(temp+1:2*temp);
    halos_final_2_2{j}.weights = halos_final_2{j}.weights(temp+1:2*temp);  
    halos_final_2_3{j}.coords = halos_final_2{j}.coords(2*temp+1:3*temp);
    halos_final_2_3{j}.weights = halos_final_2{j}.weights(2*temp+1:3*temp);  
    halos_final_2_4{j}.coords = halos_final_2{j}.coords(3*temp+1:end);
    halos_final_2_4{j}.weights = halos_final_2{j}.weights(3*temp+1:end);  
end

% 2

% generate halos
roiRadius = 5; % pixels
watershedWidth = 2; % pixels
haloWidth = 3; % multipler of roi radius for alpha of gaussian
removeOverlap = false; % remove overlapping ROI pixels
refineRois = false;
[rois, halos] = dilateCentroids_mod(centroids, roiRadius, watershedWidth, haloWidth, removeOverlap, false);

% subtract cell masks from halos
halos_final_3 = halos;
for j=1:length(halos)
    [~,these_pixels]=setdiff(halos{j}.coords,sub2ind([512,512],s2p_meta.stat{j}.ypix,s2p_meta.stat{j}.xpix));
    halos_final_3{j}.coords = halos{j}.coords(these_pixels);
    halos_final_3{j}.weights = halos{j}.weights(these_pixels);  
    for i=1:length(s2p_meta.stat{j}.xpix)
        halos_final_3{j}.image(s2p_meta.stat{j}.ypix(i),s2p_meta.stat{j}.xpix(i)) = 0;
    end
end

% 3 

% generate halos
roiRadius = 5; % pixels
watershedWidth = 2; % pixels
haloWidth = 5; % multipler of roi radius for alpha of gaussian
removeOverlap = false; % remove overlapping ROI pixels
refineRois = false;
[rois, halos] = dilateCentroids_mod(centroids, roiRadius, watershedWidth, haloWidth, removeOverlap, false);

% subtract cell masks from halos
halos_final_5 = halos;
for j=1:length(halos)
    [~,these_pixels]=setdiff(halos{j}.coords,sub2ind([512,512],s2p_meta.stat{j}.ypix,s2p_meta.stat{j}.xpix));
    halos_final_5{j}.coords = halos{j}.coords(these_pixels);
    halos_final_5{j}.weights = halos{j}.weights(these_pixels);  
    for i=1:length(s2p_meta.stat{j}.xpix)
        halos_final_5{j}.image(s2p_meta.stat{j}.ypix(i),s2p_meta.stat{j}.xpix(i)) = 0;
    end
end

% create figure
% idx = 1;
% this_img = zeros(512,512);
% for i=1:length(s2p_meta.stat{idx}.xpix)
%     this_img(s2p_meta.stat{idx}.ypix(i),s2p_meta.stat{idx}.xpix(i)) = 1;
% end
% figure;
% %imagesc(this_img + halos{idx}.image/2)
% imagesc(halos_final_2{idx}.image)
% title('haloWidth=2')


%% Halo figures

idx = 1;
figure;
imagesc(halos_final_2{idx}.image)
title('haloWidth=2')

idx = 1;
figure;
imagesc(halos_final_2_1{idx}.image)
title('haloWidth=2, part 1')

idx = 1;
figure;
imagesc(halos_final_2_2{idx}.image)
title('haloWidth=2, part 2')

idx = 1;
figure;
imagesc(halos_final_2_3{idx}.image)
title('haloWidth=2, part 3')

idx = 1;
figure;
imagesc(halos_final_2_4{idx}.image)
title('haloWidth=2, part 4')


%% Applying background mask on registered video

path.imagingFile = 'H:\Data\2021\2021-09\2021-09-24\Stanage\Imaging\registered_movie.raw';

temp = dir(path.imagingFile);
numFrames = temp.bytes/(512*512*2);
% Fhalo_all_2 = nan(numel(halos_final_2),numFrames);
Fhalo_all_2_1 = nan(numel(halos_final_2_1),numFrames);
Fhalo_all_2_2 = nan(numel(halos_final_2_2),numFrames);
Fhalo_all_2_3 = nan(numel(halos_final_2_3),numFrames);
Fhalo_all_2_4 = nan(numel(halos_final_2_4),numFrames);
% Fhalo_all_3 = nan(numel(halos_final_3),numFrames);
% Fhalo_all_5 = nan(numel(halos_final_5),numFrames);
for i=1:numFrames
    
    % load images one-by-one
    this_fid = fopen(path.imagingFile,'r');
    fseek(this_fid,(i-1)*512*512*2,'bof');
    this_frame = uint16(fread(this_fid,512*512,'uint16',0,'l'));
    frewind(this_fid);
    fclose(this_fid);
    
    % extract halo traces - 2_1
    for j = 1:numel(halos_final_2_1)
        halo_px = this_frame(halos_final_2_1{j}.coords);
        halo_w = halos_final_2_1{j}.weights;
        if all(halo_w==1)
            Fhalo_all_2_1(j,i) = mean(halo_px, 1);
        else
            Fhalo_all_2_1(j,i) = sum(halo_w .* double(halo_px), 1) ./ sum(halo_w, 1);
        end
    end
    
    % extract halo traces - 2_2
    for j = 1:numel(halos_final_2_2)
        halo_px = this_frame(halos_final_2_2{j}.coords);
        halo_w = halos_final_2_2{j}.weights;
        if all(halo_w==1)
            Fhalo_all_2_2(j,i) = mean(halo_px, 1);
        else
            Fhalo_all_2_2(j,i) = sum(halo_w .* double(halo_px), 1) ./ sum(halo_w, 1);
        end
    end
    
    % extract halo traces - 2_3
    for j = 1:numel(halos_final_2_3)
        halo_px = this_frame(halos_final_2_3{j}.coords);
        halo_w = halos_final_2_3{j}.weights;
        if all(halo_w==1)
            Fhalo_all_2_3(j,i) = mean(halo_px, 1);
        else
            Fhalo_all_2_3(j,i) = sum(halo_w .* double(halo_px), 1) ./ sum(halo_w, 1);
        end
    end
    
    % extract halo traces - 2_4
    for j = 1:numel(halos_final_2_4)
        halo_px = this_frame(halos_final_2_4{j}.coords);
        halo_w = halos_final_2_4{j}.weights;
        if all(halo_w==1)
            Fhalo_all_2_4(j,i) = mean(halo_px, 1);
        else
            Fhalo_all_2_4(j,i) = sum(halo_w .* double(halo_px), 1) ./ sum(halo_w, 1);
        end
    end
    
%     % extract halo traces - 2
%     for j = 1:numel(halos_final_2)
%         halo_px = this_frame(halos_final_2{j}.coords);
%         halo_w = halos_final_2{j}.weights;
%         if all(halo_w==1)
%             Fhalo_all_2(j,i) = mean(halo_px, 1);
%         else
%             Fhalo_all_2(j,i) = sum(halo_w .* double(halo_px), 1) ./ sum(halo_w, 1);
%         end
%     end
    
%     % extract halo traces - 3
%     for j = 1:numel(halos_final_3)
%         halo_px = this_frame(halos_final_3{j}.coords);
%         halo_w = halos_final_3{j}.weights;
%         if all(halo_w==1)
%             Fhalo_all_3(j,i) = mean(halo_px, 1);
%         else
%             Fhalo_all_3(j,i) = sum(halo_w .* double(halo_px), 1) ./ sum(halo_w, 1);
%         end
%     end
% 
%     % extract halo traces - 5
%     for j = 1:numel(halos_final_5)
%         halo_px = this_frame(halos_final_5{j}.coords);
%         halo_w = halos_final_5{j}.weights;
%         if all(halo_w==1)
%             Fhalo_all_5(j,i) = mean(halo_px, 1);
%         else
%             Fhalo_all_5(j,i) = sum(halo_w .* double(halo_px), 1) ./ sum(halo_w, 1);
%         end
%     end
    
    if mod(i,100)==0
        disp(num2str(i))
    end
end
save(['E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Fhalo_all_2_1.mat'],'Fhalo_all_2_1','-v7.3');
save(['E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Fhalo_all_2_2.mat'],'Fhalo_all_2_2','-v7.3');
save(['E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Fhalo_all_2_3.mat'],'Fhalo_all_2_3','-v7.3');
save(['E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Fhalo_all_2_4.mat'],'Fhalo_all_2_4','-v7.3');
%save(['E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Fhalo_all_2.mat'],'Fhalo_all_2','-v7.3');
%save(['E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Fhalo_all_3.mat'],'Fhalo_all_3','-v7.3');
% save(['E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Fhalo_all_5.mat'],'Fhalo_all_5','-v7.3');


%%


Fhalo = halo_traces(:,1:10000);







%%


this_F = F_beh(idx,wdw);
this_Fneu = Fneu_beh(idx,wdw);
this_FminFneu = this_F - this_Fneu;
this_Fhalo = halo_traces(idx,wdw);
this_FminFhalo = this_F - this_Fhalo;

figure;

subplot(5,1,1); hold on;
plot(this_F)
for i=1:10;
    xline(paq_beh.sync(i));
end
xlim([left,right])
title('F')

subplot(5,1,2)
plot(this_Fneu)
for i=1:10;
    xline(paq_beh.sync(i));
end
xlim([left,right])
title('Fneu')

subplot(5,1,3)
plot(this_FminFneu)
for i=1:10;
    xline(paq_beh.sync(i));
end
xlim([left,right])
title('F-Fneu')

subplot(5,1,4)
plot(this_Fhalo)
for i=1:10;
    xline(paq_beh.sync(i));
end
xlim([left,right])
title('Fhalo')

subplot(5,1,5)
plot(this_FminFhalo)
for i=1:10;
    xline(paq_beh.sync(i));
end
xlim([left,right])
title('F-Fhalo')

%%


figure; hold on;

plot(this_F,'b')
plot(this_Fhalo,'r')
plot(this_FminFhalo,'g')
for i=1:10;
    xline(paq_beh.sync(i));
end
xlim([left,right])
legend('F','Fhalo','F-Fhalo')


%%

pBaseline          = 1;  % lowest X proportion to use for baseline
limitSub           = false;  % limit subtraction of larger neuropil signals? if true, neuropil signals tha tlare larger than cell traces are set to equal the cell value (so subtraction = zero)
offsetMean         = false;  % offset mean (or baseline)of neuropil-subtrrcted to match the baseline of the raw unsubtracted trace?


haloScaleFactors = estimateNeuropilCoefficients(this_F, this_Fhalo);
processedTraces = haloSubtraction(this_F, this_Fhalo, haloScaleFactors, pBaseline, limitSub, offsetMean);


%%


figure; hold on;

plot(this_F,'b')
plot(this_Fhalo,'r')
plot(this_FminFhalo,'g')
plot(processedTraces,'k')
for i=1:10;
    xline(paq_beh.sync(i));
end
xlim([left,right])
legend('F','Fhalo','F-Fhalo','F-Fhalo*')


%%

figure;
scatter(this_Fhalo,this_F)
xlabel('Fhalo')
ylabel('F')

hold on;
[corr_r,corr_p] = fitLine(this_Fhalo',this_F','k');

%[rho,pval] = corr(this_Fhalo',this_F','Type','Pearson','Rows','Complete')














