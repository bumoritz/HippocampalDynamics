%% Load data

load('D:\SniffinHippo\Repo\Stanage\Stanage_20210924\Stanage_20210924_s2p_meta.mat');
s2p = s2p_meta;
iscell = s2p_meta.iscell(:,1);


%% Identifying background area

p.bkg.prctile = 8;
p.bkg.inset = 20;

this_meanImage = s2p.ops.meanImg;

this_meanImage_inset = this_meanImage(p.bkg.inset+1:size(this_meanImage,1)-p.bkg.inset,p.bkg.inset+1:size(this_meanImage,1)-p.bkg.inset);
this_meanImage_inset_lin = this_meanImage_inset(:);
these_background_idcs = find(this_meanImage_inset_lin<prctile(this_meanImage_inset_lin,p.bkg.prctile));
this_backgroundImage_inset_lin = zeros(size(this_meanImage_inset_lin));
this_backgroundImage_inset_lin(these_background_idcs)=1;
this_backgroundImage_inset = reshape(this_backgroundImage_inset_lin,size(this_meanImage_inset));
this_backgroundImage = zeros(size(this_meanImage));
this_backgroundImage(p.bkg.inset+1:size(this_meanImage,1)-p.bkg.inset,p.bkg.inset+1:size(this_meanImage,1)-p.bkg.inset) = this_backgroundImage_inset;
this_backgroundImage_lin = this_backgroundImage(:);

figure;
subplot(1,2,1)
imshow(this_meanImage,[nanmin(this_meanImage(:)),nanmax(this_meanImage(:))/4]);
title('average image')
subplot(1,2,2)
imshow(this_backgroundImage,[0,1]);
title('background mask')


%% Applying background mask on registered video

path.imagingFile    = 'H:\Data\2021\2021-09\2021-09-24\Stanage\Imaging\registered_movie.raw';

temp = dir(path.imagingFile);
numFrames = temp.bytes/(512*512*2)
background = nan(1,numFrames);
for i=1:numFrames
    
    % load images one-by-one
    this_fid = fopen(path.imagingFile,'r');
    fseek(this_fid,(i-1)*512*512*2,'bof');
    this_frame = uint16(fread(this_fid,512*512,'uint16',0,'l'));
    %temp = reshape(temp,512,512);
    frewind(this_fid);
    fclose(this_fid);
    
    % calculate background
    background(i) = nanmean(this_frame(find(this_backgroundImage_lin)));
    
    if mod(i,10000)==0
        disp(num2str(i))
    end
end
save(['E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Stanage_20210924_background_prctile20.mat'],'background','-v7.3');

% load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\background.mat')
% background_beh = background(108001:end-108000);

%% Load F

load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Stanage_20210924_F_beh.mat');
load('E:\SniffinHippo\RepoX\Stanage\Stanage_20210924\Stanage_20210924_Fneu_beh.mat');
load('D:\SniffinHippo\Repo\Stanage\Stanage_20210924\Stanage_20210924_paq_beh.mat');


%% --- Area-wise neuropil subtraction --- (subtracting a spatially smooth neuropil signal)

% parameters
numBins = 64;
numPixelsPerDim = 512;

% identify centroids of each area
areaBoundaries = (0:numBins:numPixelsPerDim)+0.5;
areaCentres = numBins/2+0.5:numBins:numPixelsPerDim;

% assign each neuropil trace to an area
Fneu_area_traces = cell(sqrt(numBins),sqrt(numBins));
for i=1:size(s2p_meta.iscell,1)
    if ismember(i,find(s2p_meta.iscell(:,1)))
        m = discretize(s2p_meta.stat{i}.med(1),areaBoundaries);
        n = discretize(s2p_meta.stat{i}.med(2),areaBoundaries);
        Fneu_area_traces{m,n} = [Fneu_area_traces{m,n}; Fneu_beh(i,:)];
    end
end

% calculate average neuropil trace for each area
Fneu_area_avg = cell(sqrt(numBins),sqrt(numBins)); %nan(sqrt(numBins),sqrt(numBins),size(Fneu_beh,2));
for m=1:sqrt(numBins)
    for n=1:sqrt(numBins)
        Fneu_area_avg{m,n} = nanmean(Fneu_area_traces{m,n},1);
        if isempty(Fneu_area_avg{m,n})
            Fneu_area_avg{m,n} = nan(1,size(Fneu_beh,2));
        end
    end
end


%%

Fns = nan(size(F_beh));
for i=1:size(s2p_meta.iscell,1)
    if ismember(i,find(s2p_meta.iscell(:,1)))
        m = discretize(s2p_meta.stat{i}.med(1),areaBoundaries);
        n = discretize(s2p_meta.stat{i}.med(2),areaBoundaries);             
        Fns(i,:) = F_beh(i,:) - Fneu_area_avg{m,n};
    end
end



%%

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
sync_beh.sync = paq_beh.sync; 

nrows = sqrt(numBins); ncols = sqrt(numBins);
F = default_figure([20,0.5,20,9.9]);

for m=1:sqrt(numBins)
    for n=1:sqrt(numBins)
        r=m; c=n; subplot(nrows,ncols,(r-1)*ncols+c);
        
        act = repmat(Fneu_area_avg{m,n},2000,1); 
        [~,~,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);
        this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

        plot(this_data)
        title(['row=',num2str(m),', col=',num2str(n)])
        m
        n
    end
end








%%

% Fneu_area_interpol = nan(numPixelsPerDim,numPixelsPerDim,size(Fneu_beh,2));




%%













%% --- PLOTTING MUA AVERAGES ---

act = F_beh; %repmat(background_beh,2000,1); % F_beh - 
sync_beh.sync = paq_beh.sync; 

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('F')

% 

act = Fneu_beh; %repmat(background_beh,2000,1); % F_beh - 
sync_beh.sync = paq_beh.sync; 

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('Fneu')

% 

act = F_beh - Fneu_beh; %repmat(background_beh,2000,1); % F_beh - 
sync_beh.sync = paq_beh.sync; 

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('F-Fneu')

% 

act = F_beh - 0.7*Fneu_beh; %repmat(background_beh,2000,1); % F_beh - 
sync_beh.sync = paq_beh.sync; 

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('F-0.7*Fneu')


% 

act = Fhalo_all_2(:,108001:end-108000); %repmat(background_beh,2000,1); % F_beh - 
sync_beh.sync = paq_beh.sync; 

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('Fhalo2')

% 

act = F_beh - Fhalo_all_2(:,108001:end-108000); %repmat(background_beh,2000,1); % F_beh - 
sync_beh.sync = paq_beh.sync; 

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(find(iscell==1),:,:),3),1);

figure;
plot(this_data)
title('F-Fhalo2')




%% --- PLOTTING SINGLE CELL AVERAGES ---

idx = 1;

act = F_beh; %repmat(background_beh,2000,1); % F_beh - 
sync_beh.sync = paq_beh.sync; 

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(idx,:,:),3),1);

figure;
plot(this_data)
title('F')

% 

act = Fneu_beh; %repmat(background_beh,2000,1); % F_beh - 
sync_beh.sync = paq_beh.sync; 

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(idx,:,:),3),1);

figure;
plot(this_data)
title('Fneu')

% 

act = F_beh - Fneu_beh; %repmat(background_beh,2000,1); % F_beh - 
sync_beh.sync = paq_beh.sync; 

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(idx,:,:),3),1);

figure;
plot(this_data)
title('F-Fneu')

% 

act = F_beh - 0.7*Fneu_beh; %repmat(background_beh,2000,1); % F_beh - 
sync_beh.sync = paq_beh.sync; 

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(idx,:,:),3),1);

figure;
plot(this_data)
title('F-0.7*Fneu')


% 

act = Fhalo_all_2(:,108001:end-108000); %repmat(background_beh,2000,1); % F_beh - 
sync_beh.sync = paq_beh.sync; 

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(idx,:,:),3),1);

figure;
plot(this_data)
title('Fhalo2')

% 

act = F_beh - Fhalo_all_2(:,108001:end-108000); %repmat(background_beh,2000,1); % F_beh - 
sync_beh.sync = paq_beh.sync; 

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(idx,:,:),3),1);

figure;
plot(this_data)
title('F-Fhalo2')




%% RESCALE HALO SIGNAL


idx = 1;

act = F_beh; %repmat(background_beh,2000,1); % F_beh - 
sync_beh.sync = paq_beh.sync; 

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 1;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(idx,:,:),3),1);

figure;
plot(this_data)
title('F')

% 

act = Fhalo_all_2(:,108001:end-108000); %repmat(background_beh,2000,1); % F_beh - 
sync_beh.sync = paq_beh.sync; 

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 1;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(idx,:,:),3),1);

figure;
plot(this_data)
title('Fhalo2')

% 

act = F_beh - Fhalo_all_2(:,108001:end-108000); %repmat(background_beh,2000,1); % F_beh - 
sync_beh.sync = paq_beh.sync; 

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 1;
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);

this_data = nanmean(nanmean(nft_binned(idx,:,:),3),1);

figure;
plot(this_data)
title('F-Fhalo2')










%% RESCALE HALO SIGNAL ---NEW---

idx = 5;

p = get_p;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
sync_beh.sync = paq_beh.sync; 

data_F = F_beh;
[~,~,nft_binned_F] = preprocessActivityMeasure(data_F,p.inh,p,sync_beh,iscell);
avgData_F = nanmean(nanmean(nft_binned_F(idx,:,:),3),1);

data_Fhalo = Fhalo_all_2(:,108001:end-108000);
[~,~,nft_binned_Fhalo] = preprocessActivityMeasure(data_Fhalo,p.inh,p,sync_beh,iscell);
avgData_Fhalo = nanmean(nanmean(nft_binned_Fhalo(idx,:,:),3),1);

data_F_Fhalo = F_beh - Fhalo_all_2(:,108001:end-108000);
[~,~,nft_binned_F_Fhalo] = preprocessActivityMeasure(data_F_Fhalo,p.inh,p,sync_beh,iscell);
avgData_F_Fhalo = nanmean(nanmean(nft_binned_F_Fhalo(idx,:,:),3),1);

% estimating neuropil coefficients
pBaseline          = 1;  % lowest X proportion to use for baseline
limitSub           = false;  % limit subtraction of larger neuropil signals? if true, neuropil signals tha tlare larger than cell traces are set to equal the cell value (so subtraction = zero)
offsetMean         = false;  % offset mean (or baseline)of neuropil-subtrrcted to match the baseline of the raw unsubtracted trace?
haloScaleFactors = estimateNeuropilCoefficients(data_F, data_Fhalo);
data_F_estFhalo = haloSubtraction(data_F, data_Fhalo, haloScaleFactors, pBaseline, limitSub, offsetMean);
[~,~,nft_binned_F_estFhalo] = preprocessActivityMeasure(data_F_estFhalo,p.inh,p,sync_beh,iscell);
avgData_F_estFhalo = nanmean(nanmean(nft_binned_F_estFhalo(idx,:,:),3),1);


%% Loop over subtraction factors - avg

for i=0:20

subtractionFactor = i;

figure;

subplot(3,1,1)
plot(avgData_F)
title('F')

subplot(3,1,2)
plot(avgData_Fhalo)
title('Fhalo')

subplot(3,1,3)
plot(avgData_F-subtractionFactor*avgData_Fhalo)
[rho,pval]=corr((avgData_F-subtractionFactor*avgData_Fhalo)',avgData_Fhalo','Type','Pearson','Rows','Complete')
title(['F-',num2str(subtractionFactor),'*Fhalo, corr(Fns,Fhalo)=',num2str(rho,1)]);

end
% subplot(4,1,4)
% plot(avgData_F_estFhalo)
% title(['F-Fhalo (factor=',num2str(haloScaleFactors(idx)),')'])


%% -----

%% Loop over subtraction factors - single trial

data_F = nanmean(nft_binned_F,3);
data_Fhalo = nanmean(nft_binned_Fhalo,3);

subtractionFactors = -20:0.1:20;
rho = nan(size(data_F,1),length(subtractionFactors));
for i=1:size(data_F,1)
    for j=1:length(subtractionFactors)
        subtractionFactor = subtractionFactors(j);
        rho(i,j) = corr((data_F(i,:)-subtractionFactor*data_Fhalo(i,:))',data_Fhalo(i,:)','Type','Pearson','Rows','Complete');
    end
    i
end

optimalSubtractionFactor = nan(size(data_F,1),1);
[temp1,temp2] = nanmin(abs(rho),[],2);
for i=1:size(data_F,1)
    if ~isnan(temp1(i)) && ~isnan(temp2(i))
        optimalSubtractionFactor(i) = subtractionFactors(temp2(i));
    end
end


%% Plot examples

idx = 6;

figure;

subplot(3,1,1)
plot(data_F(idx,:))
title(['idx=',num2str(idx),', F'])

subplot(3,1,2)
plot(data_Fhalo(idx,:))
title(['idx=',num2str(idx),', Fhalo'])

subplot(3,1,3)
plot(data_F(idx,:)-optimalSubtractionFactor(idx)*data_Fhalo(idx,:))
[rho,pval]=corr((data_F(idx,:)-optimalSubtractionFactor(idx)*data_Fhalo(idx,:))',data_Fhalo(idx,:)','Type','Pearson','Rows','Complete')
title(['idx=',num2str(idx),', F-',num2str(optimalSubtractionFactor(idx)),'*Fhalo']);


%% -----














%%





%% Loop over subtraction factors - single trial

rho = nan(length(0:20),1);
for i=0:20
    subtractionFactor = i;
    data_neuropil = smoothdata(data_Fhalo(idx,:),2,'gaussian', 12 *5);
    data_subtracted = smoothdata((data_F(idx,:)-subtractionFactor*data_Fhalo(idx,:)),2,'gaussian', 12 *5);
    rho(i+1) = corr(data_subtracted',data_neuropil','Type','Pearson','Rows','Complete');
end



%%

