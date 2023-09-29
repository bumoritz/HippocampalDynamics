%% Based on this stuff

% plt = struct(); plt.split = 'stimO_stimM_catch'; plt.split3_numTrials = 'same';
% F = snakePlots(traces_all,{rmmissing(temp1(:)),rmmissing(temp2(:))},[],[1;4],[],p,info,plt);


%% Get data

temp = trg.idcs_targetedCells(:,trg.sequenceClusters(:,1)); 
data_Aseq_Atrials_catch = traces_all.A_var0(rmmissing(temp(:)),p.general.bins_analysisWindow,:);
data_Aseq_Atrials_stim = traces_all.A_var1(rmmissing(temp(:)),p.general.bins_analysisWindow,:);
data_Aseq_Xtrials_catch = traces_all.X_var0(rmmissing(temp(:)),p.general.bins_analysisWindow,:);
data_Aseq_Xtrials_stim = traces_all.X_var1(rmmissing(temp(:)),p.general.bins_analysisWindow,:);

temp = trg.idcs_targetedCells(:,trg.sequenceClusters(:,2)); 
data_Xseq_Atrials_catch = traces_all.A_var0(rmmissing(temp(:)),p.general.bins_analysisWindow,:);
data_Xseq_Atrials_stim = traces_all.A_var1(rmmissing(temp(:)),p.general.bins_analysisWindow,:);
data_Xseq_Xtrials_catch = traces_all.X_var0(rmmissing(temp(:)),p.general.bins_analysisWindow,:);
data_Xseq_Xtrials_stim = traces_all.X_var1(rmmissing(temp(:)),p.general.bins_analysisWindow,:);


%%

this_data = data_Aseq_Atrials_stim;

numCells = size(this_data,1);
numBins = size(this_data,2);
numTrials = size(this_data,3);

lags = -numBins+1:numBins-1;
this_xcorr = nan(numCells,numCells,numBins*2-1,numTrials);
for k=1:numTrials
    for i=1:numCells
        for j=1:numCells
            this_xcorr(i,j,:,k) = xcorr(this_data(i,:,k),this_data(j,:,k),'coeff');
        end
    end
    k
end

xcorr_Aseq_Atrials_stim = this_xcorr;
avgxcorr_Aseq_Atrials_stim = nanmean(xcorr_Aseq_Atrials_stim,4);


this_data = data_Aseq_Atrials_catch;

numCells = size(this_data,1);
numBins = size(this_data,2);
numTrials = size(this_data,3);

lags = -numBins+1:numBins-1;
this_xcorr = nan(numCells,numCells,numBins*2-1,numTrials);
for k=1:numTrials
    for i=1:numCells
        for j=1:numCells
            this_xcorr(i,j,:,k) = xcorr(this_data(i,:,k),this_data(j,:,k),'coeff');
        end
    end
    k
end

xcorr_Aseq_Atrials_catch = this_xcorr;
avgxcorr_Aseq_Atrials_catch = nanmean(xcorr_Aseq_Atrials_catch,4);


this_data = data_Aseq_Xtrials_stim;

numCells = size(this_data,1);
numBins = size(this_data,2);
numTrials = size(this_data,3);

lags = -numBins+1:numBins-1;
this_xcorr = nan(numCells,numCells,numBins*2-1,numTrials);
for k=1:numTrials
    for i=1:numCells
        for j=1:numCells
            this_xcorr(i,j,:,k) = xcorr(this_data(i,:,k),this_data(j,:,k),'coeff');
        end
    end
    k
end

xcorr_Aseq_Xtrials_stim = this_xcorr;
avgxcorr_Aseq_Xtrials_stim = nanmean(xcorr_Aseq_Xtrials_stim,4);


this_data = data_Aseq_Xtrials_catch;

numCells = size(this_data,1);
numBins = size(this_data,2);
numTrials = size(this_data,3);

lags = -numBins+1:numBins-1;
this_xcorr = nan(numCells,numCells,numBins*2-1,numTrials);
for k=1:numTrials
    for i=1:numCells
        for j=1:numCells
            this_xcorr(i,j,:,k) = xcorr(this_data(i,:,k),this_data(j,:,k),'coeff');
        end
    end
    k
end

xcorr_Aseq_Xtrials_catch = this_xcorr;
avgxcorr_Aseq_Xtrials_catch = nanmean(xcorr_Aseq_Xtrials_catch,4);


% ---

this_data = data_Xseq_Atrials_stim;

numCells = size(this_data,1);
numBins = size(this_data,2);
numTrials = size(this_data,3);

lags = -numBins+1:numBins-1;
this_xcorr = nan(numCells,numCells,numBins*2-1,numTrials);
for k=1:numTrials
    for i=1:numCells
        for j=1:numCells
            this_xcorr(i,j,:,k) = xcorr(this_data(i,:,k),this_data(j,:,k),'coeff');
        end
    end
    k
end

xcorr_Xseq_Atrials_stim = this_xcorr;
avgxcorr_Xseq_Atrials_stim = nanmean(xcorr_Xseq_Atrials_stim,4);


this_data = data_Xseq_Atrials_catch;

numCells = size(this_data,1);
numBins = size(this_data,2);
numTrials = size(this_data,3);

lags = -numBins+1:numBins-1;
this_xcorr = nan(numCells,numCells,numBins*2-1,numTrials);
for k=1:numTrials
    for i=1:numCells
        for j=1:numCells
            this_xcorr(i,j,:,k) = xcorr(this_data(i,:,k),this_data(j,:,k),'coeff');
        end
    end
    k
end

xcorr_Xseq_Atrials_catch = this_xcorr;
avgxcorr_Xseq_Atrials_catch = nanmean(xcorr_Xseq_Atrials_catch,4);


this_data = data_Xseq_Xtrials_stim;

numCells = size(this_data,1);
numBins = size(this_data,2);
numTrials = size(this_data,3);

lags = -numBins+1:numBins-1;
this_xcorr = nan(numCells,numCells,numBins*2-1,numTrials);
for k=1:numTrials
    for i=1:numCells
        for j=1:numCells
            this_xcorr(i,j,:,k) = xcorr(this_data(i,:,k),this_data(j,:,k),'coeff');
        end
    end
    k
end

xcorr_Xseq_Xtrials_stim = this_xcorr;
avgxcorr_Xseq_Xtrials_stim = nanmean(xcorr_Xseq_Xtrials_stim,4);


this_data = data_Xseq_Xtrials_catch;

numCells = size(this_data,1);
numBins = size(this_data,2);
numTrials = size(this_data,3);

lags = -numBins+1:numBins-1;
this_xcorr = nan(numCells,numCells,numBins*2-1,numTrials);
for k=1:numTrials
    for i=1:numCells
        for j=1:numCells
            this_xcorr(i,j,:,k) = xcorr(this_data(i,:,k),this_data(j,:,k),'coeff');
        end
    end
    k
end

xcorr_Xseq_Xtrials_catch = this_xcorr;
avgxcorr_Xseq_Xtrials_catch = nanmean(xcorr_Xseq_Xtrials_catch,4);


%% Calculate normal correlation without lag

this_data = data_Aseq_Atrials_stim;
this_corr = nan(size(this_data,1),size(this_data,1),size(this_data,3));
for k=1:size(this_data,3)
    this_corr(:,:,k) = corr(this_data(:,:,k)',this_data(:,:,k)','Type','Pearson','Rows','Complete');
end
corr_Aseq_Atrials_stim = this_corr;
avgcorr_Aseq_Atrials_stim = nanmean(corr_Aseq_Atrials_stim,3);

this_data = data_Aseq_Atrials_catch;
this_corr = nan(size(this_data,1),size(this_data,1),size(this_data,3));
for k=1:size(this_data,3)
    this_corr(:,:,k) = corr(this_data(:,:,k)',this_data(:,:,k)','Type','Pearson','Rows','Complete');
end
corr_Aseq_Atrials_catch= this_corr;
avgcorr_Aseq_Atrials_catch = nanmean(corr_Aseq_Atrials_catch,3);

this_data = data_Aseq_Xtrials_stim;
this_corr = nan(size(this_data,1),size(this_data,1),size(this_data,3));
for k=1:size(this_data,3)
    this_corr(:,:,k) = corr(this_data(:,:,k)',this_data(:,:,k)','Type','Pearson','Rows','Complete');
end
corr_Aseq_Xtrials_stim = this_corr;
avgcorr_Aseq_Xtrials_stim = nanmean(corr_Aseq_Xtrials_stim,3);

this_data = data_Aseq_Xtrials_catch;
this_corr = nan(size(this_data,1),size(this_data,1),size(this_data,3));
for k=1:size(this_data,3)
    this_corr(:,:,k) = corr(this_data(:,:,k)',this_data(:,:,k)','Type','Pearson','Rows','Complete');
end
corr_Aseq_Xtrials_catch= this_corr;
avgcorr_Aseq_Xtrials_catch = nanmean(corr_Aseq_Xtrials_catch,3);

% --

this_data = data_Xseq_Atrials_stim;
this_corr = nan(size(this_data,1),size(this_data,1),size(this_data,3));
for k=1:size(this_data,3)
    this_corr(:,:,k) = corr(this_data(:,:,k)',this_data(:,:,k)','Type','Pearson','Rows','Complete');
end
corr_Xseq_Atrials_stim = this_corr;
avgcorr_Xseq_Atrials_stim = nanmean(corr_Xseq_Atrials_stim,3);

this_data = data_Xseq_Atrials_catch;
this_corr = nan(size(this_data,1),size(this_data,1),size(this_data,3));
for k=1:size(this_data,3)
    this_corr(:,:,k) = corr(this_data(:,:,k)',this_data(:,:,k)','Type','Pearson','Rows','Complete');
end
corr_Xseq_Atrials_catch= this_corr;
avgcorr_Xseq_Atrials_catch = nanmean(corr_Xseq_Atrials_catch,3);

this_data = data_Xseq_Xtrials_stim;
this_corr = nan(size(this_data,1),size(this_data,1),size(this_data,3));
for k=1:size(this_data,3)
    this_corr(:,:,k) = corr(this_data(:,:,k)',this_data(:,:,k)','Type','Pearson','Rows','Complete');
end
corr_Xseq_Xtrials_stim = this_corr;
avgcorr_Xseq_Xtrials_stim = nanmean(corr_Xseq_Xtrials_stim,3);

this_data = data_Xseq_Xtrials_catch;
this_corr = nan(size(this_data,1),size(this_data,1),size(this_data,3));
for k=1:size(this_data,3)
    this_corr(:,:,k) = corr(this_data(:,:,k)',this_data(:,:,k)','Type','Pearson','Rows','Complete');
end
corr_Xseq_Xtrials_catch= this_corr;
avgcorr_Xseq_Xtrials_catch = nanmean(corr_Xseq_Xtrials_catch,3);


%% Cross-correlation at lag zero

plt.colormap = 'jet'; % jet looks clearer than parula here

figure;

subplot(2,4,1)
imagesc(avgcorr_Aseq_Atrials_stim);
colormap(gca,plt.colormap);
title(['Aseq-Atrials-stim, ',num2str(nanmean(avgcorr_Aseq_Atrials_stim(:)),2)])

subplot(2,4,2)
imagesc(avgcorr_Aseq_Atrials_catch);
colormap(gca,plt.colormap);
title(['Aseq-Atrials-catch, ',num2str(nanmean(avgcorr_Aseq_Atrials_catch(:)),2)])

subplot(2,4,3)
imagesc(avgcorr_Aseq_Xtrials_stim);
colormap(gca,plt.colormap);
title(['Aseq-Xtrials-stim, ',num2str(nanmean(avgcorr_Aseq_Xtrials_stim(:)),2)])

subplot(2,4,4)
imagesc(avgcorr_Aseq_Xtrials_catch);
colormap(gca,plt.colormap);
title(['Aseq-Xtrials-catch, ',num2str(nanmean(avgcorr_Aseq_Xtrials_catch(:)),2)])

subplot(2,4,5)
imagesc(avgcorr_Xseq_Atrials_stim);
colormap(gca,plt.colormap);
title(['Xseq-Atrials-stim, ',num2str(nanmean(avgcorr_Xseq_Atrials_stim(:)),2)])

subplot(2,4,6)
imagesc(avgcorr_Xseq_Atrials_catch);
colormap(gca,plt.colormap);
title(['Xseq-Atrials-catch, ',num2str(nanmean(avgcorr_Xseq_Atrials_catch(:)),2)])

subplot(2,4,7)
imagesc(avgcorr_Xseq_Xtrials_stim);
colormap(gca,plt.colormap);
title(['Xseq-Xtrials-stim, ',num2str(nanmean(avgcorr_Xseq_Xtrials_stim(:)),2)])

subplot(2,4,8)
imagesc(avgcorr_Xseq_Xtrials_catch);
colormap(gca,plt.colormap);
title(['Xseq-Xtrials-catch, ',num2str(nanmean(avgcorr_Xseq_Xtrials_catch(:)),2)])


%% Things to try

% SEQUENCE IMPRINTING
% - activity of stimmed neurons higher than non-stimmed neurons?
% - field induction analysis (more new fields induced?)
% - correlations only between co-stimmed neurons
% - cross-correlations between neurons of adjacent clusters
% - decoding: train on seq stim, predict on catch
% - decoding: train on catch, predict on catch (compare this between seq stim and ctrl stim sessions)

% - (more fields than by chance? stim vs ctrl, run seq analysis on catch trials)
% - pop vector corr


%% Correlations between co-stimmed neurons





