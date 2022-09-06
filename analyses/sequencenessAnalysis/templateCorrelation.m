function [tempCorr] = templateCorrelation(prop,traces,avgTraces,tng,p)
%traces = traces_all; avgTraces = avgTraces_all; tng = tng_all;

%% Create templates

% iscells
[tempCorr.templateRanks.iscells_Atrials,tempCorr.templateMaxBins.iscells_Atrials] = createTemplate(avgTraces.A,find(prop.iscell),p);
[tempCorr.templateRanks.iscells_Xtrials,tempCorr.templateMaxBins.iscells_Xtrials] = createTemplate(avgTraces.X,find(prop.iscell),p);

% A cells
[tempCorr.templateRanks.A_Atrials,tempCorr.templateMaxBins.A_Atrials] = createTemplate(avgTraces.A,find(tng.passed.AW.A==1),p);
[tempCorr.templateRanks.A_Xtrials,tempCorr.templateMaxBins.A_Xtrials] = createTemplate(avgTraces.X,find(tng.passed.AW.A==1),p);

% X cells
[tempCorr.templateRanks.X_Atrials,tempCorr.templateMaxBins.X_Atrials] = createTemplate(avgTraces.A,find(tng.passed.AW.X==1),p);
[tempCorr.templateRanks.X_Xtrials,tempCorr.templateMaxBins.X_Xtrials] = createTemplate(avgTraces.X,find(tng.passed.AW.X==1),p);

% Aonly cells
[tempCorr.templateRanks.Aonly_Atrials,tempCorr.templateMaxBins.Aonly_Atrials] = createTemplate(avgTraces.A,find(tng.passed.AW.Aonly==1),p);
[tempCorr.templateRanks.Aonly_Xtrials,tempCorr.templateMaxBins.Aonly_Xtrials] = createTemplate(avgTraces.X,find(tng.passed.AW.Aonly==1),p);

% Xonly cells
[tempCorr.templateRanks.Xonly_Atrials,tempCorr.templateMaxBins.Xonly_Atrials] = createTemplate(avgTraces.A,find(tng.passed.AW.Xonly==1),p);
[tempCorr.templateRanks.Xonly_Xtrials,tempCorr.templateMaxBins.Xonly_Xtrials] = createTemplate(avgTraces.X,find(tng.passed.AW.Xonly==1),p);


%% Calculate rank correlations

% iscells
[tempCorr.rankCorr.Pearson.iscells_Atemp_Atrials,tempCorr.rankCorr.Spearman.iscells_Atemp_Atrials,tempCorr.rankCorr.Kendall.iscells_Atemp_Atrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.iscells_Atrials,traces.A,find(prop.iscell),p);
[tempCorr.rankCorr.Pearson.iscells_Atemp_Xtrials,tempCorr.rankCorr.Spearman.iscells_Atemp_Xtrials,tempCorr.rankCorr.Kendall.iscells_Atemp_Xtrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.iscells_Atrials,traces.X,find(prop.iscell),p);
[tempCorr.rankCorr.Pearson.iscells_Xtemp_Atrials,tempCorr.rankCorr.Spearman.iscells_Xtemp_Atrials,tempCorr.rankCorr.Kendall.iscells_Xtemp_Atrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.iscells_Xtrials,traces.A,find(prop.iscell),p);
[tempCorr.rankCorr.Pearson.iscells_Xtemp_Xtrials,tempCorr.rankCorr.Spearman.iscells_Xtemp_Xtrials,tempCorr.rankCorr.Kendall.iscells_Xtemp_Xtrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.iscells_Xtrials,traces.X,find(prop.iscell),p);

% A cells
[tempCorr.rankCorr.Pearson.A_Atemp_Atrials,tempCorr.rankCorr.Spearman.A_Atemp_Atrials,tempCorr.rankCorr.Kendall.A_Atemp_Atrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.A_Atrials,traces.A,find(tng.passed.AW.A==1),p);
[tempCorr.rankCorr.Pearson.A_Atemp_Xtrials,tempCorr.rankCorr.Spearman.A_Atemp_Xtrials,tempCorr.rankCorr.Kendall.A_Atemp_Xtrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.A_Atrials,traces.X,find(tng.passed.AW.A==1),p);
[tempCorr.rankCorr.Pearson.A_Xtemp_Atrials,tempCorr.rankCorr.Spearman.A_Xtemp_Atrials,tempCorr.rankCorr.Kendall.A_Xtemp_Atrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.A_Xtrials,traces.A,find(tng.passed.AW.A==1),p);
[tempCorr.rankCorr.Pearson.A_Xtemp_Xtrials,tempCorr.rankCorr.Spearman.A_Xtemp_Xtrials,tempCorr.rankCorr.Kendall.A_Xtemp_Xtrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.A_Xtrials,traces.X,find(tng.passed.AW.A==1),p);

% X cells
[tempCorr.rankCorr.Pearson.X_Atemp_Atrials,tempCorr.rankCorr.Spearman.X_Atemp_Atrials,tempCorr.rankCorr.Kendall.X_Atemp_Atrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.X_Atrials,traces.A,find(tng.passed.AW.X==1),p);
[tempCorr.rankCorr.Pearson.X_Atemp_Xtrials,tempCorr.rankCorr.Spearman.X_Atemp_Xtrials,tempCorr.rankCorr.Kendall.X_Atemp_Xtrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.X_Atrials,traces.X,find(tng.passed.AW.X==1),p);
[tempCorr.rankCorr.Pearson.X_Xtemp_Atrials,tempCorr.rankCorr.Spearman.X_Xtemp_Atrials,tempCorr.rankCorr.Kendall.X_Xtemp_Atrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.X_Xtrials,traces.A,find(tng.passed.AW.X==1),p);
[tempCorr.rankCorr.Pearson.X_Xtemp_Xtrials,tempCorr.rankCorr.Spearman.X_Xtemp_Xtrials,tempCorr.rankCorr.Kendall.X_Xtemp_Xtrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.X_Xtrials,traces.X,find(tng.passed.AW.X==1),p);

% Aonly cells
[tempCorr.rankCorr.Pearson.Aonly_Atemp_Atrials,tempCorr.rankCorr.Spearman.Aonly_Atemp_Atrials,tempCorr.rankCorr.Kendall.Aonly_Atemp_Atrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.Aonly_Atrials,traces.A,find(tng.passed.AW.Aonly==1),p);
[tempCorr.rankCorr.Pearson.Aonly_Atemp_Xtrials,tempCorr.rankCorr.Spearman.Aonly_Atemp_Xtrials,tempCorr.rankCorr.Kendall.Aonly_Atemp_Xtrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.Aonly_Atrials,traces.X,find(tng.passed.AW.Aonly==1),p);
[tempCorr.rankCorr.Pearson.Aonly_Xtemp_Atrials,tempCorr.rankCorr.Spearman.Aonly_Xtemp_Atrials,tempCorr.rankCorr.Kendall.Aonly_Xtemp_Atrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.Aonly_Xtrials,traces.A,find(tng.passed.AW.Aonly==1),p);
[tempCorr.rankCorr.Pearson.Aonly_Xtemp_Xtrials,tempCorr.rankCorr.Spearman.Aonly_Xtemp_Xtrials,tempCorr.rankCorr.Kendall.Aonly_Xtemp_Xtrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.Aonly_Xtrials,traces.X,find(tng.passed.AW.Aonly==1),p);

% Xonly cells
[tempCorr.rankCorr.Pearson.Xonly_Atemp_Atrials,tempCorr.rankCorr.Spearman.Xonly_Atemp_Atrials,tempCorr.rankCorr.Kendall.Xonly_Atemp_Atrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.Xonly_Atrials,traces.A,find(tng.passed.AW.Xonly==1),p);
[tempCorr.rankCorr.Pearson.Xonly_Atemp_Xtrials,tempCorr.rankCorr.Spearman.Xonly_Atemp_Xtrials,tempCorr.rankCorr.Kendall.Xonly_Atemp_Xtrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.Xonly_Atrials,traces.X,find(tng.passed.AW.Xonly==1),p);
[tempCorr.rankCorr.Pearson.Xonly_Xtemp_Atrials,tempCorr.rankCorr.Spearman.Xonly_Xtemp_Atrials,tempCorr.rankCorr.Kendall.Xonly_Xtemp_Atrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.Xonly_Xtrials,traces.A,find(tng.passed.AW.Xonly==1),p);
[tempCorr.rankCorr.Pearson.Xonly_Xtemp_Xtrials,tempCorr.rankCorr.Spearman.Xonly_Xtemp_Xtrials,tempCorr.rankCorr.Kendall.Xonly_Xtemp_Xtrials] = ...
    calculateRankCorrelations(tempCorr.templateRanks.Xonly_Xtrials,traces.X,find(tng.passed.AW.Xonly==1),p);


%% Calculate maximum bin correlations

% iscells
[tempCorr.maxBinCorr.Pearson.iscells_Atemp_Atrials,tempCorr.maxBinCorr.Spearman.iscells_Atemp_Atrials,tempCorr.maxBinCorr.Kendall.iscells_Atemp_Atrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.iscells_Atrials,traces.A,find(prop.iscell),p);
[tempCorr.maxBinCorr.Pearson.iscells_Atemp_Xtrials,tempCorr.maxBinCorr.Spearman.iscells_Atemp_Xtrials,tempCorr.maxBinCorr.Kendall.iscells_Atemp_Xtrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.iscells_Atrials,traces.X,find(prop.iscell),p);
[tempCorr.maxBinCorr.Pearson.iscells_Xtemp_Atrials,tempCorr.maxBinCorr.Spearman.iscells_Xtemp_Atrials,tempCorr.maxBinCorr.Kendall.iscells_Xtemp_Atrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.iscells_Xtrials,traces.A,find(prop.iscell),p);
[tempCorr.maxBinCorr.Pearson.iscells_Xtemp_Xtrials,tempCorr.maxBinCorr.Spearman.iscells_Xtemp_Xtrials,tempCorr.maxBinCorr.Kendall.iscells_Xtemp_Xtrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.iscells_Xtrials,traces.X,find(prop.iscell),p);

% A cells
[tempCorr.maxBinCorr.Pearson.A_Atemp_Atrials,tempCorr.maxBinCorr.Spearman.A_Atemp_Atrials,tempCorr.maxBinCorr.Kendall.A_Atemp_Atrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.A_Atrials,traces.A,find(tng.passed.AW.A==1),p);
[tempCorr.maxBinCorr.Pearson.A_Atemp_Xtrials,tempCorr.maxBinCorr.Spearman.A_Atemp_Xtrials,tempCorr.maxBinCorr.Kendall.A_Atemp_Xtrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.A_Atrials,traces.X,find(tng.passed.AW.A==1),p);
[tempCorr.maxBinCorr.Pearson.A_Xtemp_Atrials,tempCorr.maxBinCorr.Spearman.A_Xtemp_Atrials,tempCorr.maxBinCorr.Kendall.A_Xtemp_Atrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.A_Xtrials,traces.A,find(tng.passed.AW.A==1),p);
[tempCorr.maxBinCorr.Pearson.A_Xtemp_Xtrials,tempCorr.maxBinCorr.Spearman.A_Xtemp_Xtrials,tempCorr.maxBinCorr.Kendall.A_Xtemp_Xtrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.A_Xtrials,traces.X,find(tng.passed.AW.A==1),p);

% X cells
[tempCorr.maxBinCorr.Pearson.X_Atemp_Atrials,tempCorr.maxBinCorr.Spearman.X_Atemp_Atrials,tempCorr.maxBinCorr.Kendall.X_Atemp_Atrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.X_Atrials,traces.A,find(tng.passed.AW.X==1),p);
[tempCorr.maxBinCorr.Pearson.X_Atemp_Xtrials,tempCorr.maxBinCorr.Spearman.X_Atemp_Xtrials,tempCorr.maxBinCorr.Kendall.X_Atemp_Xtrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.X_Atrials,traces.X,find(tng.passed.AW.X==1),p);
[tempCorr.maxBinCorr.Pearson.X_Xtemp_Atrials,tempCorr.maxBinCorr.Spearman.X_Xtemp_Atrials,tempCorr.maxBinCorr.Kendall.X_Xtemp_Atrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.X_Xtrials,traces.A,find(tng.passed.AW.X==1),p);
[tempCorr.maxBinCorr.Pearson.X_Xtemp_Xtrials,tempCorr.maxBinCorr.Spearman.X_Xtemp_Xtrials,tempCorr.maxBinCorr.Kendall.X_Xtemp_Xtrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.X_Xtrials,traces.X,find(tng.passed.AW.X==1),p);

% Aonly cells
[tempCorr.maxBinCorr.Pearson.Aonly_Atemp_Atrials,tempCorr.maxBinCorr.Spearman.Aonly_Atemp_Atrials,tempCorr.maxBinCorr.Kendall.Aonly_Atemp_Atrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.Aonly_Atrials,traces.A,find(tng.passed.AW.Aonly==1),p);
[tempCorr.maxBinCorr.Pearson.Aonly_Atemp_Xtrials,tempCorr.maxBinCorr.Spearman.Aonly_Atemp_Xtrials,tempCorr.maxBinCorr.Kendall.Aonly_Atemp_Xtrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.Aonly_Atrials,traces.X,find(tng.passed.AW.Aonly==1),p);
[tempCorr.maxBinCorr.Pearson.Aonly_Xtemp_Atrials,tempCorr.maxBinCorr.Spearman.Aonly_Xtemp_Atrials,tempCorr.maxBinCorr.Kendall.Aonly_Xtemp_Atrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.Aonly_Xtrials,traces.A,find(tng.passed.AW.Aonly==1),p);
[tempCorr.maxBinCorr.Pearson.Aonly_Xtemp_Xtrials,tempCorr.maxBinCorr.Spearman.Aonly_Xtemp_Xtrials,tempCorr.maxBinCorr.Kendall.Aonly_Xtemp_Xtrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.Aonly_Xtrials,traces.X,find(tng.passed.AW.Aonly==1),p);

% Xonly cells
[tempCorr.maxBinCorr.Pearson.Xonly_Atemp_Atrials,tempCorr.maxBinCorr.Spearman.Xonly_Atemp_Atrials,tempCorr.maxBinCorr.Kendall.Xonly_Atemp_Atrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.Xonly_Atrials,traces.A,find(tng.passed.AW.Xonly==1),p);
[tempCorr.maxBinCorr.Pearson.Xonly_Atemp_Xtrials,tempCorr.maxBinCorr.Spearman.Xonly_Atemp_Xtrials,tempCorr.maxBinCorr.Kendall.Xonly_Atemp_Xtrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.Xonly_Atrials,traces.X,find(tng.passed.AW.Xonly==1),p);
[tempCorr.maxBinCorr.Pearson.Xonly_Xtemp_Atrials,tempCorr.maxBinCorr.Spearman.Xonly_Xtemp_Atrials,tempCorr.maxBinCorr.Kendall.Xonly_Xtemp_Atrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.Xonly_Xtrials,traces.A,find(tng.passed.AW.Xonly==1),p);
[tempCorr.maxBinCorr.Pearson.Xonly_Xtemp_Xtrials,tempCorr.maxBinCorr.Spearman.Xonly_Xtemp_Xtrials,tempCorr.maxBinCorr.Kendall.Xonly_Xtemp_Xtrials] = ...
    calculateMaxBinCorrelations(tempCorr.templateMaxBins.Xonly_Xtrials,traces.X,find(tng.passed.AW.Xonly==1),p);


%% Visualise what is actually correlated with each other

% figure;
% 
% % directly rank-correlating order of neurons
% subplot(2,1,1)
% hold on
% these_idcs = find(tng.passed.AW.Xonly==1);
% these_traces = traces.X(these_idcs,p.general.bins_analysisWindow,:);
% [~,temp1] = nanmax(these_traces(:,:,1),[],2);
% [~,temp2] = sort(temp1);
% these_ranks = these_idcs(temp2);
% plot(these_ranks,'r:')
% [~,temp1] = nanmax(these_traces(:,:,2),[],2);
% [~,temp2] = sort(temp1);
% these_ranks = these_idcs(temp2);
% plot(these_ranks,'b:')
% plot(tempCorr.templates.Xonly_Xtrials,'k-')
% xlabel('Xonly cells')
% ylabel('rank')
% title('rank-correlating order of neurons')
% 
% % rank-correlating maximal-rate bin of template to this trial
% subplot(2,1,2)
% hold on
% these_idcs = find(tng.passed.AW.Xonly==1);
% these_traces = traces.X(these_idcs,p.general.bins_analysisWindow,:);
% [~,these_bins] = nanmax(these_traces(:,:,1),[],2);
% plot(these_bins,'r:')
% [~,these_bins] = nanmax(these_traces(:,:,2),[],2);
% plot(these_bins,'b:')
% plot(tempCorr.templates.Xonly_Xtrials_bins,'k-')
% xlabel('Xonly cells')
% ylabel('bin with highest activity')
% title('(rank-)correlating maximal-rate bins')


