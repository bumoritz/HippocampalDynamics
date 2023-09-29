function [these_maxBinCorrs_Pearson,these_maxBinCorrs_Spearman,these_maxBinCorrs_Kendall] = calculateMaxBinCorrelations(these_templateMaxBins,these_traces,these_idcs,p);
% these_templateMaxBins = tempCorr.templateMaxBins.iscells_Atrials; these_traces = traces.A; these_idcs = find(prop.iscell);

% extract data
this_data = these_traces(these_idcs,p.general.bins_analysisWindow,:);

% initialise
these_maxBinCorrs_Spearman = nan(size(this_data,3),1);
these_maxBinCorrs_Kendall = nan(size(this_data,3),1);
these_maxBinCorrs_Pearson = nan(size(this_data,3),1);
for i=1:size(this_data,3)
    
    % find bins with maximum activity
    [~,these_maxBins] = nanmax(this_data(:,:,i),[],2);

    % calculate rank correlations
    [these_maxBinCorrs_Spearman(i),~] = corr(these_templateMaxBins,these_maxBins,'Type','Spearman','Rows','complete');
    [these_maxBinCorrs_Kendall(i),~] = corr(these_templateMaxBins,these_maxBins,'Type','Kendall','Rows','complete');    
    [these_maxBinCorrs_Pearson(i),~] = corr(these_templateMaxBins,these_maxBins,'Type','Pearson','Rows','complete');    
end

end
