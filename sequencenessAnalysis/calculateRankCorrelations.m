function [these_rankCorrs_Pearson,these_rankCorrs_Spearman,these_rankCorrs_Kendall] = calculateRankCorrelations(these_templateRanks,these_traces,these_idcs,p);
% these_templateRanks = tempCorr.templateRanks.iscells_Atrials; these_traces = traces.A; these_idcs = find(prop.iscell);

% extract data
this_data = these_traces(these_idcs,p.general.bins_analysisWindow,:);

% initialise
these_rankCorrs_Spearman = nan(size(this_data,3),1);
these_rankCorrs_Kendall = nan(size(this_data,3),1);
these_rankCorrs_Pearson = nan(size(this_data,3),1);
for i=1:size(this_data,3)
    
    % find bins with maximum activity
    [~,these_maxBins] = nanmax(this_data(:,:,i),[],2);
    
    % sort cells by bin with maximum activity
    [temp1,temp2] = sort(these_maxBins);
    
    % shuffle order of cells that have maximum activity in same bin
    if p.sqn.shuffleOrderForSameBins
        temp3 = temp2;
        for j=1:length(p.general.bins_analysisWindow)
            this_pool = find(temp1==j);
            temp3(this_pool) = temp2(this_pool(randperm(length(this_pool))));
        end  
    else
        temp3 = temp2;
    end
    
    % get ranks for this trial
    these_ranks = these_idcs(temp3);
    
    % calculate rank correlations
    [these_rankCorrs_Spearman(i),~] = corr(these_templateRanks,these_ranks,'Type','Spearman','Rows','complete');
    [these_rankCorrs_Kendall(i),~] = corr(these_templateRanks,these_ranks,'Type','Kendall','Rows','complete');    
    [these_rankCorrs_Pearson(i),~] = corr(these_templateRanks,these_ranks,'Type','Pearson','Rows','complete');    
end

end