function [these_templateRanks,these_templateMaxBins] = createTemplate(these_avgTraces,these_idcs,p)
% these_avgTraces = avgTraces.A; these_idcs = find(prop.iscell);

% extract data
this_data = these_avgTraces(these_idcs,p.general.bins_analysisWindow);

% find bins with maximum activity
[~,these_templateMaxBins] = nanmax(this_data,[],2);

% sort cells by bin with maximum activity
[temp1,temp2] = sort(these_templateMaxBins);

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

% assign sequence template
these_templateRanks = these_idcs(temp3);

end