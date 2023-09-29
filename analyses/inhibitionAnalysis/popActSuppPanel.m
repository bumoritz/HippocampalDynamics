function popActSuppPanel(these_avgTraces,these_idcs_unsorted,these_idcs_pos,these_idcs_zero,these_idcs_neg,p,info,plt)
%these_avgTraces = avgTraces_all.X; these_idcs_unsorted = find(iscell==1); this_testGroup = 'X_{sens}';


% these_avgTraces = nanmean(traces_all.X(:,:,1:20),3);


%% Preparations

% 
this_baselineWindow = 1:10;
this_sortingWindow = 16:20;
this_plottingWindow = 1:25;


%% Normalise traces

this_baseline_mean = nanmean(these_avgTraces(:,this_baselineWindow),2);
this_baseline_std = nanstd(these_avgTraces(:,this_baselineWindow),[],2);

these_normAvgTraces = (these_avgTraces - this_baseline_mean) ./ this_baseline_std;


%% Sort traces

these_tracesForSorting = these_normAvgTraces; % these_avgTraces or these_normAvgTraces

if false %~isempty(these_idcs_pos) && ~isempty(these_idcs_neg)
    these_idcs_other = setdiff(setdiff(these_idcs_unsorted,these_idcs_pos),these_idcs_neg);
    [~,temp] = sort(nanmean(these_normAvgTraces(these_idcs_pos,this_sortingWindow),2));
    these_idcs_pos = these_idcs_pos(flip(temp));
    [~,temp] = sort(nanmean(these_normAvgTraces(these_idcs_neg,this_sortingWindow),2));
    these_idcs_neg = these_idcs_neg(flip(temp));
    [~,temp] = sort(nanmean(these_normAvgTraces(these_idcs_other,this_sortingWindow),2));
    these_idcs_other = these_idcs_other(flip(temp));
    these_idcs_sorted = [these_idcs_pos;these_idcs_other;these_idcs_neg];
    
elseif true
    these_idcs_other = setdiff(setdiff(these_idcs_unsorted,these_idcs_pos),these_idcs_neg);
%     [~,temp] = sort(nanmean(these_normAvgTraces(these_idcs_pos,this_sortingWindow),2));
%     these_idcs_pos = these_idcs_pos(flip(temp));
%     [~,temp] = sort(nanmean(these_normAvgTraces(these_idcs_zero,this_sortingWindow),2));
%     these_idcs_zero = these_idcs_zero(flip(temp));
%     [~,temp] = sort(nanmean(these_normAvgTraces(these_idcs_neg,this_sortingWindow),2));
%     these_idcs_neg = these_idcs_neg(flip(temp));
%     [~,temp] = sort(nanmean(these_normAvgTraces(these_idcs_other,this_sortingWindow),2));
%     these_idcs_other = these_idcs_other(flip(temp));
    these_idcs_sorted = [these_idcs_pos;these_idcs_zero;these_idcs_neg;these_idcs_other];
    
else
    [~,temp] = sort(nanmean(these_normAvgTraces(these_idcs_unsorted,this_sortingWindow),2));
    these_idcs_sorted = these_idcs_unsorted(flip(temp));
end

                
%% plot              

plt.clim = [-15,15]; plt.colormap = redblue;
heatMap_task(these_normAvgTraces(these_idcs_sorted,this_plottingWindow),NaN,[],p,info,plt);

end

