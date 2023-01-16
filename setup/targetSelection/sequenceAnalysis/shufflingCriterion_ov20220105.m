function shuffling = shufflingCriterion_ov20220105(data,iscell,p)
%data = traces.traces_A(:,prop.analysisWindow,:);

%% Preparations

% only do it for iscells
compact_data = data(find(iscell),:,:);

% get basic properties
numRois = size(data,1);
numCells = size(compact_data,1);
numFrames = size(data,2);
numTrials = size(data,3);

% bring data in correct shape
data_cell = num2cell(compact_data,[1 2]);


%% Core

% generate shuffled distribution
rng(p.sca.rgnSeed);
these_maxRates_shuffled = zeros(numCells,p.sca.numShuffles);
for i=1:p.sca.numShuffles

    % circular permutation, each trial is shuffled differently
    this_shuffle = cellfun(@(x) circshift(x,randi([-round(numFrames/2),round(numFrames/2)],1,1),2),data_cell,'UniformOutput',false);

    % averaging over the shuffled intervals
    this_average = nanmean(cat(3,this_shuffle{:}),3);

    % saving the maximum mean firing rate for each cell
    these_maxRates_shuffled(:,i) = nanmax(this_average,[],2);
end

shuffling.maxRates_shuffled = nan(numRois,p.sca.numShuffles);
shuffling.maxRates_shuffled(find(iscell),:) = these_maxRates_shuffled;

% test maximum mean firing rate against shuffled distrinbution
maxRatePrctiles_shuffled = prctile(shuffling.maxRates_shuffled,p.sca.shufflingThreshold,2);
shuffling.maxRates_true =  nanmax(nanmean(data,3),[],2);
shuffling.significant = nan(numRois,1);
shuffling.significant(find(iscell)) = maxRatePrctiles_shuffled(find(iscell)) < shuffling.maxRates_true(find(iscell));

end