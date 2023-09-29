function shuffling = shufflingCriterion(prop,this_nf_binned,these_trials,these_traces,this_window,p,this_binRange)
% this_nf_binned = nf_binned; these_trials = trials_all.stimuli.A; these_traces = traces_all.A; this_window = p.general.bins_analysisWindow;
%this_nf_binned = nf_binned; these_trials = trials.stimuli.A; these_traces = traces.A; this_window = p.general.bins_analysisWindow;

%% Preparations

% strip away all non-iscells
these_traces_comp = these_traces(find(prop.iscell),this_window,:);
nf_binned_comp = this_nf_binned(find(prop.iscell),this_binRange(1):this_binRange(2));

% get basic properties
numRois = size(these_traces,1);
numCells = size(these_traces_comp,1);
numTrials = size(these_traces_comp,3);
numFramesTotal = size(nf_binned_comp,2);
numFrames = size(these_traces_comp,2);


%% Core

% generate shuffled distribution
rng(p.general.rgnSeed);
these_maxRates_shuffled = zeros(numCells,p.tng.numShuffles);
try
    parfor i=1:p.tng.numShuffles

        % circular permutation and splitting into trials
        temp = circshift(nf_binned_comp,randi([-round(numFramesTotal/2),round(numFramesTotal/2)],1,1),2);
        this_shuffle = reshape(temp(:,prop.trial_frames_binned(this_window,these_trials)-this_binRange(1)+1),[],numFrames,numTrials);

        % averaging over the shuffled intervals
        this_average = nanmean(this_shuffle,3);

        % saving the maximum mean firing rate for each cell
        these_maxRates_shuffled(:,i) = nanmax(this_average,[],2);
    end
catch
    try
        parfor i=1:p.tng.numShuffles

            % circular permutation and splitting into trials
            temp = circshift(nf_binned_comp,randi([-round(numFramesTotal/2),round(numFramesTotal/2)],1,1),2);
            this_shuffle = reshape(temp(:,prop.trial_frames_binned(this_window,these_trials(1:end-1))-this_binRange(1)+1),[],numFrames,numTrials-1);

            % averaging over the shuffled intervals
            this_average = nanmean(this_shuffle,3);

            % saving the maximum mean firing rate for each cell
            these_maxRates_shuffled(:,i) = nanmax(this_average,[],2);
        end
    catch
        parfor i=1:p.tng.numShuffles

            % circular permutation and splitting into trials
            temp = circshift(nf_binned_comp,randi([-round(numFramesTotal/2),round(numFramesTotal/2)],1,1),2);
            temp2 = these_trials(1:end-1);
            temp2 = temp2(temp2<size(prop.trial_frames_binned,2)-1);
            this_shuffle = reshape(temp(:,prop.trial_frames_binned(this_window,temp2)-this_binRange(1)+1),[],numFrames,length(temp2));

            % averaging over the shuffled intervals
            this_average = nanmean(this_shuffle,3);

            % saving the maximum mean firing rate for each cell
            these_maxRates_shuffled(:,i) = nanmax(this_average,[],2);
        end
    end
end

shuffling.maxRates_shuffled = nan(numRois,p.tng.numShuffles);
shuffling.maxRates_shuffled(find(prop.iscell),:) = these_maxRates_shuffled;

% test maximum mean firing rate against shuffled distrinbution
maxRatePrctiles_shuffled = prctile(shuffling.maxRates_shuffled,p.tng.shufflingThreshold,2);
shuffling.maxRates_true = nan(numRois,1);
shuffling.maxRates_true(find(prop.iscell)) = nanmax(nanmean(these_traces_comp,3),[],2);
shuffling.significant_raw = nan(numRois,1);
shuffling.significant_raw(find(prop.iscell)) = maxRatePrctiles_shuffled(find(prop.iscell)) < shuffling.maxRates_true(find(prop.iscell));
shuffling.pvals = nan(numRois,1);
shuffling.pvals(find(prop.iscell)) = nansum(shuffling.maxRates_shuffled(find(prop.iscell),:)>=shuffling.maxRates_true(find(prop.iscell)),2)./p.tng.numShuffles;

if p.tng.pfdrCorrectionForShuffle
    [~,shuffling.qvals,~] = mafdr(shuffling.pvals);
    shuffling.significant = nan(numRois,1);
    shuffling.significant(find(prop.iscell)) = double(shuffling.qvals(find(prop.iscell))<p.tng.qThreshold);
else
    shuffling.significant = shuffling.significant_raw;
end

shuffling = orderfields(shuffling);
end