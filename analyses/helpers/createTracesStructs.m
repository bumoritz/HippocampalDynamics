function [traces,avgTraces,normAvgTraces] = createTracesStructs(nft_binned,these_trials);
% these_trials = trials_all;

%% traces

% traces for stimuli
try
    temp = fields(these_trials.stimuli);
    for i=1:length(temp)
        traces.(temp{i}) = nft_binned(:,:,these_trials.stimuli.(temp{i}));
    end
catch
end

% traces for outcome
try
    temp = fields(these_trials.outcome);
    for i=1:length(temp)
        traces.(temp{i}) = nft_binned(:,:,these_trials.outcome.(temp{i}));
    end
catch
end

% traces for balanced
try
    temp = fields(these_trials.balanced);
    for i=1:length(temp)
        traces.(temp{i}) = nft_binned(:,:,these_trials.balanced.(temp{i}));
    end
catch
end

% try
%     % traces for stimuli_var0
%     temp = fields(these_trials.stimuli_var0);
%     for i=1:length(temp)
%         traces.([temp{i},'_var0']) = nft_binned(:,:,these_trials.stimuli_var0.(temp{i}));
%     end
% 
%     % traces for outcome_var0
%     temp = fields(these_trials.outcome_var0);
%     for i=1:length(temp)
%         traces.([temp{i},'_var0']) = nft_binned(:,:,these_trials.outcome_var0.(temp{i}));
%     end
% 
%     % traces for stimuli_var1
%     temp = fields(these_trials.stimuli_var1);
%     for i=1:length(temp)
%         traces.([temp{i},'_var1']) = nft_binned(:,:,these_trials.stimuli_var1.(temp{i}));
%     end
% 
%     % traces for outcome_var1
%     temp = fields(these_trials.outcome_var1);
%     for i=1:length(temp)
%         traces.([temp{i},'_var1']) = nft_binned(:,:,these_trials.outcome_var1.(temp{i}));
%     end
% catch
% end


%% avgTraces

% avgTraces
temp = fields(traces);
for i=1:length(temp)
    avgTraces.(temp{i}) = nanmean(traces.(temp{i}),3);
end


%% normAvgTraces

% normAvgTraces.firstOdour (e.g. for first odour selectivity index)
try
    temp = [avgTraces.A(:);avgTraces.X(:)];
    this_lower = nanmin(temp(:));
    this_upper = nanmax(temp(:));
    normAvgTraces.firstOdour.A = (avgTraces.A-this_lower) ./ (this_upper-this_lower);
    normAvgTraces.firstOdour.X = (avgTraces.X-this_lower) ./ (this_upper-this_lower);
catch
end

% normAvgTraces.firstOdour for _cv case
try
    temp = [avgTraces.A_correct_train(:);avgTraces.X_correct_train(:)];
    this_lower = nanmin(temp(:));
    this_upper = nanmax(temp(:));
    normAvgTraces.firstOdour.A_correct_train = (avgTraces.A_correct_train-this_lower) ./ (this_upper-this_lower);
    normAvgTraces.firstOdour.X_correct_train = (avgTraces.X_correct_train-this_lower) ./ (this_upper-this_lower);
catch
end

% normAvgTraces.secondOdour
try
    temp = [avgTraces.B(:);avgTraces.Y(:)];
    this_lower = nanmin(temp(:));
    this_upper = nanmax(temp(:));
    normAvgTraces.secondOdour.B = (avgTraces.B-this_lower) ./ (this_upper-this_lower);
    normAvgTraces.secondOdour.Y = (avgTraces.Y-this_lower) ./ (this_upper-this_lower);
catch
end

% normAvgTraces.trialType
try
    temp = [avgTraces.AB(:);avgTraces.XB(:);avgTraces.AY(:);avgTraces.XY(:)];
    this_lower = nanmin(temp(:));
    this_upper = nanmax(temp(:));
    normAvgTraces.trialType.AB = (avgTraces.AB-this_lower) ./ (this_upper-this_lower);
    normAvgTraces.trialType.XB = (avgTraces.XB-this_lower) ./ (this_upper-this_lower);
    normAvgTraces.trialType.XY = (avgTraces.XY-this_lower) ./ (this_upper-this_lower);
    normAvgTraces.trialType.AY = (avgTraces.AY-this_lower) ./ (this_upper-this_lower);
catch
end

% normAvgTraces.reward
try
    temp = [avgTraces.rew(:);avgTraces.norew(:)];
    this_lower = nanmin(temp(:));
    this_upper = nanmax(temp(:));
    normAvgTraces.reward.rew = (avgTraces.rew-this_lower) ./ (this_upper-this_lower);
    normAvgTraces.reward.norew = (avgTraces.norew-this_lower) ./ (this_upper-this_lower);
catch
end


%% Return

traces = orderfields(traces);
avgTraces = orderfields(avgTraces);
normAvgTraces = orderfields(normAvgTraces);
end

