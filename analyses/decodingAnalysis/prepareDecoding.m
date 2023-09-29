function [dec] = prepareDecoding(these_trials,p)
% these_trials = trials_all; 

%% Assign classes and supClasses

% create supClass fields
dec.supClasses.labels = {'A','X'};
dec.supClasses.labels_extd = {'A_{0s-5s}','X_{0s-5s}'};
dec.supClasses.numSupClasses = length(dec.supClasses.labels);
dec.supClasses.bins = {16:40,16:40};

% create bins fields
dec.bins.bins = min([dec.supClasses.bins{:}]):max([dec.supClasses.bins{:}]);
dec.bins.numBins = length(dec.bins.bins);

% create classes fields
dec.classes.supClassIdx = [];
for i=1:dec.supClasses.numSupClasses
    dec.classes.supClassIdx = [dec.classes.supClassIdx, i*ones(1,length(dec.supClasses.bins{i}))];
end
dec.classes.numClasses = length(dec.classes.supClassIdx);
dec.classes.supClass = {dec.supClasses.labels{dec.classes.supClassIdx}};
dec.classes.bin = [dec.supClasses.bins{:}];
dec.classes.t = p.general.t_binned(dec.classes.bin);
dec.classes.labels = {};
dec.classes.labels_s = {};
for i=1:dec.classes.numClasses
    dec.classes.labels{i} = [dec.classes.supClass{i},'-',num2str(dec.classes.bin(i))];
    dec.classes.labels_s{i} = [dec.classes.supClass{i},'-',num2str(dec.classes.t(i),2),'s'];
end

% calculate theoretical priors
temp = [];
for i=1:dec.supClasses.numSupClasses
    if strcmp(dec.supClasses.labels{i},'A') || strcmp(dec.supClasses.labels{i},'X') || strcmp(dec.supClasses.labels{i},'B') || strcmp(dec.supClasses.labels{i},'Y')
        temp = [temp, 0.5*ones(1,length(dec.supClasses.bins{i}))];
    elseif strcmp(dec.supClasses.labels{i},'AB') || strcmp(dec.supClasses.labels{i},'AY') || strcmp(dec.supClasses.labels{i},'XY') || strcmp(dec.supClasses.labels{i},'XB')
        temp = [temp, 0.25*ones(1,length(dec.supClasses.bins{i}))];
    end
end
dec.classes.theoreticalPriors = temp/sum(temp);


%% Splitting into sets

% assign trials for complete set
these_trialTypes = {'AB','XY','AY','XB'};
dec.cv.trialsForCompleteSet.combined = [];
for i=1:length(these_trialTypes)
    dec.cv.trialsForCompleteSet_raw.(these_trialTypes{i}) = these_trials.stimuli.(these_trialTypes{i});
    dec.cv.trialsForCompleteSet.combined = sort([dec.cv.trialsForCompleteSet.combined, dec.cv.trialsForCompleteSet_raw.(these_trialTypes{i})]);
end
dec.cv.numTrials_completeSet = length(dec.cv.trialsForCompleteSet.combined);
dec.cv.trialsForCompleteSet.A = [dec.cv.trialsForCompleteSet_raw.AB,dec.cv.trialsForCompleteSet_raw.AY];
dec.cv.trialsForCompleteSet.X = [dec.cv.trialsForCompleteSet_raw.XY,dec.cv.trialsForCompleteSet_raw.XB];

if p.dec.trainOnCompleteSet
    dec.cv.trialsForTrainingSet = dec.cv.trialsForCompleteSet;
    dec.cv.numTrials_trainingSet = dec.cv.numTrials_completeSet;
else
%     % assign trials for training set
%     these_trialTypes = {'A_H_O','A_CR_O','X_H_O','X_CR_O'};
%     dec.cv.trialsForTrainingSet.combined = [];
%     for i=1:length(these_trialTypes)
%         dec.cv.trialsForTrainingSet.(these_trialTypes{i}) = these_trials.balanced.(these_trialTypes{i});
%         dec.cv.trialsForTrainingSet.combined = sort([dec.cv.trialsForTrainingSet.combined, dec.cv.trialsForTrainingSet.(these_trialTypes{i})]);
%     end
%     dec.cv.numTrials_trainingSet = length(dec.cv.trialsForTrainingSet.combined);
% 
%     % assign trials for other sets (corr_M)
%     these_trialTypes = {'A_H_M','A_CR_M','X_H_M','X_CR_M'};
%     dec.cv.trialsForTestSet_corrM.combined = [];
%     for i=1:length(these_trialTypes)
%         dec.cv.trialsForTestSet_corrM.(these_trialTypes{i}) = these_trials.balanced.(these_trialTypes{i});
%         dec.cv.trialsForTestSet_corrM.combined = sort([dec.cv.trialsForTestSet_corrM.combined, dec.cv.trialsForTestSet_corrM.(these_trialTypes{i})]);
%     end
%     dec.cv.numTrials_testSet_corrM = length(dec.cv.trialsForTestSet_corrM.combined);
% 
%     % assign trials for other sets (incorr_M)
%     these_trialTypes = {'A_M_M','A_FA_M','X_M_M','X_FA_M'};
%     dec.cv.trialsForTestSet_incorrM.combined = [];
%     for i=1:length(these_trialTypes)
%         dec.cv.trialsForTestSet_incorrM.(these_trialTypes{i}) = these_trials.balanced.(these_trialTypes{i});
%         dec.cv.trialsForTestSet_incorrM.combined = sort([dec.cv.trialsForTestSet_incorrM.combined, dec.cv.trialsForTestSet_incorrM.(these_trialTypes{i})]);
%     end
%     dec.cv.numTrials_testSet_incorrM = length(dec.cv.trialsForTestSet_incorrM.combined);
end


%% Assigning folds to trials of training set

this_numTrialsPerSet = floor(dec.cv.numTrials_trainingSet/p.dec.cvFold);
this_pool = 1:dec.cv.numTrials_trainingSet;
dec.cv.trainingSet_foldByTrial = ones(dec.cv.numTrials_trainingSet,1);
for i=1:p.dec.cvFold
    temp = datasample(this_pool,this_numTrialsPerSet,'Replace',false);
    dec.cv.trainingSet_foldByTrial(temp) = i;
    this_pool = setdiff(this_pool,temp);
end

dec.cv.trainingSet_foldByBin = ones(dec.bins.numBins*length(dec.cv.trialsForTrainingSet.combined),1);
n=0;
for k=1:length(dec.cv.trialsForTrainingSet.combined)
    for j=1:dec.bins.numBins
        n=n+1;
        dec.cv.trainingSet_foldByBin(n) = dec.cv.trainingSet_foldByTrial(k);
    end
end

end
