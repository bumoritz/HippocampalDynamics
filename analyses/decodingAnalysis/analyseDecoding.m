function [analysis] = analyseDecoding(dec,this_set,these_trials,p,onlySequenceCells)
% dec = dec_all; this_set = 'trainingSet'; these_trials = 1:dec_all.cv.numTrials_trainingSet; onlySequenceCells = true;
%old: dec = dec_all; this_set = 'trainingSet'; these_trials = 1:dec_all.cv.numTrials_trainingSet; 
%old: dec = dec_all; this_set = 'completeSet'; these_trials = dec_all.cv.trialsForTestSet_incorrM.combined; 

if nargin < 5
    onlySequenceCells = false;
end

if onlySequenceCells
    results = dec.([this_set,'_seq']);
else
    results = dec.(this_set);
end

numTrials = length(these_trials);
numBins = dec.bins.numBins;


%% Time decoding error

analysis.timeDecodingError_s = nan(numTrials,numBins);
analysis.timeDecodingError_withinCat_s = nan(numTrials,numBins);
for t=1:numTrials
    this_trial = these_trials(t);
    for b=1:numBins
        this_label_true = results.Y(this_trial,b);
        this_label_predicted = round(results.Ypred(this_trial,b));
        this_t_true = dec.classes.t(this_label_true);
        this_t_predicted = dec.classes.t(this_label_predicted);
        
        analysis.timeDecodingError_s(t,b) = abs(this_t_true-this_t_predicted);
        if strcmp(dec.classes.supClass{this_label_true},dec.classes.supClass{this_label_predicted})
            analysis.timeDecodingError_withinCat_s(t,b) = abs(this_t_true-this_t_predicted);
        end
    end
end


%% Odour decoding

analysis.typeDecodingCorrect = nan(numTrials,numBins);
analysis.typeDecodingCorrect_withinWindow = nan(numTrials,numBins);
for t=1:numTrials
    this_trial = these_trials(t);
    for b=1:numBins
        this_label_true = results.Y(this_trial,b);
        this_label_predicted = round(results.Ypred(this_trial,b));
        this_class_true = dec.classes.supClass{this_label_true};
        this_class_predicted = dec.classes.supClass{this_label_predicted};
        
        if length(dec.classes.supClass{this_label_true})==length(dec.classes.supClass{this_label_predicted})
            if strcmp(dec.classes.supClass{this_label_true},dec.classes.supClass{this_label_predicted})
                analysis.typeDecodingCorrect(t,b) = 1;
                analysis.typeDecodingCorrect_withinWindow(t,b) = 1;
            else
                analysis.typeDecodingCorrect(t,b) = 0;
                analysis.typeDecodingCorrect_withinWindow(t,b) = 0;
            end
        else
            analysis.typeDecodingCorrect(t,b) = 0;
        end
    end
end


%% Create Ypred_within and posterior_within

Ypred_within = nan(size(results.Ypred));
posterior_within = nan(size(results.posterior));
for i=1:size(results.Ypred,1)
    for j=1:size(results.Ypred,2)
        if dec.classes.supClassIdx(results.Y(i,j)) == dec.classes.supClassIdx(round(results.Ypred(i,j)))
            Ypred_within(i,j) = results.Ypred(i,j);
        end
        posterior_within(find(dec.classes.supClassIdx==dec.classes.supClassIdx(results.Y(i,j))),j,i) = ...
            results.posterior(find(dec.classes.supClassIdx==dec.classes.supClassIdx(results.Y(i,j))),j,i);
    end
end


%% Calculate measures for first and second halves independently

% try
%     this_first = [1,length(dec.supClasses.bins{1})];
%     this_second = [length(dec.supClasses.bins{1})+1,length(dec.supClasses.bins{1})+length(dec.supClasses.bins{3})];
% 
%     analysis.first.slope_discr = nan(numTrials,1);
%     analysis.second.slope_discr = nan(numTrials,1);
%     for t=1:numTrials
%         this_trial = these_trials(t);
%         this_mdl = fitlm((this_first(1):this_first(2))',round(Ypred_within(this_trial,this_first(1):this_first(2))'),'Intercept',false);
%         analysis.first.slope_discr(t) = this_mdl.Coefficients.Estimate;
%     %     this_mdl = fitlm((this_second(1):this_second(2))'-this_first(2),Ypred_within(this_trial,this_second(1):this_second(2))','Intercept',false);
%     % 	analysis.second.slope_discr(t) = this_mdl.Coefficients.Estimate;
%     end
% catch
% end



