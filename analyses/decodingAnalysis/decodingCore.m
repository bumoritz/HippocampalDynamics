function [dec,dec_mdl] = decodingCore(dec,this_nft_binned,prop,p,normalise_data)
% dec = dec_all; this_nft_binned = nft_binned;

if nargin < 5
    normalise_data = false;
end


%% Create decoder struct

% training set
these_conditions = fields(rmfield(dec.cv.trialsForTrainingSet,'combined'));
this_data = this_nft_binned(find(prop.iscell==1),dec.bins.bins,dec.cv.trialsForTrainingSet.combined);
decs_train = createDecoderStruct(dec,this_data,dec.cv.trialsForTrainingSet,these_conditions,normalise_data);

% complete set
these_conditions = fields(rmfield(dec.cv.trialsForCompleteSet,'combined'));
this_data = this_nft_binned(find(prop.iscell==1),dec.bins.bins,dec.cv.trialsForCompleteSet.combined);
decs_complete = createDecoderStruct(dec,this_data,dec.cv.trialsForCompleteSet,these_conditions,normalise_data);


%% Train decoder and predict

decs_train.Ypred = nan(dec.bins.numBins*dec.cv.numTrials_trainingSet,1);
decs_train.posterior = nan(dec.bins.numBins*dec.cv.numTrials_trainingSet,dec.classes.numClasses);
decs_complete.Ypred = nan(dec.bins.numBins*dec.cv.numTrials_completeSet,p.dec.cvFold);
decs_complete.posterior = nan(dec.bins.numBins*dec.cv.numTrials_completeSet,dec.classes.numClasses,p.dec.cvFold);
dec_mdl = cell(1,p.dec.cvFold);
for k=1:p.dec.cvFold
    
    % identify training and test bins
    these_trainingBins = find(dec.cv.trainingSet_foldByBin~=k);
    these_testBins = find(dec.cv.trainingSet_foldByBin==k);
    
    % train decoder
    dec_mdl{k} = fitcnb(decs_train.X(these_trainingBins,:),decs_train.Y(these_trainingBins),'Prior',dec.classes.theoreticalPriors);
    
    % make predictions - trainingSet
    [this_Ypred_train,this_posterior_train,~] = predict(dec_mdl{k},decs_train.X(these_testBins,:));   
    decs_train.Ypred(these_testBins) = this_Ypred_train;
    decs_train.posterior(these_testBins,:) = this_posterior_train;

    % make predictions - completeSet
    [this_Ypred,this_posterior,~] = predict(dec_mdl{k},decs_complete.X);   
    decs_complete.Ypred(:,k) = this_Ypred;
    decs_complete.posterior(:,:,k) = this_posterior;
end
decs_train = reformatDecoderOutput(decs_train);
decs_complete.Ypred = nanmean(decs_complete.Ypred,2);
decs_complete.posterior = nanmean(decs_complete.posterior,3);
decs_complete = reformatDecoderOutput(decs_complete);


%% Add decoder output to dec struct

dec.trainingSet.trialType = decs_train.trialType_t;
dec.trainingSet.Y = decs_train.Y_tb;
dec.trainingSet.Ypred = decs_train.Ypred_tb;
dec.trainingSet.posterior = decs_train.posterior_cbt;

dec.completeSet.trialType = decs_complete.trialType_t;
dec.completeSet.Y = decs_complete.Y_tb;
dec.completeSet.Ypred = decs_complete.Ypred_tb;
dec.completeSet.posterior = decs_complete.posterior_cbt;


end
