function [decs] = createDecoderStruct(dec,this_data,these_trialsForSet,these_conditions)
%these_trialsForSet = dec.cv.trialsForTrainingSet; these_conditions = {'A_H_O','X_H_O','A_CR_O','X_CR_O'};

%% Preparations

decs.numCells = size(this_data,1);
decs.numBins = size(this_data,2);
decs.numTrials = size(this_data,3);
decs.numSamples = decs.numBins*decs.numTrials;
decs.numClasses = dec.classes.numClasses;


%% Create struct

decs.X = nan(decs.numSamples,decs.numCells);
decs.Y = nan(decs.numSamples,1);
decs.trialNumber_abs = nan(decs.numSamples,1);
decs.trialNumber_rel = nan(decs.numSamples,1);
decs.trialType = nan(decs.numSamples,1);
decs.bin = nan(decs.numSamples,1);
n=0;
for k=1:decs.numTrials
    for j=1:decs.numBins
        n=n+1;
        
        decs.X(n,:) = this_data(:,j,k)';
        
        if ismember(these_trialsForSet.combined(k),these_trialsForSet.(these_conditions{1}))
            decs.trialType(n) = 1;
            decs.Y(n) = intersect( find(dec.classes.bin==dec.bins.bins(j)), find(strcmp(dec.classes.supClass,'A')|strcmp(dec.classes.supClass,'AB')));
        elseif ismember(these_trialsForSet.combined(k),these_trialsForSet.(these_conditions{2}))
            decs.trialType(n) = 2;
            decs.Y(n) = intersect( find(dec.classes.bin==dec.bins.bins(j)), find(strcmp(dec.classes.supClass,'X')|strcmp(dec.classes.supClass,'XY')));
        elseif ismember(these_trialsForSet.combined(k),these_trialsForSet.(these_conditions{3}))
            decs.trialType(n) = 3;
            decs.Y(n) = intersect( find(dec.classes.bin==dec.bins.bins(j)), find(strcmp(dec.classes.supClass,'A')|strcmp(dec.classes.supClass,'AY')));
        elseif ismember(these_trialsForSet.combined(k),these_trialsForSet.(these_conditions{4}))
            decs.trialType(n) = 4;
            decs.Y(n) = intersect( find(dec.classes.bin==dec.bins.bins(j)), find(strcmp(dec.classes.supClass,'X')|strcmp(dec.classes.supClass,'XB')));
        else
            warning('Error. Check what happened.')
        end
        
        decs.trialNumber_abs(n) = these_trialsForSet.combined(k);
        decs.trialNumber_rel(n) = k;
        decs.bin(n) = j;
    end
end

end