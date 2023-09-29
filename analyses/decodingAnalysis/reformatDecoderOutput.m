function [decs] = reformatDecoderOutput(decs)
%decs = decs_train; decs = decs_complete; 

decs.posterior_cbt = nan(decs.numClasses,decs.numBins,decs.numTrials);
decs.Y_tb = nan(decs.numTrials,decs.numBins);
decs.Ypred_tb = nan(decs.numTrials,decs.numBins);
decs.trialType_t = nan(decs.numTrials,1);
n=0;
for k=1:decs.numTrials
    for j=1:decs.numBins
        n=n+1;
        decs.posterior_cbt(:,j,k) = decs.posterior(n,:)';
        decs.Y_tb(k,j) = decs.Y(n);
        decs.Ypred_tb(k,j) = decs.Ypred(n);
    end
    decs.trialType_t(k) = decs.trialType(n);
end
% 
% these_conditions = unique(decs.trialType_t);
% for i=1:length(these_conditions)
%     decs.(['posterior_type',num2str(these_conditions(i))]) = decs.posterior_cbt(:,:,find(decs.trialType_t==these_conditions(i)));
%     decs.(['avgPosterior_type',num2str(these_conditions(i))]) = nanmean(decs.(['posterior_type',num2str(these_conditions(i))]),3);
%     decs.(['numTrials_type',num2str(these_conditions(i))]) = size(decs.(['posterior_type',num2str(these_conditions(i))]),3);
% end

end