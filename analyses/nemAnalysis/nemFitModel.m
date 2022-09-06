function [mdl,coefs,Ypred,R2_train,R2_test,R2_all] = nemFitModel(nem,nemd,p,iscell,these_cells,X,exclude,Y)
        
% set GLM options
options = glmnetSet;
options.alpha = p.nem.alpha;
if nargin < 7
    options.exclude = [];
else
    options.exclude = exclude;
end

% extract Y
if nargin < 8
    Y = nemd.data;
end

% prepare outputs
mdl = {};
coefs = nan(length(iscell),nem.numPredictors+1);
Ypred = nan(size(Y));
R2_train = nan(length(iscell),1);
R2_test = nan(length(iscell),1);
R2_all = nan(length(iscell),1);

setsByBins = nem.setsByBins;
these_trainingBins = find(nem.setsByBins~=0);
these_testBins = find(nem.setsByBins==0);
parfor i=these_cells
    if iscell(i)
        try
            this_Y = Y(:,i);
            mdl{i} = cvglmnet(X(these_trainingBins,:),this_Y(these_trainingBins),p.nem.family,options,p.nem.cvLoss,p.nem.cvFolds,nonzeros(setsByBins),false);
            % cvglmnetPlot(nem.mdl{i}); % middle line is log(cvfit.lambda_1se)
            coefs(i,:) = cvglmnetCoef(mdl{i},p.nem.lambda)';

            this_Ypred = cvglmnetPredict(mdl{i},X,p.nem.lambda);
            Ypred(:,i) = this_Ypred;

            this_Y_train = this_Y(these_trainingBins);
            this_Y_test = this_Y(these_testBins);
            this_Ypred_train = this_Ypred(these_trainingBins);
            this_Ypred_test = this_Ypred(these_testBins);

            R2_train(i) = 1 - nansum((this_Y_train - this_Ypred_train).^2) / nansum((this_Y_train - nanmean(this_Y_train)).^2);
            R2_test(i) = 1 - nansum((this_Y_test - this_Ypred_test).^2) / nansum((this_Y_test - nanmean(this_Y_test)).^2);
            R2_all(i) = 1 - nansum((this_Y - this_Ypred).^2) / nansum((this_Y - nanmean(this_Y)).^2);
        catch
        end
    end
end

% discarded alternatives
% nem.mdl{i} = fitlm(nem.predictors.Xdsgn(find(nem.setsByBins~=0),:),nf_binned(idx,find(nem.setsByBins~=0))','linear','Intercept',true) % original
% nem.mdl{i} = lasso(nem.predictors.Xdsgn(find(nem.setsByBins~=0),:),nf_binned(idx,find(nem.setsByBins~=0))','Alpha',0.5,'CV',10)
% nem.mdl{i} = glmnet(nem.predictors.Xdsgn(find(nem.setsByBins~=0),:),nf_binned(idx,find(nem.setsByBins~=0))',p.nem.family,options); % solves for a grid of Lambda (controls overall strength of penalty)

end