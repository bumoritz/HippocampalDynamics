function [nem,nemd] = nemCore(this_nf_binned,these_events_binned,this_task,info,p,prop,these_cells)
% this_nf_binned = nf_binned; these_events_binned = events_binned; this_task = task; these_cells = 1:length(prop.iscell); %165; %1:length(prop.iscell); %165 

%% Preparations

[nem,nemd] = nemPrepareModeling(this_nf_binned,these_events_binned,this_task,info,p);


%% Fit full models
% 1.68 min for 586 cells, 800 trials

t_start = tic;
[nemd.full.mdl,nemd.full.coefs,nemd.full.Ypred,nem.full.R2_train,nem.full.R2_test,nem.full.R2_all] = ...
    nemFitModel(nem,nemd,p,prop.iscell,these_cells,nem.predictors.Xdsgn);
disp(['--- Fit full models in ',num2str(toc(t_start)/60,3),' min.'])


%% Fit shuffled models
% 32.8 min (~20 * 1.68 min) for 586 cells, 800 trials, 20 shuffles

if p.nem.do_shuffled
    t_start = tic;
    nem.shuffled.R2_train = nan(length(prop.iscell),p.nem.numShuffles);
    nem.shuffled.R2_test = nan(length(prop.iscell),p.nem.numShuffles);
    nem.shuffled.R2_all = nan(length(prop.iscell),p.nem.numShuffles);
    for j=1:p.nem.numShuffles
        
        % extract X and do shuffling
        this_X = nem.predictors.Xdsgn;
        this_X(:,1:nem.numPredictors) = circshift(this_X(:,1:nem.numPredictors),randi([-round(size(nem.predictors.Xdsgn,1)/2),round(size(nem.predictors.Xdsgn,1)/2)],1,1),1);

        % fit model
        [~,~,~,nem.shuffled.R2_train(:,j),nem.shuffled.R2_test(:,j),nem.shuffled.R2_all(:,j)] = ...
            nemFitModel(nem,nemd,p,prop.iscell,these_cells,this_X);
    end
    
    % test model significance
    nem.shuffled.significant = double(nem.full.R2_test > nanmean(nem.shuffled.R2_test,2) + p.nem.modelSignificanceStd*nanstd(nem.shuffled.R2_test,[],2));
    nem.shuffled.significant(isnan(nem.full.R2_test))=NaN;
    
    % store compressed output
    nem.cmpr.sigM = nem.shuffled.significant;

    disp(['--- Fit shuffled models in ',num2str(toc(t_start)/60,3),' min.'])
end


%% Fit test-group-shuffled models
% expected 20 testGroups * 30 min = 10 h

if p.nem.do_testGroupShuffled
    for k=1:nem.numTestGroups
        t_start = tic;
        nem.testGroupShuffled{k}.R2_train = nan(length(prop.iscell),p.nem.numShuffles);
        nem.testGroupShuffled{k}.R2_test = nan(length(prop.iscell),p.nem.numShuffles);
        nem.testGroupShuffled{k}.R2_all = nan(length(prop.iscell),p.nem.numShuffles);
        for j=1:p.nem.numShuffles
            
            % extract X and do shuffling
            this_X = nem.predictors.Xdsgn;
            if ~isempty(nem.testGroup.group)
%                 these_idcs = [];
%                 for i=1:length(nem.testGroup.group{k})
%                     these_idcs = [these_idcs,nem.predictorGroups.idcs{nem.testGroup.group{k}(i)}];
%                 end
                these_idcs = cell2mat(nem.predictorGroups.idcs(nem.testGroup.group{k}));
            elseif ~isempty(nem.testGroup.idcs)
                these_idcs = nem.testGroup.idcs{k};
            end
            this_X(:,these_idcs) = circshift(this_X(:,these_idcs),randi([-round(size(nem.predictors.Xdsgn,1)/2),round(size(nem.predictors.Xdsgn,1)/2)],1,1),1);
            
            % fit model
            [~,~,~,nem.testGroupShuffled{k}.R2_train(:,j),nem.testGroupShuffled{k}.R2_test(:,j),nem.testGroupShuffled{k}.R2_all(:,j)] = ...
                nemFitModel(nem,nemd,p,prop.iscell,these_cells,this_X);
        end
        
        % test predictor group significance
        nem.testGroupShuffled{k}.significant = double(nem.full.R2_test > nanmean(nem.testGroupShuffled{k}.R2_test,2) + p.nem.modelSignificanceStd*nanstd(nem.testGroupShuffled{k}.R2_test,[],2));
        nem.testGroupShuffled{k}.significant(isnan(nem.full.R2_test))=NaN;
        
        % store compressed output
        nem.cmpr.sigTG{k} = nem.testGroupShuffled{k}.significant;
        nem.cmpr.sigM_sigTG{k} = floor((nem.cmpr.sigM + nem.testGroupShuffled{k}.significant)/2);
   
        disp(['--- Fit test-group-shuffled models for ',nem.testGroup.label{k},' in ',num2str(toc(t_start)/60,3),' min.'])
    end
end


%% Fit test-group-out and test-group-residual models

if p.nem.do_testGroupResidual
    for k=1:nem.numTestGroups
        t_start = tic;
        
        % identify predictors from this test shuffle group
        if ~isempty(nem.testGroup.group)
%             these_idcs = [];
%             for i=1:length(nem.testGroup.group{k})
%                 these_idcs = [these_idcs,nem.predictorGroups.idcs{nem.testGroup.group{k}(i)}];
%             end
        	these_idcs = cell2mat(nem.predictorGroups.idcs(nem.testGroup.group{k}));
        elseif ~isempty(nem.testGroup.idcs)
            these_idcs = nem.testGroup.idcs{k};
        end

        % fit testGroupOut models
        this_exclude = these_idcs';
        [nemd.testGroupOut{k}.mdl,nem.testGroupOut{k}.coefs,nemd.testGroupOut{k}.Ypred,nem.testGroupOut{k}.R2_train,nem.testGroupOut{k}.R2_test,nem.testGroupOut{k}.R2_all] = ...
            nemFitModel(nem,nemd,p,prop.iscell,these_cells,nem.predictors.Xdsgn,this_exclude);
        
        % calculate residual and use it as Y
        this_Y = nemd.data - nemd.testGroupOut{k}.Ypred;
            
        % fit residual models
        this_exclude = setdiff(1:nem.numPredictors,these_idcs)';
        [nemd.testGroupResidual{k}.mdl,nem.testGroupResidual{k}.coefs,nemd.testGroupResidual{k}.Ypred,nem.testGroupResidual{k}.R2_train,nem.testGroupResidual{k}.R2_test,nem.testGroupResidual{k}.R2_all] = ...
            nemFitModel(nem,nemd,p,prop.iscell,these_cells,nem.predictors.Xdsgn,this_exclude,this_Y);
        
        % store compressed output
        nem.cmpr.sigTG_predictors{k} = cell2mat(nem.predictorGroups.idcs(nem.testGroup.group{k}));
        nem.cmpr.sigTG_influence{k} = nem.testGroupResidual{k}.coefs(:,nem.cmpr.sigTG_predictors{k}+1);
        nem.cmpr.sigTG_avgInfluence{k} = nanmean(nem.cmpr.sigTG_influence{k},2);
        nem.cmpr.sigTG_pos{k} = floor((nem.cmpr.sigTG{k} + (nem.cmpr.sigTG_avgInfluence{k}>0))/2);
        nem.cmpr.sigTG_neg{k} = floor((nem.cmpr.sigTG{k} + (nem.cmpr.sigTG_avgInfluence{k}<0))/2);
        nem.cmpr.sigM_sigTG_pos{k} = floor((nem.cmpr.sigM_sigTG{k} + (nem.cmpr.sigTG_avgInfluence{k}>0))/2);
        nem.cmpr.sigM_sigTG_neg{k} = floor((nem.cmpr.sigM_sigTG{k} + (nem.cmpr.sigTG_avgInfluence{k}<0))/2);
    
        disp(['--- Fit test-group-residual models for ',nem.testGroup.label{k},' in ',num2str(toc(t_start)/60,3),' min.'])
    end
end


%% Return

nem.info = info;
nem.p = p;
nem.prop = prop;
try
    nem.cmpr = orderfields(nem.cmpr);
catch
end
nem = orderfields(nem);
nemd = orderfields(nemd);

end


%% --- NOT UP TO DATE ---


%% Fit reduced-rank model

% t_start = tic;
% [nem.mdl,nem.coefs,nemd.Ypred_train,nemd.Ypred_test,nemd.Ypred_all,nem.R2.train,nem.R2.test,nem.R2.all] = ...
%     nemFitModel(nem,nemd,p,prop.iscell,these_cells,nemd.Xrr_train,nemd.Xrr_test,nemd.Xrr_all);
% disp(['--- Fit basic models in ',num2str(toc(t_start)/60,3),' min.'])


%% Fit shuffle-one-group model

% if p.nem.do_shuffle1group
%     for k=1:nem.numPredictorGroups
%         t_start = tic;
%         this_predictorGroup = nem.predictorGroups.names{k};
%         nem.R2.shuffle1group.(this_predictorGroup).train = nan(length(prop.iscell),p.nem.numShuffles);
%         nem.R2.shuffle1group.(this_predictorGroup).test = nan(length(prop.iscell),p.nem.numShuffles);
%         nem.R2.shuffle1group.(this_predictorGroup).all = nan(length(prop.iscell),p.nem.numShuffles);
%         for j=1:p.nem.numShuffles
%             
%             % extract X and do shuffling
%             temp = nem.predictors.Xdsgn;
%             temp(:,nem.predictorGroups.idcs{k}) = circshift(temp(:,nem.predictorGroups.idcs{k}),randi([-round(size(nem.predictors.Xdsgn,1)/2),round(size(nem.predictors.Xdsgn,1)/2)],1,1),1);
%             this_X_train = temp(find(nem.setsByBins~=0),:);
%             this_X_test = temp(find(nem.setsByBins==0),:);
%             this_X_all = temp;
%             
%             % fit model
%             [~,~,~,~,~,nem.R2.shuffle1group.(this_predictorGroup).train(:,j),nem.R2.shuffle1group.(this_predictorGroup).test(:,j),nem.R2.shuffle1group.(this_predictorGroup).all(:,j)] = ...
%                 nemFitModel(nem,nemd,p,prop.iscell,these_cells,this_X_train,this_X_test,this_X_all);
%         end
%         
%         % test predictor group significance
%         nem.significantPredictorGroups.(this_predictorGroup) = double(nem.R2.test > nanmean(nem.R2.shuffle1group.(this_predictorGroup).test,2) + p.nem.modelSignificanceStd*nanstd(nem.R2.shuffle1group.(this_predictorGroup).test,[],2));
%         nem.significantPredictorGroups.(this_predictorGroup)(isnan(nem.R2.test))=NaN;
%         
%         disp(['--- Fit shuffle-one-group models for ',this_predictorGroup,' in ',num2str(toc(t_start)/60,3),' min.'])
%     end
% end


%% Fit leave-one-group-out model

% if p.nem.do_leave1out
% 	t_start = tic;
%     for k=1:nem.numPredictorGroups
%         this_predictorGroup = nem.predictorGroups.names{k};
%         
%         % exclude predictors from this group
%         this_exclude = nem.predictorGroups.idcs{k}';
% 
%         % fit model
%         [nem.leave1out.(this_predictorGroup).mdl,nem.leave1out.(this_predictorGroup).coefs,nemd.leave1out.(this_predictorGroup).Ypred_train,nemd.leave1out.(this_predictorGroup).Ypred_test,nemd.leave1out.(this_predictorGroup).Ypred_all,nem.R2.leave1out.(this_predictorGroup).train,nem.R2.leave1out.(this_predictorGroup).test,nem.R2.leave1out.(this_predictorGroup).all] = ...
%             nemFitModel(nem,nemd,p,prop.iscell,these_cells,nemd.X_train,nemd.X_test,nemd.X_all,this_exclude);
%     end
%     disp(['--- Fit leave-one-group-out models in ',num2str(toc(t_start)/60,3),' min.'])
% end


%% Fit leave-one-group-out-residual model

% if p.nem.do_leave1outr
% 	t_start = tic;
%     for k=1:nem.numPredictorGroups
%         this_predictorGroup = nem.predictorGroups.names{k};
% 
%         % exclude predictors from this group
%         this_exclude = setdiff(1:nem.numPredictors,nem.predictorGroups.idcs{k})';
%         
%         % calculate residual and use it as Y
%         this_Y_train = nemd.Y_train - nemd.leave1out.(this_predictorGroup).Ypred_train;
%         this_Y_test = nemd.Y_test - nemd.leave1out.(this_predictorGroup).Ypred_test;
%         this_Y_all = nemd.Y_all - nemd.leave1out.(this_predictorGroup).Ypred_all;
%             
%         % fit model
%         [nem.leave1outr.(this_predictorGroup).mdl,nem.leave1outr.(this_predictorGroup).coefs,nemd.leave1outr.(this_predictorGroup).Ypred_train,nemd.leave1outr.(this_predictorGroup).Ypred_test,nemd.leave1outr.(this_predictorGroup).Ypred_all,nem.R2.leave1outr.(this_predictorGroup).train,nem.R2.leave1outr.(this_predictorGroup).test,nem.R2.leave1outr.(this_predictorGroup).all] = ...
%             nemFitModel(nem,nemd,p,prop.iscell,these_cells,nemd.X_train,nemd.X_test,nemd.X_all,this_exclude,this_Y_train,this_Y_test,this_Y_all);
%     end
%     disp(['--- Fit leave-one-group-out-residual models in ',num2str(toc(t_start)/60,3),' min.'])
% end
