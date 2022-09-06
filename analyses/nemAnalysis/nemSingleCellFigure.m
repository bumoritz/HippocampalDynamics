function F = nemSingleCellFigure(nem,nemd,nemdp_all,nem_100t,nemd_100t,trials,traces,idx,info,p,prop)
% nemdp_all=[]; nem_100t=[]; nemd_100t=[]; trials=trials_all; traces=traces_all; idx=146;
% nem=nem_all; nemd=nemd_all;

% add full model to all plots on the right

numBlocks = floor(prop.numTrials/100);
these_conditions = {'AB','AY','XY','XB'};
nrows = 3;
ncols = 4;


F = default_figure([20,0.5,20,9.9]);


%% Average traces and predictor group contributions

temp = reshape(nemd.full(prop.trial_frames_binned,idx)',[],size(prop.trial_frames_binned,1),prop.numTrials);
[traces_pred,~,~] = createTracesStructs(temp,trials);

try
    for j=1:nem.numPredictorGroups
        this_predictorGroup = nem.predictorGroups.names{j};
        temp = reshape(nemdp_all.(this_predictorGroup).Ypred_all(prop.trial_frames_binned,idx)',[],size(prop.trial_frames_binned,1),prop.numTrials);
        [traces_pred_predictorGroup{j},~,~] = createTracesStructs(temp,trials);
    end
catch
end

temp1 = zeros(4,2);
these_subplot_idcs = [1,2,5,6];
for i=1:4
    subplot(nrows,ncols,these_subplot_idcs(i))
    hold on
    
    try
        for j=1:nem.numPredictorGroups
            this_predictorGroup = nem.predictorGroups.names{j};
            plot(nanmean(traces_pred_predictorGroup{j}.(these_conditions{i})(:,:,:),3),'Color',nem.predictorGroups.cols{j},'LineWidth',1,'LineStyle','--');
        end
    catch
    end
    temp=shadedErrorBar(1:size(traces.(these_conditions{i})(idx,:,:),2),nanmean(traces.(these_conditions{i})(idx,:,:),3),nansem(traces.(these_conditions{i})(idx,:,:),3),'lineProps',p.col.black); temp.mainLine.LineWidth = 2;   
    plot(nanmean(traces_pred.(these_conditions{i})(:,:,:),3),'Color',p.col.photostim,'LineWidth',2);
    plt = struct(); traces_task([],p,info,plt);
    temp1(i,:) = get(gca,'ylim');
    ylabel('Activity (z-score)')
    title(these_conditions{i})
end
for i=1:4
    subplot(nrows,ncols,these_subplot_idcs(i))
    ylim([nanmin(temp1(:)),nanmax(temp1(:))])
end


%% Traces

this_window = 1000:2000;

subplot(nrows,ncols,9)

hold on
temp = (1:length(nemd.data(this_window,idx)))*p.general.binSize/info.scope.frameRate;
plot(temp,nemd.data(this_window,idx),'Color',p.col.black,'LineWidth',2)
plot(temp,nemd.full(this_window,idx),'Color',p.col.photostim,'LineWidth',2)

xlabel('Time (s)')
ylabel('Activity (z-score)')


%% Coefficients

subplot(nrows,ncols,10)

h=barh(diag(nem.full.coefs(idx,2:end)),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for i=1:nem.numPredictors
    set(h(i),'FaceColor',nem.predictors.cols(i,:))
    set(h(i),'EdgeColor',nem.predictors.cols(i,:))
end

set(gca,'Ydir','reverse')
xlabel('Coefficient')
ylabel('Predictor')


%% Sum of coefficients
    
this_data = nan(1,nem.numPredictorGroups);
for j=1:nem.numPredictorGroups
    this_predictorGroup = nem.predictorGroups.names{j};
    this_data(j) = nansum(nem.full.coefs(idx,1+nem.predictorGroups.idcs{j}));
end

subplot(nrows,ncols,3:4)

h=bar(diag(this_data),'stacked','BaseValue',0,'EdgeColor','none');
for j=1:nem.numPredictorGroups
    set(h(j),'FaceColor',nem.predictorGroups.cols{j})
end

xticklabels(nem.predictorGroups.labels)
ylabel('Sum of coefficients')
title('Full model')


%% Explained variance of leave-one-group-out models

this_data = nan(1,nem.numTestGroups);
for i=1:nem.numTestGroups
    this_data(i) = nem.testGroupOut{i}.R2_test(idx);
end

subplot(nrows,ncols,7:8)

h=bar(diag(this_data*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for i=1:nem.numTestGroups
    set(h(i),'FaceColor',nem.testGroup.cols{i})
end
hold on
yline(nem.full.R2_test(idx)*100,'LineStyle',':');

ylim([0,(nem.full.R2_test(idx)+0.1*nem.full.R2_test(idx))*100])
ytickformat('percentage')
xticks(1:length(this_data))
xticklabels(nem.testGroup.label)
ylabel('Explained variance')
set(gca,'box','off')
title('Leave-one-group-out models')


%% Modeling residuals from leave-one-group-out models

this_data = nan(1,nem.numTestGroups);
for i=1:nem.numTestGroups
    this_data(i) = nem.testGroupResidual{i}.R2_test(idx);
end

subplot(nrows,ncols,11:12)

h=bar(diag(this_data*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for i=1:nem.numTestGroups
    set(h(i),'FaceColor',nem.testGroup.cols{i})
end
hold on
yline(nem.full.R2_test(idx)*100,'LineStyle',':');

ylim([0,(nem.full.R2_test(idx)+0.1*nem.full.R2_test(idx))*100])
ytickformat('percentage')
xticks(1:length(this_data))
xticklabels(nem.testGroup.label)
ylabel('Unique explained variance')
set(gca,'box','off')
title('Modeling residuals from leave-one-group-out models')


%% Return

if info.stimSession
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', idx=',num2str(idx),', R2=',num2str(nem.full.R2_test(idx),2)]);
else
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-','nostim',', idx=',num2str(idx),', R2=',num2str(nem.full.R2_test(idx),2)]);
end
drawnow;
end


%% --- DISCARDED ---


% if false
% 
%     %% Std of partial prediction
% 
%     this_data = nan(1,nem.numPredictorGroups);
%     for j=1:nem.numPredictorGroups
%         this_predictorGroup = nem.predictorGroups.names{j};
%         this_data(j) = nanstd(nemd.predictorGroup.(this_predictorGroup).Ypred_all(:,idx));
%     end
% 
%     subplot(nrows,ncols,7)
% 
%     h=bar(diag(this_data),'stacked','BaseValue',0,'EdgeColor','none');
%     for j=1:nem.numPredictorGroups
%         set(h(j),'FaceColor',nem.predictorGroups.cols{j})
%     end
% 
%     xticklabels(nem.predictorGroups.labels)
%     ylabel('Std of partial prediction (both sets)')
% 
% 
%     %% Partial variance explained
% 
%     this_data = nan(1,nem.numPredictorGroups);
%     for j=1:nem.numPredictorGroups
%         this_predictorGroup = nem.predictorGroups.names{j};
%         this_data(j) = nem.R2.predictorGroup.(this_predictorGroup).test(idx) / nem.R2.test(idx);
%     end
% 
%     subplot(nrows,ncols,11)
% 
%     h=bar(diag(this_data*100),'stacked','BaseValue',0,'EdgeColor','none');
%     for j=1:nem.numPredictorGroups
%         set(h(j),'FaceColor',nem.predictorGroups.cols{j})
%     end
% 
%     xticklabels(nem.predictorGroups.labels)
%     ylim([-30,100])
%     ytickformat('percentage')
%     ylabel('Rel. partial variance explained (test set)')
% 
% 
%     %% Metrics across trials
% 
%     if ~isempty(nem_100t) && ~isempty(nemd_100t)
% 
%         %% Sum of coefficients across trials
% 
%         this_data = nan(numBlocks,nem.numPredictorGroups);
%         for j=1:nem.numPredictorGroups
%             this_predictorGroup = nem_100t{1}.predictorGroups.names{j};
%             for k=1:numBlocks
%                 this_data(k,j) = nansum(nem_100t{k}.coefs(idx,1+nem_100t{k}.predictorGroups.idcs{j}));
%             end
%         end
% 
%         subplot(nrows,ncols,4)
% 
%         for j=1:nem.numPredictorGroups
%             plot(this_data(:,j),'Color',nem_100t{1}.predictorGroups.cols{j},'LineWidth',1)
%             hold on
%         end
% 
%         xlabel('Trial block')
%         ylabel('Sum of coefficients')
% 
% 
%         %% Std of partial prediction across trials
% 
%         this_data = nan(numBlocks,nem.numPredictorGroups);
%         for j=1:nem.numPredictorGroups
%             this_predictorGroup = nem_100t{1}.predictorGroups.names{j};
%             for k=1:numBlocks
%                 this_data(k,j) = nanstd(nemd_100t{k}.predictorGroup.(this_predictorGroup).Ypred_all(:,idx));
%             end
%         end
% 
%         subplot(nrows,ncols,8)
% 
%         for j=1:nem.numPredictorGroups
%             plot(this_data(:,j),'Color',nem_100t{1}.predictorGroups.cols{j},'LineWidth',1)
%             hold on
%         end
% 
%         xlabel('Trial block')
%         ylabel('Std of partial prediction (both sets)')
% 
% 
%         %% Partial variance explained across trials
% 
%         this_data = nan(numBlocks,nem.numPredictorGroups);
%         for j=1:nem.numPredictorGroups
%             this_predictorGroup = nem_100t{1}.predictorGroups.names{j};
%             for k=1:numBlocks
%                 this_data(k,j) = nem_100t{k}.R2.predictorGroup.(this_predictorGroup).test(idx);
%             end
%         end
% 
%         subplot(nrows,ncols,12)
% 
%         for j=1:nem.numPredictorGroups
%             plot(this_data(:,j),'Color',nem_100t{1}.predictorGroups.cols{j},'LineWidth',1)
%             hold on
%         end
% 
%         xlabel('Trial block')
%         ylabel('Partial variance explained (test set)')
%     end
% 
% end