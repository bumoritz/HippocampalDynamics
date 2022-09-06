function F = nemPopulationWideEncodingFigure(nem,info,prop)

%% Preparations

numCells = nansum(prop.iscell);

nrows = 2;
ncols = 2;

F = default_figure([20,0.5,20,9.9]);


%% Significant unique explained variance

this_data = nan(1,nem.numTestGroups);
for i=1:nem.numTestGroups
    this_data(i) = nansum(nem.testGroupShuffled{i}.significant) / numCells;
end

subplot(nrows,ncols,1)
h=bar(diag(this_data*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for i=1:nem.numTestGroups
    set(h(i),'FaceColor',nem.testGroup.cols{i})
end

ylim([0,100])
ytickformat('percentage')
xticks(1:length(this_data))
xticklabels(nem.testGroup.label)
ylabel('Proportion of neurons')
set(gca,'box','off')
title('Significant unique explained variance')


%% Significant unique explained variance >=1%

this_data = nan(1,nem.numTestGroups);
for i=1:nem.numTestGroups
    this_relevant = nem.testGroupResidual{i}.R2_test >= 0.01;
    this_data(i) = nansum(floor((nem.testGroupShuffled{i}.significant + this_relevant)/2)) / numCells;
end

subplot(nrows,ncols,2)
h=bar(diag(this_data*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for i=1:nem.numTestGroups
    set(h(i),'FaceColor',nem.testGroup.cols{i})
end

ylim([0,100])
ytickformat('percentage')
xticks(1:length(this_data))
xticklabels(nem.testGroup.label)
ylabel('Proportion of neurons')
set(gca,'box','off')
title('Significant unique explained variance >=1%')


%% Significant unique explained variance (split by sign of influence)

this_data1 = nan(1,nem.numTestGroups);
for i=1:nem.numTestGroups
    this_direction = nanmean(nem.testGroupResidual{i}.coefs(:,cell2mat(nem.predictorGroups.idcs(nem.testGroup.group{i}))),2) > 0;
    this_data1(i) = nansum(floor((nem.testGroupShuffled{i}.significant + this_direction)/2)) / numCells;
end
this_data2 = nan(1,nem.numTestGroups);
for i=1:nem.numTestGroups
    this_direction = nanmean(nem.testGroupResidual{i}.coefs(:,cell2mat(nem.predictorGroups.idcs(nem.testGroup.group{i}))),2) < 0;
    this_data2(i) = nansum(floor((nem.testGroupShuffled{i}.significant + this_direction)/2)) / numCells;
end

subplot(nrows,ncols,3)
h=bar(diag(this_data1*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for i=1:nem.numTestGroups
    set(h(i),'FaceColor',nem.testGroup.cols{i})
end
hold on
h=bar(diag(-this_data2*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for i=1:nem.numTestGroups
    set(h(i),'FaceColor',nem.testGroup.cols{i})
end
yline(0);

ylim([-100,100])
ytickformat('percentage')
xticks(1:length(this_data1))
xticklabels(nem.testGroup.label)
yticks(-100:50:100)
yticklabels({'100%','50%','0%','50%','100%'})
ylabel('Proportion of neurons')
set(gca,'box','off')
title('Significant unique explained variance (split by sign of influence)')


%% Significant unique explained variance >=1% (split by sign of influence)

this_data1 = nan(1,nem.numTestGroups);
for i=1:nem.numTestGroups
    this_relevant = nem.testGroupResidual{i}.R2_test >= 0.01;
    this_direction = nanmean(nem.testGroupResidual{i}.coefs(:,cell2mat(nem.predictorGroups.idcs(nem.testGroup.group{i}))),2) > 0;
    this_data1(i) = nansum(floor((nem.testGroupShuffled{i}.significant + this_relevant + this_direction)/3)) / numCells;
end
this_data2 = nan(1,nem.numTestGroups);
for i=1:nem.numTestGroups
    this_relevant = nem.testGroupResidual{i}.R2_test >= 0.01;
    this_direction = nanmean(nem.testGroupResidual{i}.coefs(:,cell2mat(nem.predictorGroups.idcs(nem.testGroup.group{i}))),2) < 0;
    this_data2(i) = nansum(floor((nem.testGroupShuffled{i}.significant + this_relevant + this_direction)/3)) / numCells;
end

subplot(nrows,ncols,4)
h=bar(diag(this_data1*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for i=1:nem.numTestGroups
    set(h(i),'FaceColor',nem.testGroup.cols{i})
end
hold on
h=bar(diag(-this_data2*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for i=1:nem.numTestGroups
    set(h(i),'FaceColor',nem.testGroup.cols{i})
end
yline(0);

ylim([-100,100])
ytickformat('percentage')
xticks(1:length(this_data1))
xticklabels(nem.testGroup.label)
yticks(-100:50:100)
yticklabels({'100%','50%','0%','50%','100%'})
ylabel('Proportion of neurons')
set(gca,'box','off')
title('Significant unique explained variance >=1% (split by sign of influence)')


%% Return

if info.stimSession
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', Population-wide encoding']);
else
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-','nostim',', Population-wide encoding']);
end
drawnow;
end




