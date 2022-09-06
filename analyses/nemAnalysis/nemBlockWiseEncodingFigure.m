function F = nemBlockWiseEncodingFigure(nem_100t,info,p,prop)

numBlocks = length(nem_100t);
numCells = nansum(prop.iscell);

nrows = 3;
ncols = 2;
F = default_figure([20,0.5,20,9.9]);


%% Cells with significant fit

this_data = nan(1,numBlocks);
for j=1:numBlocks
	this_data(j) = nansum(nem_100t{j}.shuffled.significant) / numCells;
end

subplot(nrows,ncols,1)

plot(this_data*100,'Color',p.col.black,'LineWidth',1)
hold on
try
    this_data = nansum(nem_all.shuffled.significant) / numCells;
    plot(0,this_data*100,'o','Color',p.col.black,'LineWidth',1)
catch
end

xlim([-1,numBlocks+2])
ylim([0,100])
ytickformat('percentage')
xlabel('Trial block')
ylabel('Proportion of cells')
title('Cells with significant model fits')


%% Significantly activated cells (all cells)

these_testGroups = [2,3,11,12,13,16,19];

this_data = nan(nem_100t{1}.numTestGroups,numBlocks);
for j=1:numBlocks
    for i=1:nem_100t{1}.numTestGroups
        if ismember(i,these_testGroups)
            this_direction = nanmean(nem_100t{j}.testGroupResidual{i}.coefs(:,cell2mat(nem_100t{j}.predictorGroups.idcs(nem_100t{j}.testGroup.group{i}))),2) > 0;
            this_data(i,j) = nansum(floor((nem_100t{j}.testGroupShuffled{i}.significant + this_direction)/2)) / numCells;
        end
    end
end

subplot(nrows,ncols,3)

for i=1:nem_100t{1}.numTestGroups
    if ismember(i,these_testGroups)
        plot(this_data(i,:)*100,'Color',nem_100t{1}.testGroup.cols{i},'LineWidth',1)
        hold on
    end
end
try
    for i=1:nem_all.numTestGroups
        if ismember(i,these_testGroups)
            this_direction = nanmean(nem_all.testGroupResidual{i}.coefs(:,cell2mat(nem_all.predictorGroups.idcs(nem_all.testGroup.group{i}))),2) > 0;
            this_data = nansum(floor((nem_all.testGroupShuffled{i}.significant + this_direction)/2)) / numCells;
            plot(0,this_data*100,'o','Color',nem_all.testGroup.cols{i},'LineWidth',1)
        end
    end
catch
end

xlim([-1,numBlocks+2])
ylim([0,inf])
ytickformat('percentage')
legend({nem_100t{1}.testGroup.label{these_testGroups}})
xlabel('Trial block')
ylabel('Proportion of cells')
title('Significantly activated cells')


%% Significantly activated cells (cells with sig. fit)

these_testGroups = [2,3,11,12,13,16,19];

this_data = nan(nem_100t{1}.numTestGroups,numBlocks);
for j=1:numBlocks
    for i=1:nem_100t{1}.numTestGroups
        if ismember(i,these_testGroups)
            this_direction = nanmean(nem_100t{j}.testGroupResidual{i}.coefs(:,cell2mat(nem_100t{j}.predictorGroups.idcs(nem_100t{j}.testGroup.group{i}))),2) > 0;
            this_data(i,j) = nansum(floor((nem_100t{j}.testGroupShuffled{i}.significant + this_direction)/2)) / nansum(nem_100t{j}.shuffled.significant);
        end
    end
end

subplot(nrows,ncols,5)

for i=1:nem_100t{1}.numTestGroups
    if ismember(i,these_testGroups)
        plot(this_data(i,:)*100,'Color',nem_100t{1}.testGroup.cols{i},'LineWidth',1)
        hold on
    end
end
try
    for i=1:nem_all.numTestGroups
        if ismember(i,these_testGroups)
            this_direction = nanmean(nem_all.testGroupResidual{i}.coefs(:,cell2mat(nem_all.predictorGroups.idcs(nem_all.testGroup.group{i}))),2) > 0;
            this_data = nansum(floor((nem_all.testGroupShuffled{i}.significant + this_direction)/2)) / nansum(nem_all.shuffled.significant);
            plot(0,this_data*100,'o','Color',nem_all.testGroup.cols{i},'LineWidth',1)
        end
    end
catch
end

xlim([-1,numBlocks+2])
ylim([0,inf])
ytickformat('percentage')
legend({nem_100t{1}.testGroup.label{these_testGroups}})
xlabel('Trial block')
ylabel('Proportion of cells with sig. fit')
title('Significantly activated cells')


%% Cells with significant fit

this_data = nan(1,numBlocks);
for j=1:numBlocks
	this_data(j) = nansum(nem_100t{j}.shuffled.significant) / numCells;
end

subplot(nrows,ncols,2)

plot(this_data*100,'Color',p.col.black,'LineWidth',1)
hold on
try
    this_data = nansum(nem_all.shuffled.significant) / numCells;
    plot(0,this_data*100,'o','Color',p.col.black,'LineWidth',1)
catch
end

xlim([-1,numBlocks+2])
ylim([0,100])
ytickformat('percentage')
xlabel('Trial block')
ylabel('Proportion of cells')
title('Cells with significant model fits')


%% Significantly suppressed cells (all cells)

these_testGroups = [2,3,11,12,13,16,19];

this_data = nan(nem_100t{1}.numTestGroups,numBlocks);
for j=1:numBlocks
    for i=1:nem_100t{1}.numTestGroups
        if ismember(i,these_testGroups)
            this_direction = nanmean(nem_100t{j}.testGroupResidual{i}.coefs(:,cell2mat(nem_100t{j}.predictorGroups.idcs(nem_100t{j}.testGroup.group{i}))),2) < 0;
            this_data(i,j) = nansum(floor((nem_100t{j}.testGroupShuffled{i}.significant + this_direction)/2)) / numCells;
        end
    end
end

subplot(nrows,ncols,4)

for i=1:nem_100t{1}.numTestGroups
    if ismember(i,these_testGroups)
        plot(this_data(i,:)*100,'Color',nem_100t{1}.testGroup.cols{i},'LineWidth',1)
        hold on
    end
end
try
    for i=1:nem_all.numTestGroups
        if ismember(i,these_testGroups)
            this_direction = nanmean(nem_all.testGroupResidual{i}.coefs(:,cell2mat(nem_all.predictorGroups.idcs(nem_all.testGroup.group{i}))),2) < 0;
            this_data = nansum(floor((nem_all.testGroupShuffled{i}.significant + this_direction)/2)) / numCells;
            plot(0,this_data*100,'o','Color',nem_all.testGroup.cols{i},'LineWidth',1)
        end
    end
catch
end

xlim([-1,numBlocks+2])
ylim([0,inf])
ytickformat('percentage')
legend({nem_100t{1}.testGroup.label{these_testGroups}})
xlabel('Trial block')
ylabel('Proportion of cells')
title('Significantly suppressed cells')


%% Significantly suppressed cells (cells with sig. fit)

these_testGroups = [2,3,11,12,13,16,19];

this_data = nan(nem_100t{1}.numTestGroups,numBlocks);
for j=1:numBlocks
    for i=1:nem_100t{1}.numTestGroups
        if ismember(i,these_testGroups)
            this_direction = nanmean(nem_100t{j}.testGroupResidual{i}.coefs(:,cell2mat(nem_100t{j}.predictorGroups.idcs(nem_100t{j}.testGroup.group{i}))),2) < 0;
            this_data(i,j) = nansum(floor((nem_100t{j}.testGroupShuffled{i}.significant + this_direction)/2)) / nansum(nem_100t{j}.shuffled.significant);
        end
    end
end

subplot(nrows,ncols,6)

for i=1:nem_100t{1}.numTestGroups
    if ismember(i,these_testGroups)
        plot(this_data(i,:)*100,'Color',nem_100t{1}.testGroup.cols{i},'LineWidth',1)
        hold on
    end
end
try
    for i=1:nem_all.numTestGroups
        if ismember(i,these_testGroups)
            this_direction = nanmean(nem_all.testGroupResidual{i}.coefs(:,cell2mat(nem_all.predictorGroups.idcs(nem_all.testGroup.group{i}))),2) < 0;
            this_data = nansum(floor((nem_all.testGroupShuffled{i}.significant + this_direction)/2)) / nansum(nem_all.shuffled.significant);
            plot(0,this_data*100,'o','Color',nem_all.testGroup.cols{i},'LineWidth',1)
        end
    end
catch
end

xlim([-1,numBlocks+2])
ylim([0,inf])
ytickformat('percentage')
legend({nem_100t{1}.testGroup.label{these_testGroups}})
xlabel('Trial block')
ylabel('Proportion of cells with sig. fit')
title('Significantly suppressed cells')


%% Return

if info.stimSession
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', Block-wise population-wide encoding']);
else
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-','nostim',', Block-wise population-wide encoding']);
end
drawnow;