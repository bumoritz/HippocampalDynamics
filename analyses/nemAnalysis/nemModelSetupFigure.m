function F = nemModelSetupFigure(nem,events_binned,these_trials,info,p)
%these_trials = 50:55; %nem = nem_all;

%% Preparations

temp = find(events_binned.sync);
these_bins = temp(min(these_trials))-length(p.general.bins_pre)+1:temp(max(these_trials)+1)-length(p.general.bins_pre);

these_trialStarts = find(events_binned.sync(these_bins));
these_trialMiddles = these_trialStarts+[diff([these_trialStarts;length(events_binned.sync(these_bins))])]/2;
these_trialLabels = {};
for i=1:length(these_trials)
    these_trialLabels{i} = ['Trial ',num2str(these_trials(i))];
end

numBins = length(these_bins);


%% Figure layout

F = default_figure([20,0.5,20,9.9]);
nrows = 4;
ncols = 2;


%% Stimulus overview

this_row = 1;
this_col = 1;
subplot(nrows,ncols,(this_col-1)*ncols+this_row)
hold on

% Odour A
temp1 = find(events_binned.odourA(these_bins));
temp2 = find(events_binned.odourA_off(these_bins));
this_data = zeros(numBins,1);
for i=1:length(temp1)
    this_data(temp1(i):temp2(i))=1;
end
plot(this_data+3.5,'Color',p.col.A)

% Odour X
temp1 = find(events_binned.odourX(these_bins));
temp2 = find(events_binned.odourX_off(these_bins));
this_data = zeros(numBins,1);
for i=1:length(temp1)
    this_data(temp1(i):temp2(i))=1;
end
plot(this_data+3.5,'Color',p.col.X)

% Odour B
temp1 = find(events_binned.odourB(these_bins));
temp2 = find(events_binned.odourB_off(these_bins));
this_data = zeros(numBins,1);
for i=1:length(temp1)
    this_data(temp1(i):temp2(i))=1;
end
plot(this_data+3.5,'Color',p.col.B)

% Odour Y
temp1 = find(events_binned.odourY(these_bins));
temp2 = find(events_binned.odourY_off(these_bins));
this_data = zeros(numBins,1);
for i=1:length(temp1)
    this_data(temp1(i):temp2(i))=1;
end
plot(this_data+3.5,'Color',p.col.Y)

% Response window
temp1 = find(events_binned.responsewindow(these_bins));
temp2 = find(events_binned.responsewindow_off(these_bins));
this_data = zeros(numBins,1);
for i=1:length(temp1)
    this_data(temp1(i):temp2(i))=1;
end
plot(this_data+2,'Color',p.col.gray)

% Reward
temp1 = find(events_binned.reward(these_bins));
temp2 = find(events_binned.reward_off(these_bins));
this_data = zeros(numBins,1);
for i=1:length(temp1)
    this_data(temp1(i):temp2(i))=1;
end
plot(this_data+2,'Color',p.col.reward)

% Lick
this_data = zeros(numBins,1);
this_data(find(events_binned.lick(these_bins)))=1;
for i=1:numBins
    if this_data(i)
        plot([i,i],[0,1]+0.5,'Color',p.col.darkGray)
    else
        plot(i,0.5,'Color',p.col.darkGray)
    end
end

for i=1:length(these_trialStarts)
    xline(these_trialStarts(i),'LineStyle',':');
end
xticks(these_trialMiddles);
xticklabels(these_trialLabels);

xlim([1,numBins])
yticks([1,2.5,4]);
yticklabels({'Licks','Outcome','Odours'});
title('Task events')


%% Predictors

this_row = 1;
this_col = 2:nrows;
subplot(nrows,ncols,(this_col-1)*ncols+this_row)

imagesc(nem.predictors.Xdsgn(these_bins,:)')
set(gca,'CLim',[-3,3]);
colormap(gca,redblue);

n=0;
these_ticks = nan(1,nem.numPredictorGroups);
for i=1:nem.numPredictorGroups
    temp = length(nem.predictorGroups.idcs{i});
    these_ticks(i) = n+temp/2+0.5;
    n=n+temp;
    yline(n+0.5,'LineStyle',':');
end
yticks(these_ticks);
yticklabels(nem.predictorGroups.labels);

n=0;
these_ticks = nan(1,nem.predictorGroups.numSupGroups);
for i=1:nem.predictorGroups.numSupGroups
    for j=1:nansum(nem.predictorGroups.supGroup==i)
        temp1 = find(nem.predictorGroups.supGroup==i);
        
        temp = length(nem.predictorGroups.idcs{temp1(j)});
        n=n+temp;
    end
    these_ticks(i) = n+temp/2+0.5;
    yline(n+0.5,'LineStyle','-');
end

hold on
for i=1:length(these_trialStarts)
    xline(these_trialStarts(i),'LineStyle',':');
end
xticks(these_trialMiddles);
xticklabels(these_trialLabels);
title(['Predictors (',num2str(nem.numPredictors),')'])


%% Predictor correlations

this_row = 2;
this_col = 1:3;
subplot(nrows,ncols,(this_col-1)*ncols+this_row)

imagesc(nem.predictors.predictorCorr)
set(gca,'CLim',[-1,1]);
colormap(gca,redblue);
daspect([1,1,1])

n=0;
these_ticks = nan(1,nem.numPredictorGroups);
for i=1:nem.numPredictorGroups
    temp = length(nem.predictorGroups.idcs{i});
    these_ticks(i) = n+temp/2+0.5;
    n=n+temp;
    xline(n+0.5,'LineStyle',':');
    yline(n+0.5,'LineStyle',':');
end
xticks(these_ticks);
xticklabels(nem.predictorGroups.labels);
yticks(these_ticks);
yticklabels(nem.predictorGroups.labels);

n=0;
these_ticks = nan(1,nem.predictorGroups.numSupGroups);
for i=1:nem.predictorGroups.numSupGroups
    for j=1:nansum(nem.predictorGroups.supGroup==i)
        temp1 = find(nem.predictorGroups.supGroup==i);
        
        temp = length(nem.predictorGroups.idcs{temp1(j)});
        n=n+temp;
    end
    these_ticks(i) = n+temp/2+0.5;
    xline(n+0.5,'LineStyle','-');
    yline(n+0.5,'LineStyle','-');
end

title('Predictor correlations')


%% Cross-validation sets

this_row = 2;
this_col = 4;
subplot(nrows,ncols,(this_col-1)*ncols+this_row)

hold on
temp = (1:length(nem.setsByBins))*p.general.binSize/info.scope.frameRate/60;
plot(temp(nem.setsByBins==0),nem.setsByBins(nem.setsByBins==0),'g.')
plot(temp(nem.setsByBins>0),nem.setsByBins(nem.setsByBins>0),'b.')
ylim([-1,p.nem.cvFolds+1])

xlabel('Time (min)')
ylabel('Set')
title('Cross-validation sets')


%% Return

if info.stimSession
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', Model setup']);
else
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-','nostim',', Model setup']);
end
drawnow;
end