function F = templateCorrelationFigure(tempCorr,this_corrType,this_corrMetric,p)

% this_corrType = 'maxBinCorr'; % 'rankCorr', 'maxBinCorr'
% this_corrMetric = 'Spearman'; % 'Pearson', 'Spearman', 'Kendall'

nrows = 2; ncols = 2;
F = default_figure([20,0.5,20,9.9]);

% Correlation with A template in A trials
these_labels = {'iscells','A','X','Aonly','Xonly'};
this_data = [tempCorr.(this_corrType).(this_corrMetric).iscells_Atemp_Atrials,...
    tempCorr.(this_corrType).(this_corrMetric).A_Atemp_Atrials,...
    tempCorr.(this_corrType).(this_corrMetric).X_Atemp_Atrials,...
    tempCorr.(this_corrType).(this_corrMetric).Aonly_Atemp_Atrials,...
    tempCorr.(this_corrType).(this_corrMetric).Xonly_Atemp_Atrials];
subplot(nrows,ncols,1)
hold on
yline(0,':');
h = violinplot(this_data,these_labels);
h(1).ViolinColor = p.col.darkGray; h(2).ViolinColor = p.col.A; h(3).ViolinColor = p.col.X; h(4).ViolinColor = p.col.A; h(5).ViolinColor = p.col.X; 
h(1).BoxColor = 'k'; h(2).BoxColor = 'k'; h(3).BoxColor = 'k'; h(4).BoxColor = 'k'; h(5).BoxColor = 'k';
ylabel('Correlation')
title([this_corrType,' (',this_corrMetric,'): A template in A trials'])

% Correlation with A template in X trials
these_labels = {'iscells','A','X','Aonly','Xonly'};
this_data = [tempCorr.(this_corrType).(this_corrMetric).iscells_Atemp_Xtrials,...
    tempCorr.(this_corrType).(this_corrMetric).A_Atemp_Xtrials,...
    tempCorr.(this_corrType).(this_corrMetric).X_Atemp_Xtrials,...
    tempCorr.(this_corrType).(this_corrMetric).Aonly_Atemp_Xtrials,...
    tempCorr.(this_corrType).(this_corrMetric).Xonly_Atemp_Xtrials];
subplot(nrows,ncols,2)
hold on
yline(0,':');
h = violinplot(this_data,these_labels);
h(1).ViolinColor = p.col.darkGray; h(2).ViolinColor = p.col.A; h(3).ViolinColor = p.col.X; h(4).ViolinColor = p.col.A; h(5).ViolinColor = p.col.X; 
h(1).BoxColor = 'k'; h(2).BoxColor = 'k'; h(3).BoxColor = 'k'; h(4).BoxColor = 'k'; h(5).BoxColor = 'k';
ylabel('Correlation')
title([this_corrType,' (',this_corrMetric,'): A template in X trials'])

% Correlation with X template in A trials
these_labels = {'iscells','A','X','Aonly','Xonly'};
this_data = [tempCorr.(this_corrType).(this_corrMetric).iscells_Xtemp_Atrials,...
    tempCorr.(this_corrType).(this_corrMetric).A_Xtemp_Atrials,...
    tempCorr.(this_corrType).(this_corrMetric).X_Xtemp_Atrials,...
    tempCorr.(this_corrType).(this_corrMetric).Aonly_Xtemp_Atrials,...
    tempCorr.(this_corrType).(this_corrMetric).Xonly_Xtemp_Atrials];
subplot(nrows,ncols,3)
hold on
yline(0,':');
h = violinplot(this_data,these_labels);
h(1).ViolinColor = p.col.darkGray; h(2).ViolinColor = p.col.A; h(3).ViolinColor = p.col.X; h(4).ViolinColor = p.col.A; h(5).ViolinColor = p.col.X; 
h(1).BoxColor = 'k'; h(2).BoxColor = 'k'; h(3).BoxColor = 'k'; h(4).BoxColor = 'k'; h(5).BoxColor = 'k';
ylabel('Correlation')
title([this_corrType,' (',this_corrMetric,'): X template in A trials'])

% Correlation with X template in X trials
these_labels = {'iscells','A','X','Aonly','Xonly'};
this_data = [tempCorr.(this_corrType).(this_corrMetric).iscells_Xtemp_Xtrials,...
    tempCorr.(this_corrType).(this_corrMetric).A_Xtemp_Xtrials,...
    tempCorr.(this_corrType).(this_corrMetric).X_Xtemp_Xtrials,...
    tempCorr.(this_corrType).(this_corrMetric).Aonly_Xtemp_Xtrials,...
    tempCorr.(this_corrType).(this_corrMetric).Xonly_Xtemp_Xtrials];
subplot(nrows,ncols,4)
hold on
yline(0,':');
h = violinplot(this_data,these_labels);
h(1).ViolinColor = p.col.darkGray; h(2).ViolinColor = p.col.A; h(3).ViolinColor = p.col.X; h(4).ViolinColor = p.col.A; h(5).ViolinColor = p.col.X; 
h(1).BoxColor = 'k'; h(2).BoxColor = 'k'; h(3).BoxColor = 'k'; h(4).BoxColor = 'k'; h(5).BoxColor = 'k';
ylabel('Correlation')
title([this_corrType,' (',this_corrMetric,'): X template in X trials'])

for i=1:(nrows*ncols)
    subplot(nrows,ncols,i)
    ylim([-0.3,0.5])
end