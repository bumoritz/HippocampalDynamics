%% --- S -> P ---

%%



%% Figure

this_corrType = 'maxBinCorr'; % 'rankCorr', 'maxBinCorr'
this_corrMetric = 'Spearman'; % 'Pearson', 'Spearman', 'Kendall'

nrows = 2; ncols = 3;
F = default_figure([20,0.5,20,9.9]);

% lower_x = -0.1;
% upper_x = 0.25;
% lower_y = -0.2;
% upper_y = 1.2;


% a)

% x (behaviour, %correct)
x_AB = tempCorr.(this_corrType).(this_corrMetric).A_Atemp_Atrials(trials_all.stimuli_inA.AB)';
x_AY = tempCorr.(this_corrType).(this_corrMetric).A_Atemp_Atrials(trials_all.stimuli_inA.AY)';
x_XY = tempCorr.(this_corrType).(this_corrMetric).A_Atemp_Xtrials(trials_all.stimuli_inX.XY)';
x_XB = tempCorr.(this_corrType).(this_corrMetric).A_Atemp_Xtrials(trials_all.stimuli_inX.XB)';

% y (behaviour, %correct)
y_AB = task.response(find(task.odour1=="A"&task.odour2=="B"));
y_AB = y_AB=="H";
y_AY = task.response(find(task.odour1=="A"&task.odour2=="Y"));
y_AY = y_AY=="CR";
y_XY = task.response(find(task.odour1=="X"&task.odour2=="Y"));
y_XY = y_XY=="H";
y_XB = task.response(find(task.odour1=="X"&task.odour2=="B"));
y_XB = y_XB=="CR";

% A cells: Correlation to A template
subplot(nrows,ncols,1)
hold on
[corr_r_AB,corr_p_AB] = fitLine(x_AB',y_AB',p.col.AB);
[corr_r_XB,corr_p_XB] = fitLine(x_XB',y_XB',p.col.XB);
[corr_r_AY,corr_p_AY] = fitLine(x_AY',y_AY',p.col.AY);
[corr_r_XY,corr_p_XY] = fitLine(x_XY',y_XY',p.col.XY);
% xlim([lower_x,upper_x])
% ylim([lower_y,upper_y])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['AB: n = ',num2str(length(y_AB)),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
    'XB: n = ',num2str(length(y_XB)),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2),newline,...
    'AY: n = ',num2str(length(y_AY)),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2),newline,...
    'XY: n = ',num2str(length(y_XY)),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
legend('AB','XB','AY','XY')
xlabel([this_corrType,' (',this_corrMetric,')'])
ylabel('%Correct')
title('Performance by correlation of A cells with A sequence template')


% b)

% x (behaviour, %correct)
x_AB = tempCorr.(this_corrType).(this_corrMetric).X_Atemp_Atrials(trials_all.stimuli_inA.AB)';
x_AY = tempCorr.(this_corrType).(this_corrMetric).X_Atemp_Atrials(trials_all.stimuli_inA.AY)';
x_XY = tempCorr.(this_corrType).(this_corrMetric).X_Atemp_Xtrials(trials_all.stimuli_inX.XY)';
x_XB = tempCorr.(this_corrType).(this_corrMetric).X_Atemp_Xtrials(trials_all.stimuli_inX.XB)';

% y (behaviour, %correct)
y_AB = task.response(find(task.odour1=="A"&task.odour2=="B"));
y_AB = y_AB=="H";
y_AY = task.response(find(task.odour1=="A"&task.odour2=="Y"));
y_AY = y_AY=="CR";
y_XY = task.response(find(task.odour1=="X"&task.odour2=="Y"));
y_XY = y_XY=="H";
y_XB = task.response(find(task.odour1=="X"&task.odour2=="B"));
y_XB = y_XB=="CR";

% X cells: Correlation to A template
subplot(nrows,ncols,2)
hold on
[corr_r_AB,corr_p_AB] = fitLine(x_AB',y_AB',p.col.AB);
[corr_r_XB,corr_p_XB] = fitLine(x_XB',y_XB',p.col.XB);
[corr_r_AY,corr_p_AY] = fitLine(x_AY',y_AY',p.col.AY);
[corr_r_XY,corr_p_XY] = fitLine(x_XY',y_XY',p.col.XY);
% xlim([lower_x,upper_x])
% ylim([lower_y,upper_y])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['AB: n = ',num2str(length(y_AB)),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
    'XB: n = ',num2str(length(y_XB)),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2),newline,...
    'AY: n = ',num2str(length(y_AY)),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2),newline,...
    'XY: n = ',num2str(length(y_XY)),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
legend('AB','XB','AY','XY')
xlabel([this_corrType,' (',this_corrMetric,')'])
ylabel('%Correct')
title('Performance by correlation of X cells with A sequence template')


% c)

% x (behaviour, %correct)
x_AB = tempCorr.(this_corrType).(this_corrMetric).iscells_Atemp_Atrials(trials_all.stimuli_inA.AB)';
x_AY = tempCorr.(this_corrType).(this_corrMetric).iscells_Atemp_Atrials(trials_all.stimuli_inA.AY)';
x_XY = tempCorr.(this_corrType).(this_corrMetric).iscells_Atemp_Xtrials(trials_all.stimuli_inX.XY)';
x_XB = tempCorr.(this_corrType).(this_corrMetric).iscells_Atemp_Xtrials(trials_all.stimuli_inX.XB)';

% y (behaviour, %correct)
y_AB = task.response(find(task.odour1=="A"&task.odour2=="B"));
y_AB = y_AB=="H";
y_AY = task.response(find(task.odour1=="A"&task.odour2=="Y"));
y_AY = y_AY=="CR";
y_XY = task.response(find(task.odour1=="X"&task.odour2=="Y"));
y_XY = y_XY=="H";
y_XB = task.response(find(task.odour1=="X"&task.odour2=="B"));
y_XB = y_XB=="CR";

% all cells: Correlation to A template
subplot(nrows,ncols,3)
hold on
[corr_r_AB,corr_p_AB] = fitLine(x_AB',y_AB',p.col.AB);
[corr_r_XB,corr_p_XB] = fitLine(x_XB',y_XB',p.col.XB);
[corr_r_AY,corr_p_AY] = fitLine(x_AY',y_AY',p.col.AY);
[corr_r_XY,corr_p_XY] = fitLine(x_XY',y_XY',p.col.XY);
% xlim([lower_x,upper_x])
% ylim([lower_y,upper_y])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['AB: n = ',num2str(length(y_AB)),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
    'XB: n = ',num2str(length(y_XB)),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2),newline,...
    'AY: n = ',num2str(length(y_AY)),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2),newline,...
    'XY: n = ',num2str(length(y_XY)),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
legend('AB','XB','AY','XY')
xlabel([this_corrType,' (',this_corrMetric,')'])
ylabel('%Correct')
title('Performance by correlation of all cells with A sequence template')


% d)

% x (behaviour, %correct)
x_AB = tempCorr.(this_corrType).(this_corrMetric).A_Xtemp_Atrials(trials_all.stimuli_inA.AB)';
x_AY = tempCorr.(this_corrType).(this_corrMetric).A_Xtemp_Atrials(trials_all.stimuli_inA.AY)';
x_XY = tempCorr.(this_corrType).(this_corrMetric).A_Xtemp_Xtrials(trials_all.stimuli_inX.XY)';
x_XB = tempCorr.(this_corrType).(this_corrMetric).A_Xtemp_Xtrials(trials_all.stimuli_inX.XB)';

% y (behaviour, %correct)
y_AB = task.response(find(task.odour1=="A"&task.odour2=="B"));
y_AB = y_AB=="H";
y_AY = task.response(find(task.odour1=="A"&task.odour2=="Y"));
y_AY = y_AY=="CR";
y_XY = task.response(find(task.odour1=="X"&task.odour2=="Y"));
y_XY = y_XY=="H";
y_XB = task.response(find(task.odour1=="X"&task.odour2=="B"));
y_XB = y_XB=="CR";

% A cells: Correlation to X template
subplot(nrows,ncols,4)
hold on
[corr_r_AB,corr_p_AB] = fitLine(x_AB',y_AB',p.col.AB);
[corr_r_XB,corr_p_XB] = fitLine(x_XB',y_XB',p.col.XB);
[corr_r_AY,corr_p_AY] = fitLine(x_AY',y_AY',p.col.AY);
[corr_r_XY,corr_p_XY] = fitLine(x_XY',y_XY',p.col.XY);
% xlim([lower_x,upper_x])
% ylim([lower_y,upper_y])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['AB: n = ',num2str(length(y_AB)),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
    'XB: n = ',num2str(length(y_XB)),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2),newline,...
    'AY: n = ',num2str(length(y_AY)),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2),newline,...
    'XY: n = ',num2str(length(y_XY)),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
legend('AB','XB','AY','XY')
xlabel([this_corrType,' (',this_corrMetric,')'])
ylabel('%Correct')
title('Performance by correlation of A cells with X sequence template')


% e)

% x (behaviour, %correct)
x_AB = tempCorr.(this_corrType).(this_corrMetric).X_Xtemp_Atrials(trials_all.stimuli_inA.AB)';
x_AY = tempCorr.(this_corrType).(this_corrMetric).X_Xtemp_Atrials(trials_all.stimuli_inA.AY)';
x_XY = tempCorr.(this_corrType).(this_corrMetric).X_Xtemp_Xtrials(trials_all.stimuli_inX.XY)';
x_XB = tempCorr.(this_corrType).(this_corrMetric).X_Xtemp_Xtrials(trials_all.stimuli_inX.XB)';

% y (behaviour, %correct)
y_AB = task.response(find(task.odour1=="A"&task.odour2=="B"));
y_AB = y_AB=="H";
y_AY = task.response(find(task.odour1=="A"&task.odour2=="Y"));
y_AY = y_AY=="CR";
y_XY = task.response(find(task.odour1=="X"&task.odour2=="Y"));
y_XY = y_XY=="H";
y_XB = task.response(find(task.odour1=="X"&task.odour2=="B"));
y_XB = y_XB=="CR";

% X cells: Correlation to X template
subplot(nrows,ncols,5)
hold on
[corr_r_AB,corr_p_AB] = fitLine(x_AB',y_AB',p.col.AB);
[corr_r_XB,corr_p_XB] = fitLine(x_XB',y_XB',p.col.XB);
[corr_r_AY,corr_p_AY] = fitLine(x_AY',y_AY',p.col.AY);
[corr_r_XY,corr_p_XY] = fitLine(x_XY',y_XY',p.col.XY);
% xlim([lower_x,upper_x])
% ylim([lower_y,upper_y])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['AB: n = ',num2str(length(y_AB)),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
    'XB: n = ',num2str(length(y_XB)),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2),newline,...
    'AY: n = ',num2str(length(y_AY)),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2),newline,...
    'XY: n = ',num2str(length(y_XY)),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
legend('AB','XB','AY','XY')
xlabel([this_corrType,' (',this_corrMetric,')'])
ylabel('%Correct')
title('Performance by correlation of X cells with X sequence template')



% f)

% x (behaviour, %correct)
x_AB = tempCorr.(this_corrType).(this_corrMetric).iscells_Xtemp_Atrials(trials_all.stimuli_inA.AB)';
x_AY = tempCorr.(this_corrType).(this_corrMetric).iscells_Xtemp_Atrials(trials_all.stimuli_inA.AY)';
x_XY = tempCorr.(this_corrType).(this_corrMetric).iscells_Xtemp_Xtrials(trials_all.stimuli_inX.XY)';
x_XB = tempCorr.(this_corrType).(this_corrMetric).iscells_Xtemp_Xtrials(trials_all.stimuli_inX.XB)';

% y (behaviour, %correct)
y_AB = task.response(find(task.odour1=="A"&task.odour2=="B"));
y_AB = y_AB=="H";
y_AY = task.response(find(task.odour1=="A"&task.odour2=="Y"));
y_AY = y_AY=="CR";
y_XY = task.response(find(task.odour1=="X"&task.odour2=="Y"));
y_XY = y_XY=="H";
y_XB = task.response(find(task.odour1=="X"&task.odour2=="B"));
y_XB = y_XB=="CR";

% all cells: Correlation to X template
subplot(nrows,ncols,6)
hold on
[corr_r_AB,corr_p_AB] = fitLine(x_AB',y_AB',p.col.AB);
[corr_r_XB,corr_p_XB] = fitLine(x_XB',y_XB',p.col.XB);
[corr_r_AY,corr_p_AY] = fitLine(x_AY',y_AY',p.col.AY);
[corr_r_XY,corr_p_XY] = fitLine(x_XY',y_XY',p.col.XY);
% xlim([lower_x,upper_x])
% ylim([lower_y,upper_y])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['AB: n = ',num2str(length(y_AB)),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
    'XB: n = ',num2str(length(y_XB)),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2),newline,...
    'AY: n = ',num2str(length(y_AY)),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2),newline,...
    'XY: n = ',num2str(length(y_XY)),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
legend('AB','XB','AY','XY')
xlabel([this_corrType,' (',this_corrMetric,')'])
ylabel('%Correct')
title('Performance by correlation of all cells with X template')




