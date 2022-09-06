%% Preparations

if true
    % normal traces
    this_nft = nft_binned;
    this_suptitle = 'traces';
elseif false
    % B residuals
    this_data = nemd_all.data' - nemd_all.testGroupOut{11}'; %nemd_all.data';
	this_nft = reshape(this_data(:,prop.trial_frames_binned),[],size(prop.trial_frames_binned,1),prop.numTrials);
    this_suptitle = 'unique B residuals';
elseif false
    % Y residuals
    this_data = nemd_all.data' - nemd_all.testGroupOut{12}'; %nemd_all.data';
	this_nft = reshape(this_data(:,prop.trial_frames_binned),[],size(prop.trial_frames_binned,1),prop.numTrials);
    this_suptitle = 'unique Y residuals';
end

% idcs
idcs_B = find(nem_all_cmpr.sigM_sigTG_pos{11}==1); % B
idcs_Y = find(nem_all_cmpr.sigM_sigTG_pos{12}==1); % Y
these_idcs = setdiff(idcs_B,idcs_Y);

% X (activity, odour B cells)
B_2nd_AB = squeeze(nanmean(this_nft( these_idcs, p.general.bins_odour2Window , find(task.odour1=="A"&task.odour2=="B") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="A"&task.odour2=="B") ),2));
B_2nd_AY = squeeze(nanmean(this_nft( these_idcs, p.general.bins_odour2Window , find(task.odour1=="A"&task.odour2=="Y") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="A"&task.odour2=="Y") ),2));
B_2nd_XY = squeeze(nanmean(this_nft( these_idcs, p.general.bins_odour2Window , find(task.odour1=="X"&task.odour2=="Y") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="X"&task.odour2=="Y") ),2));
B_2nd_XB = squeeze(nanmean(this_nft( these_idcs, p.general.bins_odour2Window , find(task.odour1=="X"&task.odour2=="B") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="X"&task.odour2=="B") ),2));

% idcs
idcs_B = find(nem_all_cmpr.sigM_sigTG_pos{11}==1); % B
idcs_Y = find(nem_all_cmpr.sigM_sigTG_pos{12}==1); % Y
these_idcs = setdiff(idcs_Y,idcs_B);

% X (activity, odour Y cells)
Y_2nd_AB = squeeze(nanmean(this_nft( these_idcs, p.general.bins_odour2Window , find(task.odour1=="A"&task.odour2=="B") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="A"&task.odour2=="B") ),2));
Y_2nd_AY = squeeze(nanmean(this_nft( these_idcs, p.general.bins_odour2Window , find(task.odour1=="A"&task.odour2=="Y") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="A"&task.odour2=="Y") ),2));
Y_2nd_XY = squeeze(nanmean(this_nft( these_idcs, p.general.bins_odour2Window , find(task.odour1=="X"&task.odour2=="Y") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="X"&task.odour2=="Y") ),2));
Y_2nd_XB = squeeze(nanmean(this_nft( these_idcs, p.general.bins_odour2Window , find(task.odour1=="X"&task.odour2=="B") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="X"&task.odour2=="B") ),2));

% idcs
idcs_A = find(nem_all_cmpr.sigM_sigTG_pos{2}==1); % A
idcs_X = find(nem_all_cmpr.sigM_sigTG_pos{3}==1); % X
these_idcs = setdiff(idcs_A,idcs_X);

% X (activity, odour A cells)
this_window = p.general.bins_odour1Window;
A_2nd_AB = squeeze(nanmean(this_nft( these_idcs, this_window , find(task.odour1=="A"&task.odour2=="B") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="A"&task.odour2=="B") ),2));
A_2nd_AY = squeeze(nanmean(this_nft( these_idcs, this_window , find(task.odour1=="A"&task.odour2=="Y") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="A"&task.odour2=="Y") ),2));
A_2nd_XY = squeeze(nanmean(this_nft( these_idcs, this_window , find(task.odour1=="X"&task.odour2=="Y") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="X"&task.odour2=="Y") ),2));
A_2nd_XB = squeeze(nanmean(this_nft( these_idcs, this_window , find(task.odour1=="X"&task.odour2=="B") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="X"&task.odour2=="B") ),2));

% idcs
idcs_A = find(nem_all_cmpr.sigM_sigTG_pos{2}==1); % A
idcs_X = find(nem_all_cmpr.sigM_sigTG_pos{3}==1); % X
these_idcs = setdiff(idcs_X,idcs_A);

% X (activity, odour X cells)
this_window = p.general.bins_odour1Window;
X_2nd_AB = squeeze(nanmean(this_nft( these_idcs, this_window , find(task.odour1=="A"&task.odour2=="B") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="A"&task.odour2=="B") ),2));
X_2nd_AY = squeeze(nanmean(this_nft( these_idcs, this_window , find(task.odour1=="A"&task.odour2=="Y") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="A"&task.odour2=="Y") ),2));
X_2nd_XY = squeeze(nanmean(this_nft( these_idcs, this_window , find(task.odour1=="X"&task.odour2=="Y") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="X"&task.odour2=="Y") ),2));
X_2nd_XB = squeeze(nanmean(this_nft( these_idcs, this_window , find(task.odour1=="X"&task.odour2=="B") ),2))...
    - squeeze(nanmean(this_nft( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="X"&task.odour2=="B") ),2));

% Y (behaviour, %correct)
correct_AB = task.response(find(task.odour1=="A"&task.odour2=="B"));
correct_AB = correct_AB=="H";
correct_AY = task.response(find(task.odour1=="A"&task.odour2=="Y"));
correct_AY = correct_AY=="CR";
correct_XY = task.response(find(task.odour1=="X"&task.odour2=="Y"));
correct_XY = correct_XY=="H";
correct_XB = task.response(find(task.odour1=="X"&task.odour2=="B"));
correct_XB = correct_XB=="CR";

% Y (behaviour, presp)
presp_AB = task.response(find(task.odour1=="A"&task.odour2=="B"));
presp_AB = presp_AB=="H";
presp_AY = task.response(find(task.odour1=="A"&task.odour2=="Y"));
presp_AY = presp_AY=="FA";
presp_XY = task.response(find(task.odour1=="X"&task.odour2=="Y"));
presp_XY = presp_XY=="H";
presp_XB = task.response(find(task.odour1=="X"&task.odour2=="B"));
presp_XB = presp_XB=="FA";


%% Figure - 2nd odour

nrows = 2; ncols = 2;
F = default_figure([20,0.5,20,9.9]);

lower_x = -0.1;
upper_x = 0.25;
lower_y = -0.2;
upper_y = 1.2;

subplot(nrows,ncols,1)
hold on
[corr_r_AB,corr_p_AB] = fitLine(nanmean(B_2nd_AB,1)',correct_AB',p.col.AB);
[corr_r_XB,corr_p_XB] = fitLine(nanmean(B_2nd_XB,1)',correct_XB',p.col.XB);
[corr_r_AY,corr_p_AY] = fitLine(nanmean(B_2nd_AY,1)',correct_AY',p.col.AY);
[corr_r_XY,corr_p_XY] = fitLine(nanmean(B_2nd_XY,1)',correct_XY',p.col.XY);
xlim([lower_x,upper_x])
ylim([lower_y,upper_y])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['AB: n = ',num2str(length(B_2nd_AB(i,:))),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
    'XB: n = ',num2str(length(B_2nd_XB(i,:))),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2),newline,...
    'AY: n = ',num2str(length(B_2nd_AY(i,:))),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2),newline,...
    'XY: n = ',num2str(length(B_2nd_XY(i,:))),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
legend('AB','XB','AY','XY')
xlabel('Activity in response to 2nd odour')
ylabel('%Correct')
title('Average activity of all odour B cells')

subplot(nrows,ncols,2)
hold on
[corr_r_AB,corr_p_AB] = fitLine(nanmean(Y_2nd_AB,1)',correct_AB',p.col.AB);
[corr_r_XB,corr_p_XB] = fitLine(nanmean(Y_2nd_XB,1)',correct_XB',p.col.XB);
[corr_r_AY,corr_p_AY] = fitLine(nanmean(Y_2nd_AY,1)',correct_AY',p.col.AY);
[corr_r_XY,corr_p_XY] = fitLine(nanmean(Y_2nd_XY,1)',correct_XY',p.col.XY);
xlim([lower_x,upper_x])
ylim([lower_y,upper_y])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['AB: n = ',num2str(length(Y_2nd_AB(i,:))),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
    'XB: n = ',num2str(length(Y_2nd_XB(i,:))),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2),newline,...
    'AY: n = ',num2str(length(Y_2nd_AY(i,:))),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2),newline,...
    'XY: n = ',num2str(length(Y_2nd_XY(i,:))),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
legend('AB','XB','AY','XY')
xlabel('Activity in response to 2nd odour')
ylabel('%Correct')
title('Average activity of all odour Y cells')

subplot(nrows,ncols,3)
hold on
[corr_r_AB,corr_p_AB] = fitLine(nanmean(B_2nd_AB,1)',presp_AB',p.col.AB);
[corr_r_XB,corr_p_XB] = fitLine(nanmean(B_2nd_XB,1)',presp_XB',p.col.XB);
[corr_r_AY,corr_p_AY] = fitLine(nanmean(B_2nd_AY,1)',presp_AY',p.col.AY);
[corr_r_XY,corr_p_XY] = fitLine(nanmean(B_2nd_XY,1)',presp_XY',p.col.XY);
xlim([lower_x,upper_x])
ylim([lower_y,upper_y])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['AB: n = ',num2str(length(B_2nd_AB(i,:))),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
    'XB: n = ',num2str(length(B_2nd_XB(i,:))),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2),newline,...
    'AY: n = ',num2str(length(B_2nd_AY(i,:))),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2),newline,...
    'XY: n = ',num2str(length(B_2nd_XY(i,:))),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
legend('AB','XB','AY','XY')
xlabel('Activity in response to 2nd odour')
ylabel('P(resp)')
title('Average activity of all odour B cells')

subplot(nrows,ncols,4)
hold on
[corr_r_AB,corr_p_AB] = fitLine(nanmean(Y_2nd_AB,1)',presp_AB',p.col.AB);
[corr_r_XB,corr_p_XB] = fitLine(nanmean(Y_2nd_XB,1)',presp_XB',p.col.XB);
[corr_r_AY,corr_p_AY] = fitLine(nanmean(Y_2nd_AY,1)',presp_AY',p.col.AY);
[corr_r_XY,corr_p_XY] = fitLine(nanmean(Y_2nd_XY,1)',presp_XY',p.col.XY);
xlim([lower_x,upper_x])
ylim([lower_y,upper_y])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['AB: n = ',num2str(length(Y_2nd_AB(i,:))),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
    'XB: n = ',num2str(length(Y_2nd_XB(i,:))),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2),newline,...
    'AY: n = ',num2str(length(Y_2nd_AY(i,:))),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2),newline,...
    'XY: n = ',num2str(length(Y_2nd_XY(i,:))),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
legend('AB','XB','AY','XY')
xlabel('Activity in response to 2nd odour')
ylabel('P(resp)')
title('Average activity of all odour Y cells')

suptitle(this_suptitle)



%% Figure - 1st odour

nrows = 2; ncols = 2;
F = default_figure([20,0.5,20,9.9]);

lower_x = -0.1;
upper_x = 0.25;
lower_y = -0.2;
upper_y = 1.2;

subplot(nrows,ncols,1)
hold on
[corr_r_AB,corr_p_AB] = fitLine(nanmean(A_2nd_AB,1)',correct_AB',p.col.AB);
[corr_r_XB,corr_p_XB] = fitLine(nanmean(A_2nd_XB,1)',correct_XB',p.col.XB);
[corr_r_AY,corr_p_AY] = fitLine(nanmean(A_2nd_AY,1)',correct_AY',p.col.AY);
[corr_r_XY,corr_p_XY] = fitLine(nanmean(A_2nd_XY,1)',correct_XY',p.col.XY);
xlim([lower_x,upper_x])
ylim([lower_y,upper_y])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['AB: n = ',num2str(length(A_2nd_AB(i,:))),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
    'XB: n = ',num2str(length(A_2nd_XB(i,:))),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2),newline,...
    'AY: n = ',num2str(length(A_2nd_AY(i,:))),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2),newline,...
    'XY: n = ',num2str(length(A_2nd_XY(i,:))),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
legend('AB','XB','AY','XY')
xlabel('Activity in response to 1st odour')
ylabel('%Correct')
title('Average activity of all odour A cells')

subplot(nrows,ncols,2)
hold on
[corr_r_AB,corr_p_AB] = fitLine(nanmean(X_2nd_AB,1)',correct_AB',p.col.AB);
[corr_r_XB,corr_p_XB] = fitLine(nanmean(X_2nd_XB,1)',correct_XB',p.col.XB);
[corr_r_AY,corr_p_AY] = fitLine(nanmean(X_2nd_AY,1)',correct_AY',p.col.AY);
[corr_r_XY,corr_p_XY] = fitLine(nanmean(X_2nd_XY,1)',correct_XY',p.col.XY);
xlim([lower_x,upper_x])
ylim([lower_y,upper_y])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['AB: n = ',num2str(length(X_2nd_AB(i,:))),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
    'XB: n = ',num2str(length(X_2nd_XB(i,:))),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2),newline,...
    'AY: n = ',num2str(length(X_2nd_AY(i,:))),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2),newline,...
    'XY: n = ',num2str(length(X_2nd_XY(i,:))),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
legend('AB','XB','AY','XY')
xlabel('Activity in response to 1st odour')
ylabel('%Correct')
title('Average activity of all odour X cells')

subplot(nrows,ncols,3)
hold on
[corr_r_AB,corr_p_AB] = fitLine(nanmean(A_2nd_AB,1)',presp_AB',p.col.AB);
[corr_r_XB,corr_p_XB] = fitLine(nanmean(A_2nd_XB,1)',presp_XB',p.col.XB);
[corr_r_AY,corr_p_AY] = fitLine(nanmean(A_2nd_AY,1)',presp_AY',p.col.AY);
[corr_r_XY,corr_p_XY] = fitLine(nanmean(A_2nd_XY,1)',presp_XY',p.col.XY);
xlim([lower_x,upper_x])
ylim([lower_y,upper_y])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['AB: n = ',num2str(length(A_2nd_AB(i,:))),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
    'XB: n = ',num2str(length(A_2nd_XB(i,:))),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2),newline,...
    'AY: n = ',num2str(length(A_2nd_AY(i,:))),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2),newline,...
    'XY: n = ',num2str(length(A_2nd_XY(i,:))),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
legend('AB','XB','AY','XY')
xlabel('Activity in response to 1st odour')
ylabel('P(resp)')
title('Average activity of all odour A cells')

subplot(nrows,ncols,4)
hold on
[corr_r_AB,corr_p_AB] = fitLine(nanmean(X_2nd_AB,1)',presp_AB',p.col.AB);
[corr_r_XB,corr_p_XB] = fitLine(nanmean(X_2nd_XB,1)',presp_XB',p.col.XB);
[corr_r_AY,corr_p_AY] = fitLine(nanmean(X_2nd_AY,1)',presp_AY',p.col.AY);
[corr_r_XY,corr_p_XY] = fitLine(nanmean(X_2nd_XY,1)',presp_XY',p.col.XY);
xlim([lower_x,upper_x])
ylim([lower_y,upper_y])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['AB: n = ',num2str(length(X_2nd_AB(i,:))),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
    'XB: n = ',num2str(length(X_2nd_XB(i,:))),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2),newline,...
    'AY: n = ',num2str(length(X_2nd_AY(i,:))),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2),newline,...
    'XY: n = ',num2str(length(X_2nd_XY(i,:))),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
legend('AB','XB','AY','XY')
xlabel('Activity in response to 1st odour')
ylabel('P(resp)')
title('Average activity of all odour X cells')

suptitle(this_suptitle)


