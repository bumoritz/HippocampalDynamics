%% B_pos neurons

nrows = 2; ncols = 2;
F = default_figure([20,0.5,20,9.9]);

idcs_B = find(nem_all_cmpr.sigM_sigTG_pos{11}==1); % B
idcs_Y = find(nem_all_cmpr.sigM_sigTG_pos{12}==1); % Y

these_idcs = find(nem_all_cmpr.sigM_sigTG{13}==1);   %setdiff(idcs_B,idcs_Y);

conditions = {'AB','AY','XY','XB'};
for i=1:4
    subplot(nrows,ncols,i)
    hold on

    if strcmp(conditions{i},'AB')
        temp = nanmean(traces_all.A_M(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','r'); temp1.mainLine.LineWidth = 2;  
        temp = nanmean(traces_all.A_H(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','g'); temp1.mainLine.LineWidth = 2;  
    elseif strcmp(conditions{i},'AY')
        temp = nanmean(traces_all.A_FA(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','r'); temp1.mainLine.LineWidth = 2;  
        temp = nanmean(traces_all.A_CR(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','g'); temp1.mainLine.LineWidth = 2;  
    elseif strcmp(conditions{i},'XY')
        temp = nanmean(traces_all.X_M(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','r'); temp1.mainLine.LineWidth = 2;  
        temp = nanmean(traces_all.X_H(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','g'); temp1.mainLine.LineWidth = 2;  
    elseif strcmp(conditions{i},'XB')
        temp = nanmean(traces_all.X_FA(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','r'); temp1.mainLine.LineWidth = 2;  
        temp = nanmean(traces_all.X_CR(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','g'); temp1.mainLine.LineWidth = 2;  
    end
    
    temp = nanmean(traces_all.(conditions{i})(these_idcs,:,:),3);
    temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps',p.col.black); temp1.mainLine.LineWidth = 2;   
    plt = struct(); traces_task([],p,info,plt);
    temp2(i,:) = get(gca,'ylim');
    ylabel('Activity (z-score)')
    title(conditions{i})
end
for i=1:4
    subplot(nrows,ncols,i)
    ylim([nanmin(temp2(:)),nanmax(temp2(:))])
end

suptitle('Average of all positive/negative int neurons')


%% B_pos neurons (residuals)

% nemd_all.data
% nemd_all.full
% nemd_all.testGroupOut{11} % B
% nemd_all.testGroupResidual{11} % B
% nemd_all.data - nemd_all.testGroupOut{11};
this_data = nemd_all.data;
temp = reshape(this_data(nem_all.prop.trial_frames_binned,:)',[],size(nem_all.prop.trial_frames_binned,1),nem_all.prop.numTrials);
[this_traces_struct,~,~] = createTracesStructs(temp,trials_all);
% above: this_traces_struct = traces_all;

nrows = 2; ncols = 2;
F = default_figure([20,0.5,20,9.9]);

idcs_B = find(nem_all_cmpr.sigM_sigTG_pos{11}==1); % B
idcs_Y = find(nem_all_cmpr.sigM_sigTG_pos{12}==1); % Y
these_idcs = setdiff(idcs_B,idcs_Y);

conditions = {'AB','AY','XY','XB'};
for i=1:4
    subplot(nrows,ncols,i)
    hold on

    if strcmp(conditions{i},'AB')
        temp = nanmean(this_traces_struct.A_M(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','r'); temp1.mainLine.LineWidth = 2;  
        temp = nanmean(this_traces_struct.A_H(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','g'); temp1.mainLine.LineWidth = 2;  
    elseif strcmp(conditions{i},'AY')
        temp = nanmean(this_traces_struct.A_FA(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','r'); temp1.mainLine.LineWidth = 2;  
        temp = nanmean(this_traces_struct.A_CR(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','g'); temp1.mainLine.LineWidth = 2;  
    elseif strcmp(conditions{i},'XY')
        temp = nanmean(this_traces_struct.X_M(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','r'); temp1.mainLine.LineWidth = 2;  
        temp = nanmean(this_traces_struct.X_H(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','g'); temp1.mainLine.LineWidth = 2;  
    elseif strcmp(conditions{i},'XB')
        temp = nanmean(this_traces_struct.X_FA(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','r'); temp1.mainLine.LineWidth = 2;  
        temp = nanmean(this_traces_struct.X_CR(these_idcs,:,:),3);
        temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps','g'); temp1.mainLine.LineWidth = 2;  
    end
    
    temp = nanmean(this_traces_struct.(conditions{i})(these_idcs,:,:),3);
    temp1=shadedErrorBar(1:size(temp,2),nanmean(temp,1),nansem(temp,1),'lineProps',p.col.black); temp1.mainLine.LineWidth = 2;   
    plt = struct(); traces_task([],p,info,plt);
    temp2(i,:) = get(gca,'ylim');
    ylabel('Activity (z-score)')
    title(conditions{i})
end
for i=1:4
    subplot(nrows,ncols,i)
    ylim([nanmin(temp2(:)),nanmax(temp2(:))])
end

suptitle('Average of all positive B (and not Y) neurons')


%%

nrows = 2; ncols = 2;
F = default_figure([20,0.5,20,9.9]);





%% I: First odour integration selectivity index
% look at second odour time bins and compare activity between A and X trials (averaged over second odour identity)

% maybe better: input residuals instead of actual activity into these analyses

% e.g. take average of all odour B trials

% both as function of trial block number
% and as function of performance (trial-block-wise)
% and correlation with binary correct or incorrect (trial-wise)


% int = 13;
idx = 2501;

figure;
plot(avgTraces_all.B(idx,:))
hold on
plot(avgTraces_all.AB(idx,:),':')
plot(avgTraces_all.XB(idx,:),':')

% Here: AB higher than XB
% Hypothesis: the higher the B response, the higher P(resp)
% also within trial type?

% deciles, quartiles

%osip.i.


%%

int_2nd_AB = squeeze(nanmean(nft_binned( find(nem_all_cmpr.sigM_sigTG{13}==1), p.general.bins_odour2Window , find(task.odour1=="A"&task.odour2=="B") ),2))...
    - squeeze(nanmean(nft_binned( find(nem_all_cmpr.sigM_sigTG{13}==1), p.general.bins_baselineWindow , find(task.odour1=="A"&task.odour2=="B") ),2));
int_2nd_AY = squeeze(nanmean(nft_binned( find(nem_all_cmpr.sigM_sigTG{13}==1), p.general.bins_odour2Window , find(task.odour1=="A"&task.odour2=="Y") ),2))...
    - squeeze(nanmean(nft_binned( find(nem_all_cmpr.sigM_sigTG{13}==1), p.general.bins_baselineWindow , find(task.odour1=="A"&task.odour2=="Y") ),2));
int_2nd_XY = squeeze(nanmean(nft_binned( find(nem_all_cmpr.sigM_sigTG{13}==1), p.general.bins_odour2Window , find(task.odour1=="X"&task.odour2=="Y") ),2))...
    - squeeze(nanmean(nft_binned( find(nem_all_cmpr.sigM_sigTG{13}==1), p.general.bins_baselineWindow , find(task.odour1=="X"&task.odour2=="Y") ),2));
int_2nd_XB = squeeze(nanmean(nft_binned( find(nem_all_cmpr.sigM_sigTG{13}==1), p.general.bins_odour2Window , find(task.odour1=="X"&task.odour2=="B") ),2))...
    - squeeze(nanmean(nft_binned( find(nem_all_cmpr.sigM_sigTG{13}==1), p.general.bins_baselineWindow , find(task.odour1=="X"&task.odour2=="B") ),2));

int_2nd_AB_ranks = nan(size(int_2nd_AB));
int_2nd_AY_ranks = nan(size(int_2nd_AY));
int_2nd_XY_ranks = nan(size(int_2nd_XY));
int_2nd_XB_ranks = nan(size(int_2nd_XB));
for i=1:size(int_2nd_AB,1)
    int_2nd_AB_ranks(i,:) = quantileranks(int_2nd_AB(i,:),10);
    int_2nd_AY_ranks(i,:) = quantileranks(int_2nd_AY(i,:),10);
    int_2nd_XY_ranks(i,:) = quantileranks(int_2nd_XY(i,:),10);
    int_2nd_XB_ranks(i,:) = quantileranks(int_2nd_XB(i,:),10);
end

idcs_B = find(nem_all_cmpr.sigM_sigTG_pos{11}==1); % B
idcs_Y = find(nem_all_cmpr.sigM_sigTG_pos{12}==1); % Y
these_idcs = setdiff(idcs_B,idcs_Y);
B_2nd_AB = squeeze(nanmean(nft_binned( these_idcs, p.general.bins_odour2Window , find(task.odour1=="A"&task.odour2=="B") ),2))...
    - squeeze(nanmean(nft_binned( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="A"&task.odour2=="B") ),2));
B_2nd_AY = squeeze(nanmean(nft_binned( these_idcs, p.general.bins_odour2Window , find(task.odour1=="A"&task.odour2=="Y") ),2))...
    - squeeze(nanmean(nft_binned( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="A"&task.odour2=="Y") ),2));
B_2nd_XY = squeeze(nanmean(nft_binned( these_idcs, p.general.bins_odour2Window , find(task.odour1=="X"&task.odour2=="Y") ),2))...
    - squeeze(nanmean(nft_binned( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="X"&task.odour2=="Y") ),2));
B_2nd_XB = squeeze(nanmean(nft_binned( these_idcs, p.general.bins_odour2Window , find(task.odour1=="X"&task.odour2=="B") ),2))...
    - squeeze(nanmean(nft_binned( these_idcs, p.general.bins_baselineWindow , find(task.odour1=="X"&task.odour2=="B") ),2));

B_2nd_AB_ranks = nan(size(B_2nd_AB));
B_2nd_AY_ranks = nan(size(B_2nd_AY));
B_2nd_XY_ranks = nan(size(B_2nd_XY));
B_2nd_XB_ranks = nan(size(B_2nd_XB));
for i=1:size(B_2nd_AB,1)
    B_2nd_AB_ranks(i,:) = quantileranks(B_2nd_AB(i,:),10);
    B_2nd_AY_ranks(i,:) = quantileranks(B_2nd_AY(i,:),10);
    B_2nd_XY_ranks(i,:) = quantileranks(B_2nd_XY(i,:),10);
    B_2nd_XB_ranks(i,:) = quantileranks(B_2nd_XB(i,:),10);
end

presp_AB = task.response(find(task.odour1=="A"&task.odour2=="B"));
presp_AB = presp_AB=="H";
presp_AY = task.response(find(task.odour1=="A"&task.odour2=="Y"));
presp_AY = presp_AY=="FA";
presp_XY = task.response(find(task.odour1=="X"&task.odour2=="Y"));
presp_XY = presp_XY=="H";
presp_XB = task.response(find(task.odour1=="X"&task.odour2=="B"));
presp_XB = presp_XB=="FA";


%%

i = 73;

population = find(nem_all_cmpr.sigM_sigTG{13}==1);
[a,b,c]=intersect(population,[16,42,138,171,205,236,248,268,282,294,369,432,480,483,487,568,583,653,679]);

% for k=1:length(b) %1:73
%     i=b(k);
    
    nrows = 2; ncols = 2;
    F = default_figure([20,0.5,20,9.9]);

    % a)
    this_data_x_A = nan(1,10);
    this_data_x_X = nan(1,10);
    this_data_y_A = nan(1,10);
    this_data_y_X = nan(1,10);
    for j=1:10
        this_data_x_A(j) = nanmean(int_2nd_AB(i,int_2nd_AB_ranks(i,:)==j));
        this_data_x_X(j) = nanmean(int_2nd_XB(i,int_2nd_XB_ranks(i,:)==j));
        this_data_y_A(j) = nanmean(presp_AB(int_2nd_AB_ranks(i,:)==j));
        this_data_y_X(j) = nanmean(presp_XB(int_2nd_XB_ranks(i,:)==j));
    end

    subplot(nrows,ncols,1)
    hold on
    plot(this_data_x_A,this_data_y_A,'Color',p.col.A)
    plot(this_data_x_X,this_data_y_X,'Color',p.col.X)
    legend('AB','XB')
    xlabel('Activity in response to odour B')
    ylabel('P(response)')
    title('trials split into activity deciles')

    subplot(nrows,ncols,3)
    hold on
    [corr_r_AB,corr_p_AB] = fitLine(int_2nd_AB(i,:)',presp_AB',p.col.A);
    [corr_r_XB,corr_p_XB] = fitLine(int_2nd_XB(i,:)',presp_XB',p.col.X);
    scatter(int_2nd_AB(i,:),presp_AB,'CData',p.col.A);
    scatter(int_2nd_XB(i,:),presp_XB,'CData',p.col.X);
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),...
        ['AB: n = ',num2str(length(int_2nd_AB(i,:))),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
        'XB: n = ',num2str(length(int_2nd_XB(i,:))),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    xlabel('Activity in response to odour B')
    ylabel('P(response)')
    title('single trials')

    % a)
    this_data_x_A = nan(1,10);
    this_data_x_X = nan(1,10);
    this_data_y_A = nan(1,10);
    this_data_y_X = nan(1,10);
    for j=1:10
        this_data_x_A(j) = nanmean(int_2nd_AY(i,int_2nd_AY_ranks(i,:)==j));
        this_data_x_X(j) = nanmean(int_2nd_XY(i,int_2nd_XY_ranks(i,:)==j));
        this_data_y_A(j) = nanmean(presp_AY(int_2nd_AY_ranks(i,:)==j));
        this_data_y_X(j) = nanmean(presp_XY(int_2nd_XY_ranks(i,:)==j));
    end

    subplot(nrows,ncols,2)
    hold on
    plot(this_data_x_A,this_data_y_A,'Color',p.col.A)
    plot(this_data_x_X,this_data_y_X,'Color',p.col.X)
    legend('AY','XY')
    xlabel('Activity in response to odour Y')
    ylabel('P(response)')
    title('trials split into activity deciles')

    subplot(nrows,ncols,4)
    hold on
    [corr_r_AY,corr_p_AY] = fitLine(int_2nd_AY(i,:)',presp_AY',p.col.A);
    [corr_r_XY,corr_p_XY] = fitLine(int_2nd_XY(i,:)',presp_XY',p.col.X);
    scatter(int_2nd_AY(i,:),presp_AY,'CData',p.col.A);
    scatter(int_2nd_XY(i,:),presp_XY,'CData',p.col.X);
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),...
        ['XY: n = ',num2str(length(int_2nd_XY(i,:))),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2),newline,...
        'AY: n = ',num2str(length(int_2nd_AY(i,:))),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    xlabel('Activity in response to odour Y')
    ylabel('P(response)')
    title('single trials')
    
    this_population = find(nem_all_cmpr.sigM_sigTG{13}==1);
    suptitle(num2str(this_population(i)))
    
%     saveas(F,['C:\SniffinHippo\Temp\',num2str(i),'.png'])
%     drawnow;
% end






%% --------------------------------

% population = find(nem_all_cmpr.sigM_sigTG{13}==1);
% [a,b,c]=intersect(population,[16,42,138,171,205,236,248,268,282,294,369,432,480,483,487,568,583,653,679]);

    
nrows = 1; ncols = 2;
F = default_figure([20,0.5,20,9.9]);

% a)
this_data_x_AB = nan(length(these_idcs),10);
this_data_x_XB = nan(length(these_idcs),10);
this_data_y_AB = nan(length(these_idcs),10);
this_data_y_XB = nan(length(these_idcs),10);
this_data_x_AY = nan(length(these_idcs),10);
this_data_x_XY = nan(length(these_idcs),10);
this_data_y_AY = nan(length(these_idcs),10);
this_data_y_XY = nan(length(these_idcs),10);
for i=1:length(these_idcs)
    for j=1:10
        this_data_x_AB(i,j) = nanmean(B_2nd_AB(i,B_2nd_AB_ranks(i,:)==j));
        this_data_x_XB(i,j) = nanmean(B_2nd_XB(i,B_2nd_XB_ranks(i,:)==j));
        this_data_y_AB(i,j) = nanmean(presp_AB(B_2nd_AB_ranks(i,:)==j));
        this_data_y_XB(i,j) = nanmean(presp_XB(B_2nd_XB_ranks(i,:)==j));
        this_data_x_AY(i,j) = nanmean(B_2nd_AY(i,B_2nd_AY_ranks(i,:)==j));
        this_data_x_XY(i,j) = nanmean(B_2nd_XY(i,B_2nd_XY_ranks(i,:)==j));
        this_data_y_AY(i,j) = nanmean(presp_AY(B_2nd_AY_ranks(i,:)==j));
        this_data_y_XY(i,j) = nanmean(presp_XY(B_2nd_XY_ranks(i,:)==j));
    end
end

subplot(nrows,ncols,1)
hold on
plot(nanmean(this_data_x_AB,1),nanmean(this_data_y_AB,1),'Color',p.col.AB)
plot(nanmean(this_data_x_XB,1),nanmean(this_data_y_XB,1),'Color',p.col.XB)
plot(nanmean(this_data_x_AY,1),nanmean(this_data_y_AY,1),'Color',p.col.AY)
plot(nanmean(this_data_x_XY,1),nanmean(this_data_y_XY,1),'Color',p.col.XY)
legend('AB','XB','AY','XY')
xlabel('Activity in response to odour B')
ylabel('P(response)')
title('trials split into activity deciles')

subplot(nrows,ncols,2)
hold on
[corr_r_AB,corr_p_AB] = fitLine(nanmean(B_2nd_AB,1)',presp_AB',p.col.AB);
[corr_r_XB,corr_p_XB] = fitLine(nanmean(B_2nd_XB,1)',presp_XB',p.col.XB);
[corr_r_AY,corr_p_AY] = fitLine(nanmean(B_2nd_AY,1)',presp_AY',p.col.AY);
[corr_r_XY,corr_p_XY] = fitLine(nanmean(B_2nd_XY,1)',presp_XY',p.col.XY);
% scatter(B_2nd_AB(i,:),presp_AB,'CData',p.col.AB);
% scatter(B_2nd_XB(i,:),presp_XB,'CData',p.col.XB);
% scatter(B_2nd_AY(i,:),presp_AY,'CData',p.col.AY);
% scatter(B_2nd_XY(i,:),presp_XY,'CData',p.col.XY);
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['AB: n = ',num2str(length(B_2nd_AB(i,:))),', \rho = ',num2str(corr_r_AB,2),', p = ',num2str(corr_p_AB,2),newline,...
    'XB: n = ',num2str(length(B_2nd_XB(i,:))),', \rho = ',num2str(corr_r_XB,2),', p = ',num2str(corr_p_XB,2),newline,...
    'AY: n = ',num2str(length(B_2nd_AY(i,:))),', \rho = ',num2str(corr_r_AY,2),', p = ',num2str(corr_p_AY,2),newline,...
    'XY: n = ',num2str(length(B_2nd_XY(i,:))),', \rho = ',num2str(corr_r_XY,2),', p = ',num2str(corr_p_XY,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
xlabel('Activity in response to odour B')
ylabel('P(response)')
title('single trials')

suptitle('positive B (and not Y) neurons')







