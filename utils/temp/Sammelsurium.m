%% Sammelsurium

%popActSuppPanel()
these_avgTraces = avgTraces_all.X; this_testGroup = 'X_{sens}';
these_traces = traces_all.X;
these_trials = trials_all.stimuli.X;

% %% Preparations
% 
% these_avgTraces = these_avgTraces(:,these_bins);

%% Sort traces

these_idcs_unsorted = find(iscell==1);
%these_idcs_unsorted = these_idcs_unsorted(randperm(length(these_idcs_unsorted)));

this_testGroupIdx = find(strcmp(nem_all_cmpr.testGroup.label,this_testGroup));
these_idcs_pos = find(nem_all_cmpr.sigM_sigTG_pos{this_testGroupIdx}==1);
these_idcs_neg = find(nem_all_cmpr.sigM_sigTG_neg{this_testGroupIdx}==1);
these_idcs_other = setdiff(setdiff(these_idcs_unsorted,these_idcs_pos),these_idcs_neg);


[~,temp] = sort(nanmean(these_avgTraces(these_idcs_pos,16:25),2));
these_idcs_pos = these_idcs_pos(flip(temp));
[~,temp] = sort(nanmean(these_avgTraces(these_idcs_neg,16:25),2));
these_idcs_neg = these_idcs_neg(flip(temp));
[~,temp] = sort(nanmean(these_avgTraces(these_idcs_other,16:25),2));
these_idcs_other = these_idcs_other(flip(temp));

these_idcs_sorted = [these_idcs_pos;these_idcs_other;these_idcs_neg];


%% Normalise traces

% this_normalisationData = these_avgTraces;
% this_normalisationData_bl = these_avgTraces(:,1:5);
% 
% this_zero = nanmedian(this_normalisationData_bl,2);
% this_lower = nan(size(this_normalisationData,1),1);
% this_upper = nan(size(this_normalisationData,1),1);
% for i=1:size(this_normalisationData,1)
%     if abs(nanmin(this_normalisationData(i,:))-this_zero(i)) > abs(nanmax(this_normalisationData(i,:))-this_zero(i))
%         this_lower(i) = nanmin(this_normalisationData(i,:));
%         this_upper(i) = abs(nanmin(this_normalisationData(i,:)))+2*this_zero(i);
%     else
%         this_lower(i) = -nanmax(this_normalisationData(i,:))+2*this_zero(i);
%         this_upper(i) = nanmax(this_normalisationData(i,:));            
%     end
% end
% 
% these_normAvgTraces = nan(size(these_avgTraces,1),size(these_avgTraces,2)+3);
% for i=1:size(these_avgTraces,1)
%     these_normAvgTraces(i,:) = rescale([these_avgTraces(i,:),[this_lower(i),this_zero(i),this_upper(i)]],-1,1);
% end
% these_normAvgTraces = these_normAvgTraces(:,1:end-3);


                
%% plot              

figure;
these_normAvgTraces = these_avgTraces;
plt.clim = [-1,1]; plt.colormap = redblue;
heatMap_task(these_normAvgTraces(these_idcs_sorted,:),NaN,[],p,info,plt);


%% traces plot
      
these_data_pos = nanmean(these_avgTraces(these_idcs_pos,:),1);
these_data_pos_first = nanmean(these_avgTraces(these_idcs_pos(1:5),:),1);
these_data_neg = nanmean(these_avgTraces(these_idcs_neg,:),1);

figure;
hold on
plot(these_data_pos,'r')
plot(these_data_pos_first,'m')
plot(these_data_neg,'b')
                
                
%% bl sub plot


these_bins_bl = 1:10;
temp = nanmean(these_traces(:,these_bins_bl,:),2);
this_bl_mean = nanmean(temp,3);
this_bl_std = nanstd(temp,[],3);

figure;
these_normAvgTraces = these_avgTraces - this_bl_mean;
plt.clim = [-1,1]; plt.colormap = redblue;
heatMap_task(these_normAvgTraces(these_idcs_sorted,:),NaN,[],p,info,plt);
      


%% The more activation, the more suppression? 

these_data_pos = squeeze(nanmean(nanmean(these_traces(these_idcs_pos,16:20,:),2),1));
these_data_neg = squeeze(nanmean(nanmean(these_traces(these_idcs_neg,16:20,:),2),1));
[rho,pval] = corr(these_data_pos,these_data_neg,'Type','Pearson','Rows','Complete')

figure
scatter(these_data_pos,these_data_neg)
xlabel('Activity (z) of activated cells')
ylabel('Activity (z) of suppressed cells')
title(['each dot is a trial'])
% weird. Do baseline-subtraction!


%%

these_data_pos = squeeze(nanmean(nanmean(these_traces(these_idcs_pos,16:20,:),2) - nanmean(these_traces(these_idcs_pos,6:10,:),2),1));
these_data_neg = squeeze(nanmean(nanmean(these_traces(these_idcs_neg,16:20,:),2) - nanmean(these_traces(these_idcs_neg,6:10,:),2),1));
[rho,pval] = corr(these_data_pos,these_data_neg,'Type','Pearson','Rows','Complete')

figure
scatter(these_data_pos,these_data_neg)
xlabel('Activity (z) of activated cells')
ylabel('Activity (z) of suppressed cells')
title(['each dot is a trial'])

figure
hold on
plot(these_data_pos,'r')
plot(these_data_neg,'b')

%%

these_data_pos = squeeze(nanmean(nanmean(these_traces(these_idcs_pos,16:25,:),2) - nanmean(these_traces(these_idcs_pos,1:10,:),2),1));
these_data_neg = squeeze(nanmean(nanmean(these_traces(these_idcs_neg,16:25,:),2) - nanmean(these_traces(these_idcs_neg,1:10,:),2),1));
[rho,pval] = corr(these_data_pos,these_data_neg,'Type','Pearson','Rows','Complete')

figure
scatter(these_data_pos,these_data_neg)
xlabel('Activity (z) of activated cells')
ylabel('Activity (z) of suppressed cells')
title(['each dot is a trial'])


figure
hold on
plot(these_data_pos,'r')
plot(these_data_neg,'b')





%% avg traces plot - correct vs incorrect

these_data_neg_corr = nanmean(traces_all.X_correct(these_idcs_neg,:,:),1);
these_data_neg_incorr = nanmean(traces_all.X_incorrect(these_idcs_neg,:,:),1);

figure;
hold on 
temp=shadedErrorBar(1:size(these_data_neg_corr,2),nanmean(these_data_neg_corr,3),nansem(these_data_neg_corr,3),'lineProps','g'); temp.mainLine.LineWidth = 2; temp.mainLine.LineStyle = '-';  
temp=shadedErrorBar(1:size(these_data_neg_incorr,2),nanmean(these_data_neg_incorr,3),nansem(these_data_neg_incorr,3),'lineProps','r'); temp.mainLine.LineWidth = 2; temp.mainLine.LineStyle = '-';  


%% individual trial heat plot

figure;
plt.clim = [-3,3]; plt.colormap = redblue;
heatMap_task(these_traces(these_idcs_sorted,:,10),NaN,[],p,info,plt);



%% Trial-wise relationship between activation and suppression

this_activity = squeeze(nanmean(these_traces(these_idcs_sorted,16:20,:),2) - nanmean(these_traces(these_idcs_sorted,6:10,:),2));
these_data_pos = squeeze(nanmean(nanmean(these_traces(these_idcs_pos,16:20,:),2) - nanmean(these_traces(these_idcs_pos,6:10,:),2),1));
these_data_neg = squeeze(nanmean(nanmean(these_traces(these_idcs_neg,16:20,:),2) - nanmean(these_traces(these_idcs_neg,6:10,:),2),1));

activity_threshold = 3;
this_pos = this_activity > activity_threshold;
this_neg = this_activity < -activity_threshold;
this_num_pos = nansum(this_pos,1);
this_num_neg = nansum(this_neg,1);


figure;
hold on
plot(this_num_pos,'r')
plot(this_num_neg,'b')
xlabel('trial')
ylabel('number of cells')
legend('activated cells','suppressed cells')
title('Trial-wise relationship between activation and suppression')

% Different ways of identifying activation and suppression

figure;
scatter(this_num_pos,these_data_pos)
xlabel('Number of activated cells')
ylabel('Activation (z) of the GLM-identified activated neurons')
[rho,pval] = corr(this_num_pos',these_data_pos,'Type','Pearson','Rows','Complete');
title(['Different ways of identifying activation, rho=',num2str(rho,2),', p=',num2str(pval,2)])

figure;
scatter(this_num_neg,these_data_neg)
xlabel('Number of suppressed cells')
ylabel('Activation (z) of the GLM-identified suppressed neurons')
[rho,pval] = corr(this_num_neg',these_data_neg,'Type','Pearson','Rows','Complete');
title(['Different ways of identifying suppression, rho=',num2str(rho,2),', p=',num2str(pval,2)])

% Relationship between activation and suppression

figure;
scatter(this_num_pos,this_num_neg)
xlabel('Number of activated cells')
ylabel('Number of suppressed cells')
[rho,pval] = corr(this_num_pos',this_num_neg','Type','Pearson','Rows','Complete');
title(['Relationship between activation and suppression, rho=',num2str(rho,2),', p=',num2str(pval,2)])

figure;
scatter(these_data_pos,these_data_neg)
xlabel('Activation (z) of the GLM-identified activated neurons')
ylabel('Activation (z) of the GLM-identified suppressed neurons')
[rho,pval] = corr(these_data_pos,these_data_neg,'Type','Pearson','Rows','Complete');
title(['Relationship between activation and suppression, rho=',num2str(rho,2),', p=',num2str(pval,2)])

figure;
scatter(this_num_pos,these_data_neg)
xlabel('Number of activated cells')
ylabel('Activation (z) of the GLM-identified suppressed neurons')
[rho,pval] = corr(this_num_pos',these_data_neg,'Type','Pearson','Rows','Complete');
title(['Relationship between activation and suppression, rho=',num2str(rho,2),', p=',num2str(pval,2)])

figure;
scatter(these_data_pos,this_num_neg)
xlabel('Activation (z) of the GLM-identified activated neurons')
ylabel('Number of suppressed cells')
[rho,pval] = corr(these_data_pos,this_num_neg','Type','Pearson','Rows','Complete');
title(['Relationship between activation and suppression, rho=',num2str(rho,2),', p=',num2str(pval,2)])


%% Startle analysis

these_data_pos = nanmean(these_avgTraces(these_idcs_pos,:),1);
these_data_pos_first = nanmean(these_avgTraces(these_idcs_pos(1:5),:),1);
these_data_neg = nanmean(these_avgTraces(these_idcs_neg,:),1);

velocity = reshape(events_binned.velocity(prop.trial_frames_binned),prop.numTrials,size(prop.trial_frames_binned,1));
acceleration = reshape(events_binned.acceleration(prop.trial_frames_binned),prop.numTrials,size(prop.trial_frames_binned,1));

figure
subplot(3,1,1)
hold on
plot(these_data_pos,'r')
plot(these_data_neg,'b')
subplot(3,1,2)
temp=shadedErrorBar(1:size(velocity,2),nanmean(velocity,1),nansem(velocity,1),'lineProps','k'); temp.mainLine.LineWidth = 2; temp.mainLine.LineStyle = '-';  
subplot(3,1,3)
temp=shadedErrorBar(1:size(acceleration,2),nanmean(acceleration,1),nansem(acceleration,1),'lineProps','k'); temp.mainLine.LineWidth = 2; temp.mainLine.LineStyle = '-';  


%% Velocity vs activation/suppression

this_activity = squeeze(nanmean(these_traces(these_idcs_sorted,16:20,:),2) - nanmean(these_traces(these_idcs_sorted,6:10,:),2));
these_data_pos = squeeze(nanmean(nanmean(these_traces(these_idcs_pos,16:20,:),2) - nanmean(these_traces(these_idcs_pos,6:10,:),2),1));
these_data_neg = squeeze(nanmean(nanmean(these_traces(these_idcs_neg,16:20,:),2) - nanmean(these_traces(these_idcs_neg,6:10,:),2),1));

activity_threshold = 3;
this_pos = this_activity > activity_threshold;
this_neg = this_activity < -activity_threshold;
this_num_pos = nansum(this_pos,1);
this_num_neg = nansum(this_neg,1);

this_velocity_onset = nanmean(velocity(these_trials,6:10),2);

figure;
scatter(this_velocity_onset,these_data_pos)
xlabel('Velocity at trial onset (cm/s)')
ylabel('Activation (z) of the GLM-identified activated neurons')
[rho,pval] = corr(this_velocity_onset,these_data_pos,'Type','Pearson','Rows','Complete');
title(['Velocity and activation, rho=',num2str(rho,2),', p=',num2str(pval,2)])

figure;
scatter(this_velocity_onset,these_data_neg)
xlabel('Velocity at trial onset (cm/s)')
ylabel('Activation (z) of the GLM-identified suppressed neurons')
[rho,pval] = corr(this_velocity_onset,these_data_neg,'Type','Pearson','Rows','Complete');
title(['Velocity and suppression, rho=',num2str(rho,2),', p=',num2str(pval,2)])

figure;
scatter(this_velocity_onset,this_num_pos)
xlabel('Velocity at trial onset (cm/s)')
ylabel('Number of activated cells')
[rho,pval] = corr(this_velocity_onset,this_num_pos','Type','Pearson','Rows','Complete');
title(['Velocity and activation, rho=',num2str(rho,2),', p=',num2str(pval,2)])

figure;
scatter(this_velocity_onset,this_num_neg)
xlabel('Velocity at trial onset (cm/s)')
ylabel('Number of suppressed cells')
[rho,pval] = corr(this_velocity_onset,this_num_neg','Type','Pearson','Rows','Complete');
title(['Velocity and suppression, rho=',num2str(rho,2),', p=',num2str(pval,2)])


%% Velocity vs activation/suppression

this_activity = squeeze(nanmean(these_traces(these_idcs_sorted,16:20,:),2) - nanmean(these_traces(these_idcs_sorted,6:10,:),2));
these_data_pos = squeeze(nanmean(nanmean(these_traces(these_idcs_pos,16:20,:),2) - nanmean(these_traces(these_idcs_pos,6:10,:),2),1));
these_data_neg = squeeze(nanmean(nanmean(these_traces(these_idcs_neg,16:20,:),2) - nanmean(these_traces(these_idcs_neg,6:10,:),2),1));

activity_threshold = 3;
this_pos = this_activity > activity_threshold;
this_neg = this_activity < -activity_threshold;
this_num_pos = nansum(this_pos,1);
this_num_neg = nansum(this_neg,1);

this_velocity_onset = nanmean(velocity(these_trials,6:10),2);
this_velocity_stim = nanmean(velocity(these_trials,16:20),2);
this_velocity_delta = this_velocity_stim - this_velocity_onset;

figure;
scatter(this_velocity_delta,these_data_pos)
xlabel('Velocity delta (cm/s)')
ylabel('Activation (z) of the GLM-identified activated neurons')
[rho,pval] = corr(this_velocity_delta,these_data_pos,'Type','Pearson','Rows','Complete');
title(['Velocity and activation, rho=',num2str(rho,2),', p=',num2str(pval,2)])

figure;
scatter(this_velocity_delta,these_data_neg)
xlabel('Velocity delta (cm/s)')
ylabel('Activation (z) of the GLM-identified suppressed neurons')
[rho,pval] = corr(this_velocity_delta,these_data_neg,'Type','Pearson','Rows','Complete');
title(['Velocity and suppression, rho=',num2str(rho,2),', p=',num2str(pval,2)])

figure;
scatter(this_velocity_delta,this_num_pos)
xlabel('Velocity delta (cm/s)')
ylabel('Number of activated cells')
[rho,pval] = corr(this_velocity_delta,this_num_pos','Type','Pearson','Rows','Complete');
title(['Velocity and activation, rho=',num2str(rho,2),', p=',num2str(pval,2)])

figure;
scatter(this_velocity_delta,this_num_neg)
xlabel('Velocity delta (cm/s)')
ylabel('Number of suppressed cells')
[rho,pval] = corr(this_velocity_delta,this_num_neg','Type','Pearson','Rows','Complete');
title(['Velocity and suppression, rho=',num2str(rho,2),', p=',num2str(pval,2)])

