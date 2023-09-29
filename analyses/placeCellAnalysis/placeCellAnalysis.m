function [pca] = placeCellAnalysis(info,iscell,ops,p,path,task,sync_beh,act,tng_all,bcon)
% p = get_p; sync_beh = paq_beh; 

%% Preparations

disp('--- Preparations')

% pre-process activity measure
[prop,~,nft_binned] = preprocessActivityMeasure(act,p.tng,p,sync_beh,iscell);

% create trials, traces, avgTraces and normAvgTraces structs
trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
% [traces_all,~,~] = createTracesStructs(nft_binned,trials_all);

% pre-process distance measure
[events_binned] = binPaqEvents(sync_beh,task,p,prop);
try
    pca.distance = reshape(events_binned.distance(prop.trial_frames_binned'),size(prop.trial_frames_binned,2),size(prop.trial_frames_binned,1));
catch
    pca.distance = reshape(events_binned.distance(prop.trial_frames_binned(:,1:end-1)'),size(prop.trial_frames_binned,2)-1,size(prop.trial_frames_binned,1));
end
pca.distance_trunc = pca.distance;
pca.distance_trunc(:,42:end) = NaN;


%% Make a place map

pca.p_bin_size = 5;
pca.p_occupancy = 10;
pca.p_smoothing = 5;

% bin activity by position
pca.p_edges = 0:pca.p_bin_size:2000; %0:20:2000; % cm
temp = discretize(pca.distance_trunc,pca.p_edges);
nft_space = nan(size(nft_binned,1),length(pca.p_edges)-1,size(nft_binned,3));
for k=1:size(nft_binned,3)
    for j=1:length(pca.p_edges)-1
        nft_space(:,j,k) = nanmean(nft_binned(:,find(temp(k,:)==j),k),2);
    end
end

% smoothing
if pca.p_smoothing ~=0
    nft_space_smo = smoothdata(nft_space,2,'gaussian',pca.p_smoothing);
else
    nft_space_smo = nft_space;
end

% make average map (only at bins with sufficient occupancy)
nft_space_A = nft_space_smo(:,:,trials_all.stimuli.A);
nft_space_X = nft_space_smo(:,:,trials_all.stimuli.X);
pca.nft_space_A_avg = nan(size(nft_space_A,1),size(nft_space_A,2));
pca.nft_space_X_avg = nan(size(nft_space_X,1),size(nft_space_X,2));
for j=1:length(pca.p_edges)-1
    if max(sum(~isnan(squeeze(nft_space_A(:,j,:))),2)) > pca.p_occupancy*2.5 %125 %50 % 25 is 10% occupancy
        pca.nft_space_A_avg(:,j) = nanmean(nft_space_A(:,j,:),3);
    end
    if max(sum(~isnan(squeeze(nft_space_X(:,j,:))),2)) > pca.p_occupancy*2.5 %125 %50 % 25 is 10% occupancy
        pca.nft_space_X_avg(:,j) = nanmean(nft_space_X(:,j,:),3);
    end
end


%% Plot rate map (with time on x-axis, like for Fig2_ExampleSnakePlots_1tile)

% load data
inputType = 'tng'; %'tngn';
% act = getActivityMeasure(p.(inputType),act_struct);
load([path.filepart_out,[inputType,'_all.mat']]);
this_struct = eval([inputType,'_all']);

% process data
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.(inputType),p,sync_beh,iscell);
bl_trial_frames_binned = [prop.trial_frames_binned(1,:)-[1:5]';prop.trial_frames_binned(1:10,:)];
nft_bl_binned = reshape(nf_binned(:,bl_trial_frames_binned),[],size(bl_trial_frames_binned,1),prop.numTrials);
nft_bl_avg_binned = nanmean(nft_bl_binned,3);
[traces_all,avgTraces_all,~] = createTracesStructs(nft_binned,trials_all);

% pool data from A and X
these_refs_sort = 1; 
this_tracesStruct = traces_all; these_idcs_unsorted = {find(this_struct.passed.AW.Aonly_sel==1),find(this_struct.passed.AW.Xonly_sel==1)}; these_idcs_sorted = []; these_refs_sort = [1;4]; these_refs_norm = [];
these_avgTraces = {}; these_avgTraces{1} = avgTraces_all.A; these_avgTraces{2} = avgTraces_all.X;
this_data_A = these_avgTraces{1}(these_idcs_unsorted{1},:);
this_data_X = these_avgTraces{2}(these_idcs_unsorted{2},:);
this_data = [this_data_A;this_data_X];

% normalise data
this_normalisationData = this_data(:,p.general.bins_analysisWindow);
this_normalisationData_bl = this_data(:,p.general.bins_baselineWindow);
this_lower = nanmean(nft_bl_avg_binned([these_idcs_unsorted{1};these_idcs_unsorted{2}],:),2);
this_upper = nanmax(this_normalisationData,[],2);
this_data_norm = (this_data-this_lower) ./ (this_upper-this_lower);

% sort data
[~,temp] = nanmax(this_data_norm(:,p.general.bins_analysisWindow),[],2);
[~,this_order] = sort(temp);
this_order = flip(this_order);

% plot data
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; heatMap_task(this_data_norm(this_order,:),NaN,[],p,info,plt);
xlim([6,53])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off

savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','placeCellAnalysis_t.fig']);
saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','placeCellAnalysis_t.png']);


%% Plot rate map (with distance on x-axis, equivalent to time plot)
% always run right after time on x-axis plot

% select data to plot (and keep old indexing and order)
this_data_A = pca.nft_space_A_avg(these_idcs_unsorted{1},:); %pca_ub.pca.nft_space_A_avg(these_idcs_unsorted{1},:); %pca.nft_space_A_avg(these_idcs_unsorted{1},:);
this_data_X = pca.nft_space_X_avg(these_idcs_unsorted{2},:); %pca_ub.pca.nft_space_X_avg(these_idcs_unsorted{2},:); %pca.nft_space_X_avg(these_idcs_unsorted{2},:);
this_data = [this_data_A;this_data_X];
this_data = rmmissing(this_data,2);

% normalise data
% this_normalisationData = this_data(:,p.general.bins_analysisWindow);
% this_lower = nanmean(nft_bl_avg_binned([these_idcs_unsorted{1};these_idcs_unsorted{2}],:),2);
% this_upper = nanmax(this_normalisationData,[],2);
this_upper = nanmax(this_data,[],2);
this_data_norm = (this_data-this_lower) ./ (this_upper-this_lower);

% plot data
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; %heatMap_task(this_data_norm(this_order,:),NaN,[],p,info,plt);

% work around heatMap_task
data = this_data_norm(this_order,:);
imagesc(data);
xlim([0,size(data,2)]+0.5);
ylim([0,size(data,1)]+0.5);
set(gca,'CLim',plt.clim);
colormap(gca,'parula');

title([num2str(pca.p_bin_size),' cm bins, ',num2str(pca.p_occupancy),' occ., ',num2str(pca.p_smoothing),' smoothing, on pre-proc.'])

% make plot nice
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off

savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','placeCellAnalysis_x_tsort.fig']);
saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','placeCellAnalysis_x_tsort.png']);


%%

[~,pca.x_peakBin] = nanmax(this_data_norm,[],2);
[~,this_order] = sort(pca.x_peakBin);
this_order = flip(this_order);
pca.x_peakBin_cm = pca.p_edges(pca.x_peakBin);

% plot data
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; %heatMap_task(this_data_norm(this_order,:),NaN,[],p,info,plt);

% work around heatMap_task
data = this_data_norm(this_order,:);
imagesc(data);
xlim([0,size(data,2)]+0.5);
ylim([0,size(data,1)]+0.5);
set(gca,'CLim',plt.clim);
colormap(gca,'parula');

title([num2str(pca.p_bin_size),' cm bins, ',num2str(pca.p_occupancy),' occ., ',num2str(pca.p_smoothing),' smoothing, on pre-proc.'])

% make plot nice
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off

savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','placeCellAnalysis_x_xsort.fig']);
saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','placeCellAnalysis_x_xsort.png']);


%% Same for unbinned

% % get unbinned distances
% distance_raw = nan(size(sync_beh.position));
% position_at_trial_start = NaN;
% for i=1:length(distance_raw)
%     if ismember(i,sync_beh.sync) %sync_beh.sync(i)==1
%         position_at_trial_start = sync_beh.position(i);
%         distance_raw(i-15*p.general.binSize:i-1) = sync_beh.position(i-15*p.general.binSize:i-1) - position_at_trial_start;
%     end
%     distance_raw(i) = sync_beh.position(i)-position_at_trial_start;
% end
% 
% % identify trial bounds
% these_trial_bounds(1,:) = sync_beh.sync-length(p.general.bins_pre)*p.general.binSize;
% these_trial_bounds(2,:) = these_trial_bounds(1,:)+p.general.numBins*p.general.binSize-1;
% for i=1:length(sync_beh.sync)
%     pca_ub.trial_frames(:,i) = these_trial_bounds(1,i) : these_trial_bounds(2,i);
% end
% 
% % split unbinned distances into trials
% try
%     pca_ub.distance = reshape(distance_raw(pca_ub.trial_frames'),size(prop.trial_frames_binned,2),size(prop.trial_frames_binned,1)*p.general.binSize);
% catch
%     pca_ub.distance = reshape(distance_raw(pca_ub.trial_frames(:,1:end-1)'),size(prop.trial_frames_binned,2)-1,size(prop.trial_frames_binned,1)*p.general.binSize);
% end
% pca_ub.distance_trunc = pca_ub.distance;
% pca_ub.distance_trunc(:,41*p.general.binSize+1:end) = NaN;
% 
% % get nft_unbinned
% nft_unbinned = reshape(act(:,pca_ub.trial_frames),[],size(pca_ub.trial_frames,1),prop.numTrials);
% 
% % bin activity by position
% pca_ub.pca.p_edges = 0:20:2000; % cm
% temp = discretize(pca_ub.distance_trunc,pca_ub.pca.p_edges);
% nft_space = nan(size(nft_unbinned,1),length(pca_ub.pca.p_edges)-1,size(nft_unbinned,3));
% for k=1:size(nft_unbinned,3)
%     for j=1:length(pca_ub.pca.p_edges)-1
%         pca_ub.nft_space(:,j,k) = nanmean(nft_unbinned(:,find(temp(k,:)==j),k),2);
%     end
% end
% 
% % make average map (only at bins with sufficient occupancy)
% pca_ub.nft_space_A = pca_ub.nft_space(:,:,trials_all.stimuli.A);
% pca_ub.nft_space_X = pca_ub.nft_space(:,:,trials_all.stimuli.X);
% pca_ub.pca.nft_space_A_avg = nan(size(pca_ub.nft_space_A,1),size(pca_ub.nft_space_A,2));
% pca_ub.pca.nft_space_X_avg = nan(size(pca_ub.nft_space_X,1),size(pca_ub.nft_space_X,2));
% for j=1:length(pca_ub.pca.p_edges)-1
%     if max(sum(~isnan(squeeze(pca_ub.nft_space_A(:,j,:))),2)) > 125 %50 % 25 is 10% occupancy
%         pca_ub.pca.nft_space_A_avg(:,j) = nanmean(pca_ub.nft_space_A(:,j,:),3);
%     end
%     if max(sum(~isnan(squeeze(pca_ub.nft_space_X(:,j,:))),2)) > 125 %50 % 25 is 10% occupancy
%         pca_ub.pca.nft_space_X_avg(:,j) = nanmean(pca_ub.nft_space_X(:,j,:),3);
%     end
% end


%% Calculate spatial metrics

% % active trials
% pca.trials.A = trials_all.stimuli.A;
% pca.trials.X = trials_all.stimuli.X;
% pca.trials.active_A = tng_all.firingField.A_AW.activeTrials_ipsi;
% pca.trials.active_X = tng_all.firingField.X_AW.activeTrials_ipsi;
% 
% % speed
% pca.speed.sess = bcon.trialwise_full.velocity;
% pca.speed.AW = bcon.trialwise_AW.velocity;
% pca.speed.AW_A = pca.speed.AW(pca.trials.A);
% pca.speed.AW_X = pca.speed.AW(pca.trials.X);
% 
% % all A trials
% pca.actAtPeak_allTrials_A = nan(size(pca.trials.active_A));
% pca.binAtPeak_allTrials_A = nan(size(pca.trials.active_A));
% pca.timeAtPeak_allTrials_A = nan(size(pca.trials.active_A));
% pca.placeAtPeak_allTrials_A = nan(size(pca.trials.active_A));
% pca.speed.AW_allTrials_A = nan(size(pca.trials.active_A));
% pca.corr_speed_timeAtPeak_allTrials_A = nan(size(pca.trials.active_A,1),2);
% for i=1:length(iscell)
%     if iscell(i)==1
%         this_data = traces_all.A(i,p.general.bins_analysisWindow,:);
%         if any(~isnan(this_data(:)))
%             [temp1,temp2] = nanmax(this_data,[],2);
%             pca.actAtPeak_allTrials_A(i,1:length(temp1)) = squeeze(temp1);
%             pca.binAtPeak_allTrials_A(i,1:length(temp2)) = squeeze(temp2)+p.general.bins_analysisWindow(1)-1;
%             pca.timeAtPeak_allTrials_A(i,1:length(temp2)) = p.general.t_binned(squeeze(temp2)+p.general.bins_analysisWindow(1)-1);
%             for k=1:length(pca.trials.A)
%                 pca.placeAtPeak_allTrials_A(i,k) = pca.distance(pca.trials.A(k),pca.binAtPeak_allTrials_A(i,k));
%             end
%             pca.speed.AW_allTrials_A(i,:) = pca.speed.AW(pca.trials.A);
%             [this_rho,this_pval] = corr(pca.speed.AW_allTrials_A(i,:)',pca.timeAtPeak_allTrials_A(i,:)','Type','Pearson','Rows','Complete');
%             pca.corr_speed_timeAtPeak_allTrials_A(i,:) = [this_rho,this_pval];
%         end
%     end
% end
% 
% % all X trials
% pca.actAtPeak_allTrials_X = nan(size(pca.trials.active_X));
% pca.binAtPeak_allTrials_X = nan(size(pca.trials.active_X));
% pca.timeAtPeak_allTrials_X = nan(size(pca.trials.active_X));
% pca.placeAtPeak_allTrials_X = nan(size(pca.trials.active_X));
% pca.speed.AW_allTrials_X = nan(size(pca.trials.active_X));
% pca.corr_speed_timeAtPeak_allTrials_X = nan(size(pca.trials.active_X,1),2);
% for i=1:length(iscell)
%     if iscell(i)==1
%         this_data = traces_all.X(i,p.general.bins_analysisWindow,:);
%         if any(~isnan(this_data(:)))
%             [temp1,temp2] = nanmax(this_data,[],2);
%             pca.actAtPeak_allTrials_X(i,1:length(temp1)) = squeeze(temp1);
%             pca.binAtPeak_allTrials_X(i,1:length(temp2)) = squeeze(temp2)+p.general.bins_analysisWindow(1)-1;
%             pca.timeAtPeak_allTrials_X(i,1:length(temp2)) = p.general.t_binned(squeeze(temp2)+p.general.bins_analysisWindow(1)-1);
%             for k=1:length(pca.trials.X)
%                 pca.placeAtPeak_allTrials_X(i,k) = pca.distance(pca.trials.X(k),pca.binAtPeak_allTrials_X(i,k));
%             end
%             pca.speed.AW_allTrials_X(i,:) = pca.speed.AW(pca.trials.X);
%             [this_rho,this_pval] = corr(pca.speed.AW_allTrials_X(i,:)',pca.timeAtPeak_allTrials_X(i,:)','Type','Pearson','Rows','Complete');
%             pca.corr_speed_timeAtPeak_allTrials_X(i,:) = [this_rho,this_pval];
%         end
%     end
% end
% 
% % active A trials
% pca.actAtPeak_activeTrials_A = nan(size(pca.trials.active_A));
% pca.binAtPeak_activeTrials_A = nan(size(pca.trials.active_A));
% pca.timeAtPeak_activeTrials_A = nan(size(pca.trials.active_A));
% pca.placeAtPeak_activeTrials_A = nan(size(pca.trials.active_A));
% pca.speed.AW_activeTrials_A = nan(size(pca.trials.active_A));
% pca.corr_speed_timeAtPeak_activeTrials_A = nan(size(pca.trials.active_A,1),2);
% for i=1:length(iscell)
%     if iscell(i)==1
%         this_data = traces_all.A(i,p.general.bins_analysisWindow,find(pca.trials.active_A(i,:)));
%         if any(~isnan(this_data(:)))
%             [temp1,temp2] = nanmax(this_data,[],2);
%             pca.actAtPeak_activeTrials_A(i,1:length(temp1)) = squeeze(temp1);
%             pca.binAtPeak_activeTrials_A(i,1:length(temp2)) = squeeze(temp2)+p.general.bins_analysisWindow(1)-1;
%             pca.timeAtPeak_activeTrials_A(i,1:length(temp2)) = p.general.t_binned(squeeze(temp2)+p.general.bins_analysisWindow(1)-1);
%             temp1 = pca.trials.active_A(i,:); temp2 = find(temp1);
%             for k=1:length(temp2)
%                 pca.placeAtPeak_activeTrials_A(i,k) = pca.distance(pca.trials.A(temp2(k)),pca.binAtPeak_activeTrials_A(i,k));
%             end
%             pca.speed.AW_activeTrials_A(i,1:length(temp2)) = pca.speed.AW(pca.trials.A(temp2));
%             [this_rho,this_pval] = corr(pca.speed.AW_activeTrials_A(i,:)',pca.timeAtPeak_activeTrials_A(i,:)','Type','Pearson','Rows','Complete');
%             pca.corr_speed_timeAtPeak_activeTrials_A(i,:) = [this_rho,this_pval];
%         end
%     end
% end
% 
% % active X trials
% pca.actAtPeak_activeTrials_X = nan(size(pca.trials.active_X));
% pca.binAtPeak_activeTrials_X = nan(size(pca.trials.active_X));
% pca.timeAtPeak_activeTrials_X = nan(size(pca.trials.active_X));
% pca.placeAtPeak_activeTrials_X = nan(size(pca.trials.active_X));
% pca.speed.AW_activeTrials_X = nan(size(pca.trials.active_X));
% pca.corr_speed_timeAtPeak_activeTrials_X = nan(size(pca.trials.active_X,1),2);
% for i=1:length(iscell)
%     if iscell(i)==1
%         this_data = traces_all.X(i,p.general.bins_analysisWindow,find(pca.trials.active_X(i,:)));
%         if any(~isnan(this_data(:)))
%             [temp1,temp2] = nanmax(this_data,[],2);
%             pca.actAtPeak_activeTrials_X(i,1:length(temp1)) = squeeze(temp1);
%             pca.binAtPeak_activeTrials_X(i,1:length(temp2)) = squeeze(temp2)+p.general.bins_analysisWindow(1)-1;
%             pca.timeAtPeak_activeTrials_X(i,1:length(temp2)) = p.general.t_binned(squeeze(temp2)+p.general.bins_analysisWindow(1)-1);
%             temp1 = pca.trials.active_X(i,:); temp2 = find(temp1);
%             for k=1:length(temp2)
%                 pca.placeAtPeak_activeTrials_X(i,k) = pca.distance(pca.trials.X(temp2(k)),pca.binAtPeak_activeTrials_X(i,k));
%             end
%             pca.speed.AW_activeTrials_X(i,1:length(temp2)) = pca.speed.AW(pca.trials.X(temp2));
%             [this_rho,this_pval] = corr(pca.speed.AW_activeTrials_X(i,:)',pca.timeAtPeak_activeTrials_X(i,:)','Type','Pearson','Rows','Complete');
%             pca.corr_speed_timeAtPeak_activeTrials_X(i,:) = [this_rho,this_pval];
%         end
%     end
% end


%% Save data

pca_all = pca;
pca_all = orderfields(pca_all);
save([path.filepart_out,'pca_all.mat'],'pca_all','-v7.3');
disp(['--- Saved pca_all file as ',[path.filepart_out,'pca_all.mat'],'.'])


%%

% these_idcs = [90,260,261,264,266,196,41,254,289,274,267,216,204]; %find(passed_all.AW.Xonly==1);
% for i=1:length(these_idcs)
%     
% 	idx = these_idcs(i); %431
% 
%     figure
% 
%     subplot(2,3,1); hold on;
%     histogram(pca.timeAtPeak_activeTrials_X(idx,:))
%     xlabel('Time (s)')
% 
%     subplot(2,3,2); hold on;
%     histogram(pca.placeAtPeak_activeTrials_X(idx,:))
%     xlabel('Position (cm)')
% 
%     subplot(2,3,3); hold on;
%     scatter(pca.timeAtPeak_activeTrials_X(idx,:),pca.placeAtPeak_activeTrials_X(idx,:))
%     fitLine(pca.timeAtPeak_activeTrials_X(idx,:)',pca.placeAtPeak_activeTrials_X(idx,:)');
%     xlabel('Time (s)')
%     ylabel('Position (cm)')
% 
%     subplot(2,3,4); hold on;
%     scatter(pca.timeAtPeak_activeTrials_X(idx,:),pca.actAtPeak_activeTrials_X(idx,:))
%     fitLine(pca.timeAtPeak_activeTrials_X(idx,:)',pca.actAtPeak_activeTrials_X(idx,:)');
%     xlabel('Time (s)')
%     ylabel('Activity (z)')
% 
%     subplot(2,3,5); hold on;
%     scatter(pca.placeAtPeak_activeTrials_X(idx,:),pca.actAtPeak_activeTrials_X(idx,:))
%     fitLine(pca.placeAtPeak_activeTrials_X(idx,:)',pca.actAtPeak_activeTrials_X(idx,:)');
%     xlabel('Position (cm)')
%     ylabel('Activity (z)')
% 
%     subplot(2,3,6); hold on;
%     scatter(pca.speed.AW_activeTrials_X(idx,:),pca.timeAtPeak_activeTrials_X(idx,:))
%     fitLine(pca.speed.AW_activeTrials_X(idx,:)',pca.timeAtPeak_activeTrials_X(idx,:)');
%     xlabel('Speed (cm/s)')
%     ylabel('Peak time (s)')
% 
%     suptitle([num2str(idx),', active X trials'])
% 
% end


%% Mutual information

% % identify relevant frames per trial
% mia.frames_AW = 91:249;
% mia.numFrames_AW = length(mia.frames_AW);
% mia.t_AW = p.general.t_unbinned(91:249);
% 
% % identify relevant frames overall
% mia.frames_allInAW = pca_ub.trial_frames(mia.frames_AW,:);
% mia.frames_allInAW_flat = mia.frames_allInAW(:);
% 
% % get time
% mia.t_flat = repmat(mia.t_AW',size(mia.frames_allInAW,2),1)';
% mia.t_flat_z = zscore(mia.t_flat);
% 
% % get space
% space = nan(size(sync_beh.position));
% position_at_trial_start = NaN;
% for i=1:length(space)
%     if ismember(i,sync_beh.sync) %sync_beh.sync(i)==1
%         position_at_trial_start = sync_beh.position(i);
%     end
%     space(i) = sync_beh.position(i)-position_at_trial_start;
% end
% mia.x_flat = space(mia.frames_allInAW_flat);
% mia.x_flat_z = zscore(mia.x_flat);
% 
% % get activity
% mia.act = act(:,mia.frames_allInAW_flat);
% mia.act_z = zscore(mia.act,[],2);
% 
% % calculate mutual information
% mia.mi_t_act = nan(size(mia.act,1),1);
% mia.mi_x_act = nan(size(mia.act,1),1);
% for i=1:size(mia.act,1)
%     mia.mi_t_act(i) = mi_cont_cont(mia.t_flat_z,mia.act_z(i,:));
%     mia.mi_x_act(i) = mi_cont_cont(mia.x_flat_z,mia.act_z(i,:));
% end
% 
% % calculate correlation
% mia.corr_t_act = nan(size(mia.act,1),1);
% mia.corr_x_act = nan(size(mia.act,1),1);
% for i=1:size(mia.act,1)
%     mia.corr_t_act(i) = corr(mia.t_flat',mia.act(i,:)','Type','Pearson','Rows','Complete');
%     mia.corr_x_act(i) = corr(mia.x_flat',mia.act(i,:)','Type','Pearson','Rows','Complete');
% end

