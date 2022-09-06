function F = singleCellFigure(info,p,tng,idx,paq_beh)

%  calculate mean activity in firing field during active trials
% figure; tngtter(tng.firingField.A_AW.peakLocation_s(find(tng.passed.AW.A==1)),tng.firingField.A_AW.meanAmplitude_blSub_A(find(tng.passed.AW.A==1)))

%temp = find(tng.passed.AW.A==1);
%idx = 38  %temp(5)  % 23,  103          %67 %757 % 943 %temp(7) % e.g. 7
%idx = idcs_sorted_Aseq_dFF(1);
% tng.firingField.A_AW.diffBl_q(idx)
% tng.firingField.A_AW.diffBl_size(idx)

compact_version = 1;

if compact_version
    F = default_figure([20,0.5,5,9.8]);
    nrows = 8;
    ncols = 1;
else
    F = default_figure([20,0.5,20,9.9]);
    nrows = 8;
    ncols = 7;
end

ops.normMode_trialHeatmap = 'one-sided'; % one-sided, two-sided
ops.normAvg_heatmap = 'avg_minmax';
ops.normTrials_heatmap = 'trials_5_95';
ops.norm_traces = 'avg_trialTypes_minmax';


if ~isempty(find(idx==find(tng.prop.iscell)))
    this_iscell = find(idx==find(tng.prop.iscell));
else
    this_iscell = NaN;
end
if ~isempty(find(idx==find(tng.passed.AW.A==1)))
    this_seqA = find(idx==find(tng.passed.AW.A==1));
else
    this_seqA = NaN;
end
if ~isempty(find(idx==find(tng.passed.AW.X==1)))
    this_seqX = find(idx==find(tng.passed.AW.X==1));
else
    this_seqX = NaN;
end


%% Extract data

% originally used data and scaling for sequence cell analysis
origAvgData_A = permute(squeeze(tng.avgTraces.A(idx,:)),[2,1]);
origAvgData_X = permute(squeeze(tng.avgTraces.X(idx,:)),[2,1]);
origAvgData_B = permute(squeeze(tng.avgTraces.B(idx,:)),[2,1]);
origAvgData_Y = permute(squeeze(tng.avgTraces.Y(idx,:)),[2,1]);
temp = [origAvgData_A;origAvgData_X;origAvgData_B;origAvgData_Y];
origAvgData_lower = nanmin(temp(:));
origAvgData_upper = nanmax(temp(:));

% extract data
data_AB = permute(squeeze(tng.traces.A(idx,:,tng.trials.stimuli_inA.AB)),[2,1]);
data_AY = permute(squeeze(tng.traces.A(idx,:,tng.trials.stimuli_inA.AY)),[2,1]);
data_XY = permute(squeeze(tng.traces.X(idx,:,tng.trials.stimuli_inX.XY)),[2,1]);
data_XB = permute(squeeze(tng.traces.X(idx,:,tng.trials.stimuli_inX.XB)),[2,1]);
data_AB_H = permute(squeeze(tng.traces.A(idx,:,tng.trials.outcome_inA.A_H)),[2,1]);
data_AY_CR = permute(squeeze(tng.traces.A(idx,:,tng.trials.outcome_inA.A_CR)),[2,1]);
data_XY_H = permute(squeeze(tng.traces.X(idx,:,tng.trials.outcome_inX.X_H)),[2,1]);
data_XB_CR = permute(squeeze(tng.traces.X(idx,:,tng.trials.outcome_inX.X_CR)),[2,1]);
data_AB_M = permute(squeeze(tng.traces.A(idx,:,tng.trials.outcome_inA.A_M)),[2,1]);
data_AY_FA = permute(squeeze(tng.traces.A(idx,:,tng.trials.outcome_inA.A_FA)),[2,1]);
data_XY_M = permute(squeeze(tng.traces.X(idx,:,tng.trials.outcome_inX.X_M)),[2,1]);
data_XB_FA = permute(squeeze(tng.traces.X(idx,:,tng.trials.outcome_inX.X_FA)),[2,1]);

if size(data_AB_M,2)==1
    data_AB_M = data_AB_M';
end
if size(data_XY_M,2)==1
    data_XY_M = data_XY_M';
end
if size(data_AY_FA,2)==1
    data_AY_FA = data_AY_FA';
end
if size(data_XB_FA,2)==1
    data_XB_FA = data_XB_FA';
end
if size(data_AY_CR,2)==1
    data_AY_CR = data_AY_CR';
end
if size(data_XB_CR,2)==1
    data_XB_CR = data_XB_CR';
end

% average data
avgData_AB = nanmean(data_AB,1);
avgData_AY = nanmean(data_AY,1);
avgData_XY = nanmean(data_XY,1);
avgData_XB = nanmean(data_XB,1);
avgData_AB_H = nanmean(data_AB_H,1);
avgData_AY_CR = nanmean(data_AY_CR,1);
avgData_XY_H = nanmean(data_XY_H,1);
avgData_XB_CR = nanmean(data_XB_CR,1);
avgData_AB_M = nanmean(data_AB_M,1);
avgData_AY_FA = nanmean(data_AY_FA,1);
avgData_XY_M = nanmean(data_XY_M,1);
avgData_XB_FA = nanmean(data_XB_FA,1);

% normalise average data (only used for average heatmap)
if strcmp(ops.normAvg_heatmap,'avg_minmax')
    temp = [avgData_AB;avgData_AY;avgData_XY;avgData_XB];
    this_lower = nanmin(temp(:));
    this_upper = nanmax(temp(:));
    normAvgData_heatmap1s_AB = (avgData_AB-this_lower) ./ (this_upper-this_lower);
    normAvgData_heatmap1s_AY = (avgData_AY-this_lower) ./ (this_upper-this_lower);
    normAvgData_heatmap1s_XY = (avgData_XY-this_lower) ./ (this_upper-this_lower);
    normAvgData_heatmap1s_XB = (avgData_XB-this_lower) ./ (this_upper-this_lower);
end
if strcmp(ops.normAvg_heatmap,'avg_minmax')
    temp = [avgData_AB;avgData_AY;avgData_XY;avgData_XB];
    temp2 = temp(:,(1:p.general.bins_pre));
    this_zero = nanmedian(temp2(:));
    if abs(nanmin(temp(:))-this_zero) > abs(nanmax(temp(:))-this_zero)
        this_lower = nanmin(temp(:));
        this_upper = abs(nanmin(temp(:)))+2*this_zero;
    else
        this_lower = -nanmax(temp(:))+2*this_zero;
        this_upper = nanmax(temp(:));            
    end
    normAvgData_heatmap2s_AB = rescale([avgData_AB,this_lower,this_zero,this_upper],-1,1); normAvgData_heatmap2s_AB = normAvgData_heatmap2s_AB(1:end-3);
    normAvgData_heatmap2s_AY = rescale([avgData_AY,this_lower,this_zero,this_upper],-1,1); normAvgData_heatmap2s_AY = normAvgData_heatmap2s_AY(1:end-3);
    normAvgData_heatmap2s_XY = rescale([avgData_XY,this_lower,this_zero,this_upper],-1,1); normAvgData_heatmap2s_XY = normAvgData_heatmap2s_XY(1:end-3);
    normAvgData_heatmap2s_XB = rescale([avgData_XB,this_lower,this_zero,this_upper],-1,1); normAvgData_heatmap2s_XB = normAvgData_heatmap2s_XB(1:end-3);        
end

% normalise data (only used for trial-wise heatmap)
if strcmp(ops.normMode_trialHeatmap,'one-sided')
    if strcmp(ops.normTrials_heatmap,'avg_minmax')
        temp = [avgData_AB;avgData_AY;avgData_XY;avgData_XB];
        this_lower = nanmin(temp(:));
        this_upper = nanmax(temp(:));
        normData_heatmap_AB = (data_AB-this_lower) ./ (this_upper-this_lower);
        normData_heatmap_AY = (data_AY-this_lower) ./ (this_upper-this_lower);
        normData_heatmap_XY = (data_XY-this_lower) ./ (this_upper-this_lower);
        normData_heatmap_XB = (data_XB-this_lower) ./ (this_upper-this_lower);
    elseif strcmp(ops.normTrials_heatmap,'trials_minmax')
        temp = [data_AB;data_AY;data_XY;data_XB];
        this_lower = nanmin(temp(:));
        this_upper = nanmax(temp(:));
        normData_heatmap_AB = (data_AB-this_lower) ./ (this_upper-this_lower);
        normData_heatmap_AY = (data_AY-this_lower) ./ (this_upper-this_lower);
        normData_heatmap_XY = (data_XY-this_lower) ./ (this_upper-this_lower);
        normData_heatmap_XB = (data_XB-this_lower) ./ (this_upper-this_lower);
    elseif strcmp(ops.normTrials_heatmap,'trials_1_99')
        temp = [data_AB;data_AY;data_XY;data_XB];
        temp2 = temp(:,(1:p.general.bins_pre));
        this_zero = nanmedian(temp2(:));
        if abs(prctile(rmmissing(temp(:)),1)-this_zero) > abs(prctile(rmmissing(temp(:)),99)-this_zero)
            this_lower = prctile(rmmissing(temp(:)),1);
            this_upper = abs(prctile(rmmissing(temp(:)),1))+2*this_zero;
        else
            this_lower = -prctile(rmmissing(temp(:)),99)+2*this_zero;
            this_upper = prctile(rmmissing(temp(:)),99);            
        end
        normData_heatmap_AB = rescale([data_AB,ones(size(data_AB,1),1)*[this_lower,this_zero,this_upper]],-1,1,'InputMin',this_lower,'InputMax',this_upper); normData_heatmap_AB = normData_heatmap_AB(:,1:end-3);
        normData_heatmap_AY = rescale([data_AY,ones(size(data_AY,1),1)*[this_lower,this_zero,this_upper]],-1,1,'InputMin',this_lower,'InputMax',this_upper); normData_heatmap_AY = normData_heatmap_AY(:,1:end-3);
        normData_heatmap_XY = rescale([data_XY,ones(size(data_XY,1),1)*[this_lower,this_zero,this_upper]],-1,1,'InputMin',this_lower,'InputMax',this_upper); normData_heatmap_XY = normData_heatmap_XY(:,1:end-3);
        normData_heatmap_XB = rescale([data_XB,ones(size(data_XB,1),1)*[this_lower,this_zero,this_upper]],-1,1,'InputMin',this_lower,'InputMax',this_upper); normData_heatmap_XB = normData_heatmap_XB(:,1:end-3);         
    elseif strcmp(ops.normTrials_heatmap,'trials_5_95')
        temp = [data_AB;data_AY;data_XY;data_XB];
        temp2 = temp(:,(1:p.general.bins_pre));
        this_zero = nanmedian(temp2(:));
        if abs(prctile(rmmissing(temp(:)),5)-this_zero) > abs(prctile(rmmissing(temp(:)),95)-this_zero)
            this_lower = prctile(rmmissing(temp(:)),5);
            this_upper = abs(prctile(rmmissing(temp(:)),5))+2*this_zero;
        else
            this_lower = -prctile(rmmissing(temp(:)),95)+2*this_zero;
            this_upper = prctile(rmmissing(temp(:)),95);            
        end
        normData_heatmap_AB = rescale([data_AB,ones(size(data_AB,1),1)*[this_lower,this_zero,this_upper]],-1,1,'InputMin',this_lower,'InputMax',this_upper); normData_heatmap_AB = normData_heatmap_AB(:,1:end-3);
        normData_heatmap_AY = rescale([data_AY,ones(size(data_AY,1),1)*[this_lower,this_zero,this_upper]],-1,1,'InputMin',this_lower,'InputMax',this_upper); normData_heatmap_AY = normData_heatmap_AY(:,1:end-3);
        normData_heatmap_XY = rescale([data_XY,ones(size(data_XY,1),1)*[this_lower,this_zero,this_upper]],-1,1,'InputMin',this_lower,'InputMax',this_upper); normData_heatmap_XY = normData_heatmap_XY(:,1:end-3);
        normData_heatmap_XB = rescale([data_XB,ones(size(data_XB,1),1)*[this_lower,this_zero,this_upper]],-1,1,'InputMin',this_lower,'InputMax',this_upper); normData_heatmap_XB = normData_heatmap_XB(:,1:end-3);         
    end
elseif strcmp(ops.normMode_trialHeatmap,'two-sided')
    if strcmp(ops.normTrials_heatmap,'avg_minmax')
        disp('Normalisation not yet implemented')
    elseif strcmp(ops.normTrials_heatmap,'trials_minmax')
        disp('Normalisation not yet implemented')
    end
end

% normalise data (used for all traces, except for originals)
if strcmp(ops.norm_traces,'avg_trialTypes_minmax')
    temp = [avgData_AB;avgData_AY;avgData_XY;avgData_XB];
    this_lower = nanmin(temp(:));
    this_upper = nanmax(temp(:));
    normData_traces_AB = (data_AB-this_lower) ./ (this_upper-this_lower);
    normData_traces_AY = (data_AY-this_lower) ./ (this_upper-this_lower);
    normData_traces_XY = (data_XY-this_lower) ./ (this_upper-this_lower);
    normData_traces_XB = (data_XB-this_lower) ./ (this_upper-this_lower);
    normData_traces_AB_H = (data_AB_H-this_lower) ./ (this_upper-this_lower);
    normData_traces_AY_CR = (data_AY_CR-this_lower) ./ (this_upper-this_lower);
    normData_traces_XY_H = (data_XY_H-this_lower) ./ (this_upper-this_lower);
    normData_traces_XB_CR = (data_XB_CR-this_lower) ./ (this_upper-this_lower);
    normData_traces_AB_M = (data_AB_M-this_lower) ./ (this_upper-this_lower);
    normData_traces_AY_FA = (data_AY_FA-this_lower) ./ (this_upper-this_lower);
    normData_traces_XY_M = (data_XY_M-this_lower) ./ (this_upper-this_lower);
    normData_traces_XB_FA = (data_XB_FA-this_lower) ./ (this_upper-this_lower);
    temp = [normData_traces_AB;normData_traces_AY;normData_traces_XY;normData_traces_XB];
    normData_lower = nanmin(temp(:));
    normData_upper = nanmax(temp(:));
end
avgNormData_traces_AB = nanmean(normData_traces_AB,1);
avgNormData_traces_AY = nanmean(normData_traces_AY,1);
avgNormData_traces_XY = nanmean(normData_traces_XY,1);
avgNormData_traces_XB = nanmean(normData_traces_XB,1);
avgNormData_traces_AB_H = nanmean(normData_traces_AB_H,1);
avgNormData_traces_AY_CR = nanmean(normData_traces_AY_CR,1);
avgNormData_traces_XY_H = nanmean(normData_traces_XY_H,1);
avgNormData_traces_XB_CR = nanmean(normData_traces_XB_CR,1);
avgNormData_traces_AB_M = nanmean(normData_traces_AB_M,1);
avgNormData_traces_AY_FA = nanmean(normData_traces_AY_FA,1);
avgNormData_traces_XY_M = nanmean(normData_traces_XY_M,1);
avgNormData_traces_XB_FA = nanmean(normData_traces_XB_FA,1);
temp = [avgNormData_traces_AB;avgNormData_traces_AY;avgNormData_traces_XY;avgNormData_traces_XB;...
    avgNormData_traces_AB_H;avgNormData_traces_AY_CR;avgNormData_traces_XY_H;avgNormData_traces_XB_CR;...
    avgNormData_traces_AB_M;avgNormData_traces_AY_FA;avgNormData_traces_XY_M;avgNormData_traces_XB_FA];
avgNormData_inclSplit_lower = nanmin(temp(:));
avgNormData_inclSplit_upper = nanmax(temp(:));
temp = [avgNormData_traces_AB;avgNormData_traces_AY;avgNormData_traces_XY;avgNormData_traces_XB];
avgNormData_exclSplit_lower = nanmin(temp(:));
avgNormData_exclSplit_upper = nanmax(temp(:));


%% Heatmaps

if compact_version
    subplot(nrows,ncols,sort([1:6]))
else
    subplot(nrows,ncols,sort([([1:6]-1)*ncols+1,([1:6]-1)*ncols+2]))
end
temp = movmean([0;cumsum([size(normData_heatmap_AB,1);size(normData_heatmap_AY,1);size(normData_heatmap_XY,1);size(normData_heatmap_XB,1)])],2);
plt = struct(); plt.illustrator = true; plt.yticks = temp(2:end); plt.yticklabels = {'AB','AY','XY','XB'};
if strcmp(ops.normMode_trialHeatmap,'one-sided')
    plt.colormap = parula; plt.clim = [0,1]; % comment out clim for smoother look
elseif strcmp(ops.normMode_trialHeatmap,'two-sided')
    plt.colormap = redblue; plt.clim = [-1,1];
end
heatMap_task([normData_heatmap_AB;normData_heatmap_AY;normData_heatmap_XY;normData_heatmap_XB],cumsum([size(normData_heatmap_AB,1);size(normData_heatmap_AY,1);size(normData_heatmap_XY,1)]),tng,p,info,plt);

if ~compact_version
    subplot(nrows,ncols,(7-1)*ncols+[1:2])
    plt = struct(); plt.illustrator = false; plt.yticks = 1:4; plt.yticklabels = {'AB','AY','XY','XB'}; plt.clim = [0,1]; plt.colormap = parula;
    heatMap_task([normAvgData_heatmap1s_AB;normAvgData_heatmap1s_AY;normAvgData_heatmap1s_XY;normAvgData_heatmap1s_XB],[],tng,p,info,plt);

    subplot(nrows,ncols,(8-1)*ncols+[1:2])
    plt = struct(); plt.illustrator = false; plt.yticks = 1:4; plt.yticklabels = {'AB','AY','XY','XB'}; plt.clim = [-1,1]; plt.colormap = redblue;
    heatMap_task([normAvgData_heatmap2s_AB;normAvgData_heatmap2s_AY;normAvgData_heatmap2s_XY;normAvgData_heatmap2s_XB],[],tng,p,info,plt);
end

%% Licking and running

% if exist('paq_beh','var')
%     
%     subplot(nrows,ncols,(6-1)*ncols+[3:4])
%     if nanmax(tng.prop.trial_frames(:))>length(paq_beh.speed)
%         temp = paq_beh.speed;
%         this_speed = NaN(1,nanmax(tng.prop.trial_frames(:)));
%         this_speed(1,1:size(temp,2)) = temp;
%     else
%         this_speed = paq_beh.speed;
%     end
%     if true % smoothing
%         this_speed = smoothdata(this_speed,2,'gaussian',5*5);
%     end
%     this_speed = reshape(this_speed(tng.prop.trial_frames),tng.prop.numTrials,tng.prop.numFrames_raw);
%     
%     % binning
%     temp = movmean(this_speed,p.general.binSize,2,'omitnan');
%     this_speed_binned = temp(:,floor(p.general.binSize/2)+1:p.general.binSize:end,:);
% 
%     hold on
%     yline(0,'Color',p.col.black,'LineStyle',':','LineWidth',1);
%     these_trials = intersect(trials.trials_X,trials.trials_B);
%     temp=shadedErrorBar(1:size(this_speed_binned(these_trials,:),2),nanmean(this_speed_binned(these_trials,:),1),nansem(this_speed_binned(these_trials,:),1),'lineProps',p.col.XB); temp.mainLine.LineWidth = 2;   
%     hold on
%     these_trials = intersect(trials.trials_A,trials.trials_Y);
%     temp=shadedErrorBar(1:size(this_speed_binned(these_trials,:),2),nanmean(this_speed_binned(these_trials,:),1),nansem(this_speed_binned(these_trials,:),1),'lineProps',p.col.AY); temp.mainLine.LineWidth = 2;   
%     these_trials = intersect(trials.trials_X,trials.trials_Y);
%     temp=shadedErrorBar(1:size(this_speed_binned(these_trials,:),2),nanmean(this_speed_binned(these_trials,:),1),nansem(this_speed_binned(these_trials,:),1),'lineProps',p.col.XY); temp.mainLine.LineWidth = 2;   
%     these_trials = intersect(trials.trials_A,trials.trials_B);
%     temp=shadedErrorBar(1:size(this_speed_binned(these_trials,:),2),nanmean(this_speed_binned(these_trials,:),1),nansem(this_speed_binned(these_trials,:),1),'lineProps',p.col.AB); temp.mainLine.LineWidth = 2;   
%     plt = struct(); plt.ylabel = 'Speed (cm/s)'; plt.ylim = [-10,160]; traces_task(tng,p,info,plt); 
%     hold off
% 
% end

%% Traces

% subplot(nrows,ncols,(3-1)*ncols+3)
% hold on
% plot(normData_traces_AB','Color',p.col.gray,'LineWidth',0.5)
% plot(avgNormData_traces_AB,'Color',p.col.AB,'LineWidth',2)
% plt = struct(); plt.ylim = [normData_lower,normData_upper]; traces_task(tng,p,info,plt); 
% 
% subplot(nrows,ncols,(4-1)*ncols+3)
% hold on
% plot(normData_traces_AY','Color',p.col.gray,'LineWidth',0.5)
% plot(avgNormData_traces_AY,'Color',p.col.AY,'LineWidth',2)
% plt = struct(); plt.ylim = [normData_lower,normData_upper]; traces_task(tng,p,info,plt); 
% 
% subplot(nrows,ncols,(4-1)*ncols+4)
% hold on
% plot(normData_traces_XY','Color',p.col.gray,'LineWidth',0.5)
% plot(avgNormData_traces_XY,'Color',p.col.XY,'LineWidth',2)
% plt = struct(); plt.ylim = [normData_lower,normData_upper]; traces_task(tng,p,info,plt); 
% 
% subplot(nrows,ncols,(3-1)*ncols+4)
% hold on
% plot(normData_traces_XB','Color',p.col.gray,'LineWidth',0.5)
% plot(avgNormData_traces_XB,'Color',p.col.XB,'LineWidth',2)
% plt = struct(); plt.ylim = [normData_lower,normData_upper]; traces_task(tng,p,info,plt); 


%%

if compact_version
    subplot(nrows,ncols,sort([7:8]))
else
	subplot(nrows,ncols,sort([([1:2]-1)*ncols+3,([1:2]-1)*ncols+4]))
end
hold on
temp=shadedErrorBar(1:size(normData_traces_XB,2),nanmean(normData_traces_XB,1),nansem(normData_traces_XB,1),'lineProps',p.col.XB); temp.mainLine.LineWidth = 2;   
temp=shadedErrorBar(1:size(normData_traces_AY,2),nanmean(normData_traces_AY,1),nansem(normData_traces_AY,1),'lineProps',p.col.AY); temp.mainLine.LineWidth = 2;   
temp=shadedErrorBar(1:size(normData_traces_XY,2),nanmean(normData_traces_XY,1),nansem(normData_traces_XY,1),'lineProps',p.col.XY); temp.mainLine.LineWidth = 2;   
temp=shadedErrorBar(1:size(normData_traces_AB,2),nanmean(normData_traces_AB,1),nansem(normData_traces_AB,1),'lineProps',p.col.AB); temp.mainLine.LineWidth = 2;   
plt = struct(); plt.ylim = [avgNormData_exclSplit_lower-0.2,avgNormData_exclSplit_upper+0.2]; traces_task(tng,p,info,plt); 
hold off

%%

if ~compact_version
    subplot(nrows,ncols,(3-1)*ncols+[3:4])
    hold on
    temp=shadedErrorBar(1:size(normData_traces_AY_FA,2),nanmean(normData_traces_AY_FA,1),nansem(normData_traces_AY_FA,1),'lineProps',p.col.AY); temp.mainLine.LineWidth = 2; temp.mainLine.LineStyle = ':';  
    temp=shadedErrorBar(1:size(normData_traces_AB_M,2),nanmean(normData_traces_AB_M,1),nansem(normData_traces_AB_M,1),'lineProps',p.col.AB); temp.mainLine.LineWidth = 2; temp.mainLine.LineStyle = ':';  
    temp=shadedErrorBar(1:size(normData_traces_AY_CR,2),nanmean(normData_traces_AY_CR,1),nansem(normData_traces_AY_CR,1),'lineProps',p.col.AY); temp.mainLine.LineWidth = 2;
    temp=shadedErrorBar(1:size(normData_traces_AB_H,2),nanmean(normData_traces_AB_H,1),nansem(normData_traces_AB_H,1),'lineProps',p.col.AB); temp.mainLine.LineWidth = 2; 
    plt = struct(); plt.ylim = [avgNormData_inclSplit_lower-0.2,avgNormData_inclSplit_upper+0.2]; traces_task(tng,p,info,plt); 
    hold off

    subplot(nrows,ncols,(4-1)*ncols+[3:4])
    hold on
    temp=shadedErrorBar(1:size(normData_traces_XB_FA,2),nanmean(normData_traces_XB_FA,1),nansem(normData_traces_XB_FA,1),'lineProps',p.col.XB); temp.mainLine.LineWidth = 2; temp.mainLine.LineStyle = ':';  
    temp=shadedErrorBar(1:size(normData_traces_XY_M,2),nanmean(normData_traces_XY_M,1),nansem(normData_traces_XY_M,1),'lineProps',p.col.XY); temp.mainLine.LineWidth = 2; temp.mainLine.LineStyle = ':';  
    temp=shadedErrorBar(1:size(normData_traces_XB_CR,2),nanmean(normData_traces_XB_CR,1),nansem(normData_traces_XB_CR,1),'lineProps',p.col.XB); temp.mainLine.LineWidth = 2;   
    temp=shadedErrorBar(1:size(normData_traces_XY_H,2),nanmean(normData_traces_XY_H,1),nansem(normData_traces_XY_H,1),'lineProps',p.col.XY); temp.mainLine.LineWidth = 2;   
    plt = struct(); plt.ylim = [avgNormData_inclSplit_lower-0.2,avgNormData_inclSplit_upper+0.2]; traces_task(tng,p,info,plt); 
    hold off
end

%% Odour B trials vs Odour Y trials traces and selectivity

if ~compact_version

    subplot(nrows,ncols,(7-1)*ncols+3)
    hold on
    % yline(tng.prop.baseline(idx),'Color',p.col.black,'LineStyle',':','LineWidth',1);
    % yline(tng.prop.threshold(idx),'Color',p.col.photostim,'LineStyle','-','LineWidth',1);
    plot(origAvgData_B,'Color',nanmean([p.col.AB;p.col.XB],1),'LineWidth',2)
    ylim([origAvgData_lower,origAvgData_upper])
    plt = struct(); traces_task(tng,p,info,plt); 
    title('Odour B trials')

    subplot(nrows,ncols,(7-1)*ncols+4)
    hold on
    % yline(tng.prop.baseline(idx),'Color',p.col.black,'LineStyle',':','LineWidth',1);
    % yline(tng.prop.threshold(idx),'Color',p.col.photostim,'LineStyle','-','LineWidth',1);
    plot(origAvgData_Y,'Color',nanmean([p.col.XY;p.col.AY],1),'LineWidth',2)
    ylim([origAvgData_lower,origAvgData_upper])
    plt = struct(); traces_task(tng,p,info,plt); 
    title('Odour Y trials')
end
% 
% % Selectivity B
% subplot(nrows,ncols,(8-1)*ncols+3)
% plt = struct(); plt.xlabel = 'Selectivity'; plt.xlim = [-1,1]; plt.ylabel = 'Cells'; plt.norm = 'count'; plt.edges = [-1:0.1:1];
% plt.dataPoint = tng.odour2.selectivity_AB_XB(idx);
% myHistogram(tng.odour2.selectivity_AB_XB(find(tng.prop.iscell)),p,plt); temp1 = get(gca,'ylim');
% 
% % Selectivity Y
% subplot(nrows,ncols,(8-1)*ncols+4)
% plt = struct(); plt.xlabel = 'Selectivity'; plt.xlim = [-1,1]; plt.ylabel = 'Cells'; plt.norm = 'count'; plt.edges = [-1:0.1:1];
% plt.dataPoint = tng.odour2.selectivity_XY_AY(idx);
% myHistogram(tng.odour2.selectivity_XY_AY(find(tng.prop.iscell)),p,plt); temp2 = get(gca,'ylim');
% 
% subplot(nrows,ncols,(8-1)*ncols+3)
% ylim([0,nanmax([temp1(2),temp2(2)])])
% subplot(nrows,ncols,(8-1)*ncols+4)
% ylim([0,nanmax([temp1(2),temp2(2)])])


%% Odour A trials vs Odour X trials traces

if ~compact_version

    subplot(nrows,ncols,(1-1)*ncols+5)
    hold on
    yline(tng.firingField.A_AW.baseline_meanOfTrialwise(idx),'Color',p.col.black,'LineStyle',':','LineWidth',1);
    yline(tng.firingField.A_AW.activeInFiringFieldThreshold(idx),'Color',p.col.photostim,'LineStyle',':','LineWidth',1);
    if ~isnan(this_seqA)
        %yline(tng.firingField.A_AW.peakAmplitude(idx),'Color',p.col.black,'LineStyle','-','LineWidth',1);
        line([tng.firingField.A_AW.peakLocation(idx),tng.firingField.A_AW.peakLocation(idx)],[tng.firingField.A_AW.baseline_meanOfTrialwise(idx),tng.firingField.A_AW.peakAmplitude(idx)],'Color',p.col.black,'LineStyle','-','LineWidth',2)
        temp = tng.firingField.A_AW.baseline_meanOfTrialwise + p.tng.firingFieldBoundaries*tng.firingField.A_AW.peakAmplitude_blSub;
        line([tng.firingField.A_AW.boundaries(idx,1),tng.firingField.A_AW.boundaries(idx,2)],[temp(idx),temp(idx)],'Color',p.col.black,'LineStyle','-','LineWidth',2)
    end
    plot(origAvgData_A,'Color',nanmean([p.col.AB;p.col.AY],1),'LineWidth',2)
    ylim([origAvgData_lower,nanmax([origAvgData_upper,tng.firingField.A_AW.activeInFiringFieldThreshold(idx),tng.firingField.X_AW.activeInFiringFieldThreshold(idx)])])
    plt = struct(); traces_task(tng,p,info,plt);
    title('Odour A trials')

    subplot(nrows,ncols,(1-1)*ncols+6)
    hold on
    yline(tng.firingField.X_AW.baseline_meanOfTrialwise(idx),'Color',p.col.black,'LineStyle',':','LineWidth',1);
    yline(tng.firingField.X_AW.activeInFiringFieldThreshold(idx),'Color',p.col.photostim,'LineStyle',':','LineWidth',1);
    if ~isnan(this_seqX)
        %yline(tng.firingField.X_AW.peakAmplitude(idx),'Color',p.col.black,'LineStyle','-','LineWidth',1);
        line([tng.firingField.X_AW.peakLocation(idx),tng.firingField.X_AW.peakLocation(idx)],[tng.firingField.X_AW.baseline(idx),tng.firingField.X_AW.peakAmplitude(idx)],'Color',p.col.black,'LineStyle','-','LineWidth',2)
        temp = tng.firingField.X_AW.baseline + p.tng.firingFieldBoundaries*tng.firingField.X_AW.peakAmplitude_blSub;
        line([tng.firingField.X_AW.boundaries(idx,1),tng.firingField.X_AW.boundaries(idx,2)],[temp(idx),temp(idx)],'Color',p.col.black,'LineStyle','-','LineWidth',2)
    end
    plot(origAvgData_X,'Color',nanmean([p.col.XY;p.col.XB],1),'LineWidth',2)
    ylim([origAvgData_lower,nanmax([origAvgData_upper,tng.firingField.A_AW.activeInFiringFieldThreshold(idx),tng.firingField.X_AW.activeInFiringFieldThreshold(idx)])])
    plt = struct(); traces_task(tng,p,info,plt); 
    title('Odour X trials')
end

%% Sequence cell criteria and properties

if ~compact_version

    % Trial-wise amplitude A
    if tng.shuffling.A_AW.significant(idx)==1 && (~isnan(nanmean(tng.firingField.A_AW.diffBl_diff(idx,:))))
        subplot(nrows,ncols,(3-1)*ncols+5)
        temp = tng.firingField.A_AW.diffBl_diff(idx,:);
        plt = struct(); plt.xlabel = 'Trial-wise amplitude'; plt.xlim = [nanmin(temp),nanmax(temp)]; plt.ylabel = 'Trials'; plt.numBins = 20; plt.criterion = 0;
        plt.dataPoint = nanmean(tng.firingField.A_AW.diffBl_diff(idx,:));
        myHistogram(tng.firingField.A_AW.diffBl_diff(idx,:),p,plt); temp1 = get(gca,'ylim');
        if tng.firingField.A_AW.diffBl(idx)
            title('*')
        end
    end

    % Trial-wise amplitude X
    if tng.shuffling.X_AW.significant(idx)==1 && (~isnan(nanmean(tng.firingField.X_AW.diffBl_diff(idx,:))))
        subplot(nrows,ncols,(3-1)*ncols+6)
        temp = tng.firingField.X_AW.diffBl_diff(idx,:);
        plt = struct(); plt.xlabel = 'Trial-wise amplitude'; plt.xlim = [nanmin(temp),nanmax(temp)]; plt.ylabel = 'Trials'; plt.numBins = 20; plt.criterion = 0;
        plt.dataPoint = nanmean(tng.firingField.X_AW.diffBl_diff(idx,:));
        myHistogram(tng.firingField.X_AW.diffBl_diff(idx,:),p,plt); temp2 = get(gca,'ylim');
        if tng.firingField.X_AW.diffBl(idx)
            title('*')
        end
    end

    if (tng.shuffling.A_AW.significant(idx)==1 && tng.shuffling.X_AW.significant(idx)==1) && ((~isnan(nanmean(tng.firingField.A_AW.diffBl_diff(idx,:)))) && (~isnan(nanmean(tng.firingField.X_AW.diffBl_diff(idx,:)))))
        subplot(nrows,ncols,(3-1)*ncols+5)
        ylim([0,nanmax([temp1(2),temp2(2)])])
        subplot(nrows,ncols,(3-1)*ncols+6)
        ylim([0,nanmax([temp1(2),temp2(2)])])
    end

    % Shuffling A
    subplot(nrows,ncols,(2-1)*ncols+5)
    temp = [tng.shuffling.A_AW.maxRates_true(idx),tng.shuffling.A_AW.maxRates_shuffled(idx,:)];
    plt = struct(); plt.xlabel = 'Shuffling'; plt.xlim = [nanmin(temp),nanmax(temp)]; plt.ylabel = 'Shuffles'; plt.numBins = 20; plt.criterion = prctile(tng.shuffling.A_AW.maxRates_shuffled(idx,:),tng.p.tng.shufflingThreshold,2);
    plt.dataPoint = tng.shuffling.A_AW.maxRates_true(idx);
    myHistogram(tng.shuffling.A_AW.maxRates_shuffled(idx,:),p,plt); temp1 = get(gca,'ylim');

    % Shuffling X
    subplot(nrows,ncols,(2-1)*ncols+6)
    temp = [tng.shuffling.X_AW.maxRates_true(idx),tng.shuffling.X_AW.maxRates_shuffled(idx,:)];
    plt = struct(); plt.xlabel = 'Shuffling'; plt.xlim = [nanmin(temp),nanmax(temp)]; plt.ylabel = 'Shuffles'; plt.numBins = 20; plt.criterion = prctile(tng.shuffling.X_AW.maxRates_shuffled(idx,:),tng.p.tng.shufflingThreshold,2);
    plt.dataPoint = tng.shuffling.X_AW.maxRates_true(idx);
    myHistogram(tng.shuffling.X_AW.maxRates_shuffled(idx,:),p,plt); temp2 = get(gca,'ylim');

    subplot(nrows,ncols,(2-1)*ncols+5)
    ylim([0,nanmax([temp1(2),temp2(2)])])
    subplot(nrows,ncols,(2-1)*ncols+6)
    ylim([0,nanmax([temp1(2),temp2(2)])])

    % Reliability A
    subplot(nrows,ncols,(4-1)*ncols+5)
    plt = struct(); plt.xlabel = 'Reliability'; plt.xlim = [0,1]; plt.ylabel = 'A cells'; plt.norm = 'count'; plt.edges = [0:0.05:1]; plt.criterion = tng.p.tng.reliabilityThreshold;
    plt.dataPoint = tng.firingField.A_AW.activationProbability_ipsi(idx);
    hold on
    myHistogram(tng.firingField.A_AW.activationProbability_ipsi(find(tng.shuffling.A_AW.significant==1)),p,plt); temp1 = get(gca,'ylim');
    plt.color = p.col.darkGray;
    myHistogram(tng.firingField.A_AW.activationProbability_ipsi(find(tng.passed.AW.A==1)),p,plt); temp1 = get(gca,'ylim');
    hold off

    % Reliability X
    subplot(nrows,ncols,(4-1)*ncols+6)
    plt = struct(); plt.xlabel = 'Reliability'; plt.xlim = [0,1]; plt.ylabel = 'X cells'; plt.norm = 'count'; plt.edges = [0:0.05:1]; plt.criterion = tng.p.tng.reliabilityThreshold;
    plt.dataPoint = tng.firingField.X_AW.activationProbability_ipsi(idx);
    hold on
    myHistogram(tng.firingField.X_AW.activationProbability_ipsi(find(tng.shuffling.X_AW.significant==1)),p,plt); temp2 = get(gca,'ylim');
    plt.color = p.col.darkGray;
    myHistogram(tng.firingField.X_AW.activationProbability_ipsi(find(tng.passed.AW.X==1)),p,plt); temp2 = get(gca,'ylim');
    hold off

    subplot(nrows,ncols,(4-1)*ncols+5)
    ylim([0,nanmax([temp1(2),temp2(2)])])
    subplot(nrows,ncols,(4-1)*ncols+6)
    ylim([0,nanmax([temp1(2),temp2(2)])])

    % Selectivity A
    subplot(nrows,ncols,(5-1)*ncols+5)
    plt = struct(); plt.xlabel = 'Selectivity'; plt.xlim = [-1,1]; plt.ylabel = 'A cells'; plt.norm = 'count'; plt.edges = [-1:0.1:1];
    plt.dataPoint = tng.firingField.A_AW.selectivity_ipsi(idx);
    hold on
    myHistogram(tng.firingField.A_AW.selectivity_ipsi(find(tng.shuffling.A_AW.significant==1)),p,plt); temp1 = get(gca,'ylim');
    plt.color = p.col.darkGray;
    myHistogram(tng.firingField.A_AW.selectivity_ipsi(find(tng.passed.AW.A==1)),p,plt); temp1 = get(gca,'ylim');
    hold off

    % Selectivity X
    subplot(nrows,ncols,(5-1)*ncols+6)
    plt = struct(); plt.xlabel = 'Selectivity'; plt.xlim = [-1,1]; plt.ylabel = 'X cells'; plt.norm = 'count'; plt.edges = [-1:0.1:1];
    plt.dataPoint = tng.firingField.X_AW.selectivity_ipsi(idx);
    hold on
    myHistogram(tng.firingField.X_AW.selectivity_ipsi(find(tng.shuffling.X_AW.significant==1)),p,plt); temp2 = get(gca,'ylim');
    plt.color = p.col.darkGray;
    myHistogram(tng.firingField.X_AW.selectivity_ipsi(find(tng.passed.AW.X==1)),p,plt); temp2 = get(gca,'ylim');
    hold off

    subplot(nrows,ncols,(5-1)*ncols+5)
    ylim([0,nanmax([temp1(2),temp2(2)])])
    subplot(nrows,ncols,(5-1)*ncols+6)
    ylim([0,nanmax([temp1(2),temp2(2)])])

    % Amplitude A
    subplot(nrows,ncols,(6-1)*ncols+5)
    temp = [tng.firingField.A_AW.meanAmplitude_blSub_active_ipsi(idx);tng.firingField.A_AW.meanAmplitude_blSub_active_ipsi(find(tng.passed.AW.A==1));tng.firingField.X_AW.meanAmplitude_blSub_active_ipsi(idx);tng.firingField.X_AW.meanAmplitude_blSub_active_ipsi(find(tng.passed.AW.X==1))];
    plt = struct(); plt.xlabel = 'Mean amplitude (active)';  plt.xlim = [0,nanmax(temp)]; plt.ylabel = 'A cells'; plt.norm = 'count'; plt.edges = [0:0.05:10];
    plt.dataPoint = tng.firingField.A_AW.meanAmplitude_blSub_active_ipsi(idx);
    hold on
    myHistogram(tng.firingField.A_AW.meanAmplitude_blSub_active_ipsi(find(tng.shuffling.A_AW.significant==1)),p,plt); temp1 = get(gca,'ylim');
    plt.color = p.col.darkGray;
    myHistogram(tng.firingField.A_AW.meanAmplitude_blSub_active_ipsi(find(tng.passed.AW.A==1)),p,plt); temp1 = get(gca,'ylim');
    hold off

    % Amplitude X
    subplot(nrows,ncols,(6-1)*ncols+6)
    temp = [tng.firingField.A_AW.meanAmplitude_blSub_active_ipsi(idx);tng.firingField.A_AW.meanAmplitude_blSub_active_ipsi(find(tng.passed.AW.A==1));tng.firingField.X_AW.meanAmplitude_blSub_active_ipsi(idx);tng.firingField.X_AW.meanAmplitude_blSub_active_ipsi(find(tng.passed.AW.X==1))];
    plt = struct(); plt.xlabel = 'Mean amplitude (active)';  plt.xlim = [0,nanmax(temp)]; plt.ylabel = 'X cells'; plt.norm = 'count'; plt.edges = [0:0.05:10];
    plt.dataPoint = tng.firingField.X_AW.meanAmplitude_blSub_active_ipsi(idx);
    hold on
    myHistogram(tng.firingField.X_AW.meanAmplitude_blSub_active_ipsi(find(tng.shuffling.X_AW.significant==1)),p,plt); temp2 = get(gca,'ylim');
    plt.color = p.col.darkGray;
    myHistogram(tng.firingField.X_AW.meanAmplitude_blSub_active_ipsi(find(tng.passed.AW.X==1)),p,plt); temp2 = get(gca,'ylim');
    hold off

    subplot(nrows,ncols,(6-1)*ncols+5)
    ylim([0,nanmax([temp1(2),temp2(2)])])
    subplot(nrows,ncols,(6-1)*ncols+6)
    ylim([0,nanmax([temp1(2),temp2(2)])])

    % Location A
    subplot(nrows,ncols,(7-1)*ncols+5)
    plt = struct(); plt.xlabel = 'Location (s)'; plt.xlim = [0,info.task.trialStructure.tOdour1+info.task.trialStructure.tGap]; plt.ylabel = 'A cells'; plt.norm = 'count'; plt.edges = [0:0.3:5.3];
    plt.dataPoint = tng.firingField.A_AW.peakLocation_s(idx);
    hold on
    myHistogram(tng.firingField.A_AW.peakLocation_s(find(tng.shuffling.A_AW.significant==1)),p,plt); temp1 = get(gca,'ylim');
    plt.color = p.col.darkGray;
    myHistogram(tng.firingField.A_AW.peakLocation_s(find(tng.passed.AW.A==1)),p,plt); temp1 = get(gca,'ylim');
    hold off

    % Location X
    subplot(nrows,ncols,(7-1)*ncols+6)
    plt = struct(); plt.xlabel = 'Location (s)'; plt.xlim = [0,info.task.trialStructure.tOdour1+info.task.trialStructure.tGap]; plt.ylabel = 'X cells'; plt.norm = 'count'; plt.edges = [0:0.3:5.3];
    plt.dataPoint = tng.firingField.X_AW.peakLocation_s(idx);
    hold on
    myHistogram(tng.firingField.X_AW.peakLocation_s(find(tng.shuffling.X_AW.significant==1)),p,plt); temp2 = get(gca,'ylim');
    plt.color = p.col.darkGray;
    myHistogram(tng.firingField.X_AW.peakLocation_s(find(tng.passed.AW.X==1)),p,plt); temp2 = get(gca,'ylim');
    hold off

    subplot(nrows,ncols,(7-1)*ncols+5)
    ylim([0,nanmax([temp1(2),temp2(2)])])
    subplot(nrows,ncols,(7-1)*ncols+6)
    ylim([0,nanmax([temp1(2),temp2(2)])])

    % % Width A
    % subplot(nrows,ncols,(7-1)*ncols+5)
    % plt = struct(); plt.xlabel = 'Width (s)'; plt.xlim = [0,info.task.trialStructure.tOdour1+info.task.trialStructure.tGap]; plt.ylabel = 'A cells'; plt.norm = 'count'; plt.edges = [0:0.3:5.3];
    % plt.dataPoint = tng.firingField.A_AW.width_s(idx);
    % hold on
    % myHistogram(tng.firingField.A_AW.width_s(find(tng.shuffling.A_AW.significant==1)),p,plt); temp1 = get(gca,'ylim');
    % plt.color = p.col.darkGray;
    % myHistogram(tng.firingField.A_AW.width_s(find(tng.passed.AW.A==1)),p,plt); temp1 = get(gca,'ylim');
    % hold off
    % 
    % % Width X
    % subplot(nrows,ncols,(7-1)*ncols+6)
    % plt = struct(); plt.xlabel = 'Width (s)'; plt.xlim = [0,info.task.trialStructure.tOdour1+info.task.trialStructure.tGap]; plt.ylabel = 'X cells'; plt.norm = 'count'; plt.edges = [0:0.3:5.3];
    % plt.dataPoint = tng.firingField.X_AW.width_s(idx);
    % hold on
    % myHistogram(tng.firingField.X_AW.width_s(find(tng.shuffling.X_AW.significant==1)),p,plt); temp2 = get(gca,'ylim');
    % plt.color = p.col.darkGray;
    % myHistogram(tng.firingField.X_AW.width_s(find(tng.passed.AW.X==1)),p,plt); temp2 = get(gca,'ylim');
    % hold off
    % 
    % subplot(nrows,ncols,(7-1)*ncols+5)
    % ylim([0,nanmax([temp1(2),temp2(2)])])
    % subplot(nrows,ncols,(7-1)*ncols+6)
    % ylim([0,nanmax([temp1(2),temp2(2)])])

    % % Peak time variance
    % subplot(nrows,ncols,(8-1)*ncols+5)
    % temp = [tng.firingField.A_AW.peakTimeSd_s(idx);tng.firingField.A_AW.peakTimeSd_s(find(tng.passed.AW.A==1));tng.firingField.A_AW.peakTimeSd_s(idx);tng.firingField.A_AW.peakTimeSd_s(find(tng.passed.AW.X==1))];
    % plt = struct(); plt.xlabel = 'Peak time variance (s)'; plt.xlim = [nanmin(temp),nanmax(temp)]; plt.ylabel = 'A cells'; plt.norm = 'count'; plt.edges = [0:0.05:5.3];
    % plt.dataPoint = tng.firingField.A_AW.peakTimeSd_s(idx);
    % hold on
    % myHistogram(tng.firingField.A_AW.peakTimeSd_s(find(tng.shuffling.A_AW.significant==1)),p,plt); temp1 = get(gca,'ylim');
    % plt.color = p.col.darkGray;
    % myHistogram(tng.firingField.A_AW.peakTimeSd_s(find(tng.passed.AW.A==1)),p,plt); temp1 = get(gca,'ylim');
    % hold off
    % 
    % % Peak time variance
    % subplot(nrows,ncols,(8-1)*ncols+6)
    % temp = [tng.firingField.A_AW.peakTimeSd_s(idx);tng.firingField.A_AW.peakTimeSd_s(find(tng.passed.AW.A==1));tng.firingField.A_AW.peakTimeSd_s(idx);tng.firingField.A_AW.peakTimeSd_s(find(tng.passed.AW.X==1))];
    % plt = struct(); plt.xlabel = 'Peak time variance (s)'; plt.xlim = [nanmin(temp),nanmax(temp)]; plt.ylabel = 'X cells'; plt.norm = 'count'; plt.edges = [0:0.05:5.3];
    % plt.dataPoint = tng.firingField.A_AW.peakTimeSd_s(idx);
    % hold on
    % myHistogram(tng.firingField.A_AW.peakTimeSd_s(find(tng.shuffling.X_AW.significant==1)),p,plt); temp2 = get(gca,'ylim');
    % plt.color = p.col.darkGray;
    % myHistogram(tng.firingField.A_AW.peakTimeSd_s(find(tng.passed.AW.X==1)),p,plt); temp2 = get(gca,'ylim');
    % hold off

    % Amplitude A
    subplot(nrows,ncols,(8-1)*ncols+5)
    temp = [tng.firingField.A_AW.meanAmplitude_blSub_ipsi(idx);tng.firingField.A_AW.meanAmplitude_blSub_ipsi(find(tng.passed.AW.A==1));tng.firingField.X_AW.meanAmplitude_blSub_ipsi(idx);tng.firingField.X_AW.meanAmplitude_blSub_ipsi(find(tng.passed.AW.X==1))];
    plt = struct(); plt.xlabel = 'Mean amplitude (all)';  plt.xlim = [0,nanmax(temp)]; plt.ylabel = 'A cells'; plt.norm = 'count'; plt.edges = [0:0.05:10];
    plt.dataPoint = tng.firingField.A_AW.meanAmplitude_blSub_ipsi(idx);
    hold on
    myHistogram(tng.firingField.A_AW.meanAmplitude_blSub_ipsi(find(tng.shuffling.A_AW.significant==1)),p,plt); temp1 = get(gca,'ylim');
    plt.color = p.col.darkGray;
    myHistogram(tng.firingField.A_AW.meanAmplitude_blSub_ipsi(find(tng.passed.AW.A==1)),p,plt); temp1 = get(gca,'ylim');
    hold off

    % Amplitude X
    subplot(nrows,ncols,(8-1)*ncols+6)
    temp = [tng.firingField.A_AW.meanAmplitude_blSub_ipsi(idx);tng.firingField.A_AW.meanAmplitude_blSub_ipsi(find(tng.passed.AW.A==1));tng.firingField.X_AW.meanAmplitude_blSub_ipsi(idx);tng.firingField.X_AW.meanAmplitude_blSub_ipsi(find(tng.passed.AW.X==1))];
    plt = struct(); plt.xlabel = 'Mean amplitude (all)';  plt.xlim = [0,nanmax(temp)]; plt.ylabel = 'X cells'; plt.norm = 'count'; plt.edges = [0:0.05:10];
    plt.dataPoint = tng.firingField.X_AW.meanAmplitude_blSub_ipsi(idx);
    hold on
    myHistogram(tng.firingField.X_AW.meanAmplitude_blSub_ipsi(find(tng.shuffling.X_AW.significant==1)),p,plt); temp2 = get(gca,'ylim');
    plt.color = p.col.darkGray;
    myHistogram(tng.firingField.X_AW.meanAmplitude_blSub_ipsi(find(tng.passed.AW.X==1)),p,plt); temp2 = get(gca,'ylim');
    hold off

    subplot(nrows,ncols,(8-1)*ncols+5)
    ylim([0,nanmax([temp1(2),temp2(2)])])
    subplot(nrows,ncols,(8-1)*ncols+6)
    ylim([0,nanmax([temp1(2),temp2(2)])])
end


%% Properties


%% Return

% ; dFF (o=100ms, 10ms bins)

if compact_version
    subplot(nrows,ncols,sort([1:6]))
else
    subplot(nrows,ncols,sort([([1:2]-1)*ncols+3,([1:2]-1)*ncols+4]))
end
%set(gca,'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[],'xcolor',[1,1,1],'ycolor',[1,1,1])
if info.stimSession
    title([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,newline,...
        'idx=',num2str(idx),', iscell=',num2str(this_iscell),', seqA=',num2str(this_seqA),', seqX=',num2str(this_seqX),...
        '']);
else
    title([info.animal,'-',info.date,'-d',num2str(info.expDay),'-','nostim',newline,...
        'idx=',num2str(idx),', iscell=',num2str(this_iscell),', seqA=',num2str(this_seqA),', seqX=',num2str(this_seqX),...
        '']);
end
drawnow;
%end


