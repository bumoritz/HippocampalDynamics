function [sequenceCellActivity] = sequenceCellActivityByOutcome(traces,avgTraces,normAvgTraces,firingField,idcs,timewindow,ipsi,contra,prop,p)
%i=1; traces=traces_cv{i}; avgTraces=avgTraces_cv{i}; normAvgTraces=normAvgTraces_cv{i}; firingField=firingField_cv{i}.A_AW; idcs=find(passed_cv{i}.AW.A==1); timewindow='AW'; ipsi='A'; contra='X';
%i=1; traces=traces_cv{i}; avgTraces=avgTraces_cv{i}; normAvgTraces=normAvgTraces_cv{i}; firingField=firingField_cv{i}.A_AW; idcs=find(passed_cv{i}.AW.A==1); timewindow='FF'; ipsi='A'; contra='X';

%% Preparations

all_fields = fields(traces);

% extract relevant time windows
if strcmp(timewindow,'AW')
    this_window = nan(prop.numRois,2);
    for i=1:prop.numRois
        this_window(i,:) = [p.general.bins_analysisWindow(1),p.general.bins_analysisWindow(end)];
    end
	this_window(~prop.iscell==1,:) = NaN;
elseif strcmp(timewindow,'FF')
	this_window = firingField.boundaries_wdw;
end


%% Peak amplitude

for f=1:length(all_fields)
    this_field = all_fields{f};
    
    sequenceCellActivity.(this_field).peakAmplitude_blSub = nan(prop.numRois,1);
    for i=1:prop.numRois
        if ~isnan(this_window(i,1)) && ~isnan(this_window(i,2))
            sequenceCellActivity.(this_field).peakAmplitude_blSub(i) = nanmax(avgTraces.(this_field)(i,this_window(i,1):this_window(i,2))) - firingField.baseline_meanOfTrialwise(i);
        end
    end
end

figure;
yline(0,':');
hold on
v = violinplot([sequenceCellActivity.A_correct_train.peakAmplitude_blSub(idcs),...
    sequenceCellActivity.A_correct_test.peakAmplitude_blSub(idcs),...
    sequenceCellActivity.A_incorrect.peakAmplitude_blSub(idcs),...
    sequenceCellActivity.X_correct_train.peakAmplitude_blSub(idcs),...
    sequenceCellActivity.X_correct_test.peakAmplitude_blSub(idcs),...
    sequenceCellActivity.X_incorrect.peakAmplitude_blSub(idcs)],...
    {char('A-corr-train'),char('A-corr-test'),char('A-incorr'),char('X-corr-train'),char('X-corr-test'),char('X-incorr')});
ylim([-inf,inf])
ylabel(['Peak amplitude of sequence A cells'])
title('Peak amplitude by outcome')


%% Mean amplitude

for f=1:length(all_fields)
    this_field = all_fields{f};
    
    sequenceCellActivity.(this_field).meanAmplitude_blSub = nan(prop.numRois,1);
    for i=1:prop.numRois
        if ~isnan(this_window(i,1)) && ~isnan(this_window(i,2))
            sequenceCellActivity.(this_field).meanAmplitude_blSub(i) = nanmean(avgTraces.(this_field)(i,this_window(i,1):this_window(i,2))) - firingField.baseline_meanOfTrialwise(i);
        end
    end
end

figure;
yline(0,':');
hold on
v = violinplot([sequenceCellActivity.A_correct_train.meanAmplitude_blSub(idcs),...
    sequenceCellActivity.A_correct_test.meanAmplitude_blSub(idcs),...
    sequenceCellActivity.A_incorrect.meanAmplitude_blSub(idcs),...
    sequenceCellActivity.X_correct_train.meanAmplitude_blSub(idcs),...
    sequenceCellActivity.X_correct_test.meanAmplitude_blSub(idcs),...
    sequenceCellActivity.X_incorrect.meanAmplitude_blSub(idcs)],...
    {char('A-corr-train'),char('A-corr-test'),char('A-incorr'),char('X-corr-train'),char('X-corr-test'),char('X-incorr')});
ylim([-inf,inf])
ylabel(['Mean amplitude of sequence A cells'])
title('Mean amplitude by outcome')


%% Reliability

for f=1:length(all_fields)
    this_field = all_fields{f};
    
    sequenceCellActivity.(this_field).activeTrials = nan(prop.numRois,size(traces.(this_field),3));
    for i=1:prop.numRois
        if ~isnan(this_window(i,1)) && ~isnan(this_window(i,2))
            sequenceCellActivity.(this_field).activeTrials(i,:) = double(nanmean(squeeze(traces.(this_field)(i,this_window(i,1):this_window(i,2),:)),1) > firingField.activeInFiringFieldThreshold(i));
        end
    end
    sequenceCellActivity.(this_field).activeTrials(~prop.iscell==1,:) = NaN;
    sequenceCellActivity.(this_field).activationProbability = nanmean(sequenceCellActivity.(this_field).activeTrials,2);
    sequenceCellActivity.(this_field).reliability = double(sequenceCellActivity.(this_field).activationProbability > p.ett.reliabilityThreshold);
    sequenceCellActivity.(this_field).reliability(~prop.iscell==1) = NaN;
end

figure;
v = violinplot([sequenceCellActivity.A_correct_train.activationProbability(idcs),...
    sequenceCellActivity.A_correct_test.activationProbability(idcs),...
    sequenceCellActivity.A_incorrect.activationProbability(idcs),...
    sequenceCellActivity.X_correct_train.activationProbability(idcs),...
    sequenceCellActivity.X_correct_test.activationProbability(idcs),...
    sequenceCellActivity.X_incorrect.activationProbability(idcs)],...
    {char('A-corr-train'),char('A-corr-test'),char('A-incorr'),char('X-corr-train'),char('X-corr-test'),char('X-incorr')});
ylim([0,inf])
ylabel(['Activation probability of sequence A cells'])
title('Reliability by outcome')


%% Peak-time variability

% wrong!
% for f=1:length(all_fields)
%     this_field = all_fields{f};
%     
%     sequenceCellActivity.(this_field).peakTimeSd_s = nan(prop.numRois,1);
%     for i=1:prop.numRois
%         if ~isnan(this_window(i,1)) && ~isnan(this_window(i,2)) && this_window(i,1)~=this_window(i,2)
%             [~,temp] = nanmax(squeeze(traces.(this_field)(i,this_window(i,1):this_window(i,2),:)),[],1);
%             temp
%             sequenceCellActivity.(this_field).peakTimeSd_s(i) = nanstd(p.general.t_binned(temp+this_window(i,1)-1));
%         end
%     end
%     sequenceCellActivity.(this_field).peakTimeSd_s(~prop.iscell==1) = NaN;
% end
% 
% figure;
% v = violinplot([sequenceCellActivity.A_correct_train.peakTimeSd_s(idcs),...
%     sequenceCellActivity.A_correct_test.peakTimeSd_s(idcs),...
%     sequenceCellActivity.A_incorrect.peakTimeSd_s(idcs),...
%     sequenceCellActivity.X_correct_train.peakTimeSd_s(idcs),...
%     sequenceCellActivity.X_correct_test.peakTimeSd_s(idcs),...
%     sequenceCellActivity.X_incorrect.peakTimeSd_s(idcs)],...
%     {char('A-corr-train'),char('A-corr-test'),char('A-incorr'),char('X-corr-train'),char('X-corr-test'),char('X-incorr')});
% ylim([0,inf])
% ylabel(['Standard deviation of peak time of sequence A cells (s)'])
% title('Peak-time variability by outcome')


%% Return

end







