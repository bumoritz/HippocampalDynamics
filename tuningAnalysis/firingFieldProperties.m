function firingField = firingFieldProperties(traces,avgTraces,normAvgTraces,timewindow,ipsi,contra,prop,p)
%traces = traces_all; avgTraces = avgTraces_all; normAvgTraces = normAvgTraces_all; timewindow = 'AW'; ipsi = 'A'; contra = 'X';
%traces = traces_60t{1}; avgTraces = avgTraces_60t{1}; normAvgTraces = normAvgTraces_60t{1}; timewindow = 'AW'; ipsi = 'A'; contra = 'X';

%% Preparations

% extract relevant traces, avgTraces and normAvgTraces
these_traces_ipsi = traces.(ipsi);
these_traces_contra = traces.(contra);
these_avgTraces_ipsi = avgTraces.(ipsi);
these_avgTraces_contra = avgTraces.(contra);
if strcmp(timewindow,'AW')
    these_normAvgTraces_ipsi = normAvgTraces.firstOdour.(ipsi);
    these_normAvgTraces_contra = normAvgTraces.firstOdour.(contra);
elseif strcmp(timewindow,'odour1')
    these_normAvgTraces_ipsi = normAvgTraces.firstOdour.(ipsi);
    these_normAvgTraces_contra = normAvgTraces.firstOdour.(contra);
elseif strcmp(timewindow,'odour2')
    these_normAvgTraces_ipsi = normAvgTraces.secondOdour.(ipsi);
    these_normAvgTraces_contra = normAvgTraces.secondOdour.(contra);
elseif strcmp(timewindow,'reward')
    these_normAvgTraces_ipsi = normAvgTraces.reward.(ipsi);
    these_normAvgTraces_contra = normAvgTraces.reward.(contra);
end

% extract relevant time windows
if strcmp(timewindow,'AW')
    this_window_core = p.general.bins_analysisWindow;
    this_window_base = p.general.bins_baselineWindow;
elseif strcmp(timewindow,'odour1')
    this_window_core = p.general.bins_odour1Window;
    this_window_base = p.general.bins_odour1Window_bl;
elseif strcmp(timewindow,'odour2')
    this_window_core = p.general.bins_odour2Window;
    this_window_base = p.general.bins_odour2Window_bl;
elseif strcmp(timewindow,'reward')
    this_window_core = p.general.bins_rewardWindow;
    this_window_base = p.general.bins_rewardWindow_bl;
end


%% Identify firing field peak

% identify firing field peak
[firingField.peakAmplitude,temp] = nanmax(these_avgTraces_ipsi(:,this_window_core),[],2);
firingField.peakLocation = temp + this_window_core(1)-1;
firingField.peakLocation_s = p.general.t_binned(firingField.peakLocation)';
firingField.peakLocation(~prop.iscell==1) = NaN;
firingField.peakLocation_s(~prop.iscell==1) = NaN;

% early vs. late sequence cell
if strcmp(timewindow,'AW')
    firingField.binary.early = double(firingField.peakLocation_s <= p.tng.earlyVsLateThreshold);
    firingField.binary.late = double(firingField.peakLocation_s > p.tng.earlyVsLateThreshold);
    firingField.binary.early(~prop.iscell==1) = NaN;
    firingField.binary.late(~prop.iscell==1) = NaN;
end

% calculate peak-time variability
firingField.peakTimeSd_s = nan(prop.numRois,1);
for i=1:prop.numRois
    [~,temp] = nanmax(squeeze(these_traces_ipsi(i,this_window_core,:)),[],1);
    firingField.peakTimeSd_s(i) = nanstd(p.general.t_binned(temp+this_window_core(1)-1));
end
firingField.peakTimeSd_s(~prop.iscell==1) = NaN;


%% Calculate baseline and sd threshold

% baseline
all_traces = cat(3,these_traces_ipsi,these_traces_contra);
firingField.baseline_trialwise = squeeze(nanmean(all_traces(:,this_window_base,:),2));
firingField.baseline_meanOfTrialwise = nanmean(firingField.baseline_trialwise,2);
firingField.baseline_sdOfTrialwise = nanstd(firingField.baseline_trialwise,[],2);

% sd above baseline threshold
if strcmp(p.tng.sdAboveBaselineType,'sd_of_entire_recording')
    firingField.activeInFiringFieldThreshold = firingField.baseline_meanOfTrialwise + p.tng.activeInFiringFieldSd * prop.nf_sd;
elseif strcmp(p.tng.sdAboveBaselineType,'sd_of_trialwise_baseline')
    firingField.activeInFiringFieldThreshold = firingField.baseline_meanOfTrialwise + p.tng.activeInFiringFieldSd * firingField.baseline_sdOfTrialwise;
end

% correct peak amplitude by baseline subtraction
firingField.peakAmplitude_blSub = firingField.peakAmplitude - firingField.baseline_meanOfTrialwise;


%% Firing field extent

% identify firing field boundaries
firingField.boundaries = nan(prop.numRois,2);
firingField.width_s = nan(prop.numRois,1);
temp = these_avgTraces_ipsi > firingField.baseline_meanOfTrialwise + p.tng.firingFieldBoundaries*firingField.peakAmplitude_blSub;
for i=1:prop.numRois
    n=0;
    try
        while temp(i,firingField.peakLocation(i)-n)
            firingField.boundaries(i,1) = firingField.peakLocation(i)-n;
            n=n+1;
        end
    catch
    end
    n=0;
    try
        while temp(i,firingField.peakLocation(i)+n)
            firingField.boundaries(i,2) = firingField.peakLocation(i)+n;
            n=n+1;
        end
    catch
    end
    if ~isnan(firingField.boundaries(i,1)) && ~isnan(firingField.boundaries(i,2))
        firingField.width_s(i) = p.general.t_binned(firingField.boundaries(i,2))-p.general.t_binned(firingField.boundaries(i,1));
    end
end
firingField.boundaries_wdw = firingField.boundaries;
firingField.boundaries_wdw(firingField.boundaries_wdw>max(this_window_core)) = max(this_window_core);
firingField.boundaries_wdw(firingField.boundaries_wdw<min(this_window_core)) = min(this_window_core);

% identify firing field centre
firingField.centreLocation_s = nan(prop.numRois,1);
for i=1:prop.numRois
    try
        firingField.centreLocation_s(i) = nanmean([p.general.t_binned(firingField.boundaries(i,1)),p.general.t_binned(firingField.boundaries(i,2))]);
    catch
    end
end

% peak in window
firingField.binary.peakInWindow = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField.boundaries_wdw(i,1)) && ~isnan(firingField.boundaries_wdw(i,2))
        firingField.binary.peakInWindow(i) = double(nanmax(these_avgTraces_ipsi(i,firingField.boundaries_wdw(i,1):firingField.boundaries_wdw(i,2))) >= ...
            nanmax(these_avgTraces_ipsi(i,firingField.boundaries(i,1):firingField.boundaries(i,2))));
    end
end
firingField.binary.peakInWindow(~prop.iscell==1) = NaN;


%% Firing field properties (reliability, mean amplitude, selectivity)

% reliability index
firingField.activeTrials_ipsi = nan(prop.numRois,size(these_traces_ipsi,3));
firingField.activeTrials_contra = nan(prop.numRois,size(these_traces_contra,3));
firingField.trialwiseMeanActivity_ipsi = nan(prop.numRois,1);
firingField.trialwiseMeanActivity_contra = nan(prop.numRois,1);
firingField.trialwiseSdActivity_ipsi = nan(prop.numRois,1);
firingField.trialwiseSdActivity_contra = nan(prop.numRois,1);
firingField.comLocation_s = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField.boundaries_wdw(i,1)) && ~isnan(firingField.boundaries_wdw(i,2))
        firingField.activeTrials_ipsi(i,:) = double(nanmean(squeeze(these_traces_ipsi(i,firingField.boundaries_wdw(i,1):firingField.boundaries_wdw(i,2),:)),1) > firingField.activeInFiringFieldThreshold(i));
        firingField.activeTrials_contra(i,:) = double(nanmean(squeeze(these_traces_contra(i,firingField.boundaries_wdw(i,1):firingField.boundaries_wdw(i,2),:)),1) > firingField.activeInFiringFieldThreshold(i));
        firingField.trialwiseMeanActivity_ipsi(i) = nanmean(nanmean(squeeze(these_traces_ipsi(i,firingField.boundaries_wdw(i,1):firingField.boundaries_wdw(i,2),:)),1));
        firingField.trialwiseMeanActivity_contra(i) = nanmean(nanmean(squeeze(these_traces_contra(i,firingField.boundaries_wdw(i,1):firingField.boundaries_wdw(i,2),:)),1));
        firingField.trialwiseSdActivity_ipsi(i) = nanstd(nanmean(squeeze(these_traces_ipsi(i,firingField.boundaries_wdw(i,1):firingField.boundaries_wdw(i,2),:)),1));
        firingField.trialwiseSdActivity_contra(i) = nanstd(nanmean(squeeze(these_traces_contra(i,firingField.boundaries_wdw(i,1):firingField.boundaries_wdw(i,2),:)),1));
        
        % calculate centre-of-mass
        temp = nanmean(these_traces_ipsi(i,firingField.boundaries_wdw(i,1):firingField.boundaries_wdw(i,2),:),3) + 1000; % + large number to avoid error in next line
        firingField.comLocation_s(i) = p.general.t_binned(firingField.boundaries_wdw(i,1)+round(nansum((1:length(temp)).*temp)/nansum(temp)) - 1);
    end
end
firingField.activeTrials_ipsi(~prop.iscell==1,:) = NaN;
firingField.activeTrials_contra(~prop.iscell==1,:) = NaN;
firingField.activationProbability_ipsi = nanmean(firingField.activeTrials_ipsi,2);
firingField.activationProbability_contra = nanmean(firingField.activeTrials_contra,2);
firingField.activationProbability_ipsi(~prop.iscell==1) = NaN;
firingField.activationProbability_contra(~prop.iscell==1) = NaN;
firingField.binary.reliability_ipsi = double(firingField.activationProbability_ipsi > p.tng.reliabilityThreshold);
firingField.binary.reliability_contra = double(firingField.activationProbability_contra > p.tng.reliabilityThreshold);
firingField.binary.reliability_ipsi(~prop.iscell==1) = NaN;
firingField.binary.reliability_contra(~prop.iscell==1) = NaN;

% mean activity in firing field (incl. bl-sub)
firingField.meanAmplitude_ipsi = nan(prop.numRois,1);
firingField.meanAmplitude_contra = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField.boundaries_wdw(i,1)) && ~isnan(firingField.boundaries_wdw(i,2))
        firingField.meanAmplitude_ipsi(i) = nanmean(these_avgTraces_ipsi(i,firingField.boundaries_wdw(i,1):firingField.boundaries_wdw(i,2)));
        firingField.meanAmplitude_contra(i) = nanmean(these_avgTraces_contra(i,firingField.boundaries_wdw(i,1):firingField.boundaries_wdw(i,2)));
    end
end
firingField.meanAmplitude_blSub_ipsi = firingField.meanAmplitude_ipsi - firingField.baseline_meanOfTrialwise;
firingField.meanAmplitude_blSub_contra = firingField.meanAmplitude_contra - firingField.baseline_meanOfTrialwise;
firingField.binary.meanAmplitude_blSub_ipsi = double(firingField.meanAmplitude_blSub_ipsi > p.tng.amplitudeThreshold);
firingField.binary.meanAmplitude_blSub_contra = double(firingField.meanAmplitude_blSub_contra > p.tng.amplitudeThreshold);
firingField.binary.meanAmplitude_blSub_ipsi(~prop.iscell==1) = NaN;
firingField.binary.meanAmplitude_blSub_contra(~prop.iscell==1) = NaN;

% mean bl-sub activity in firing field in active trials only
firingField.meanAmplitude_blSub_active_ipsi = nan(prop.numRois,1);
firingField.meanAmplitude_blSub_active_contra = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField.boundaries_wdw(i,1)) && ~isnan(firingField.boundaries_wdw(i,2))
        temp = squeeze(these_traces_ipsi(i,firingField.boundaries_wdw(i,1):firingField.boundaries_wdw(i,2),find(firingField.activeTrials_ipsi(i,:)==1)));
        firingField.meanAmplitude_blSub_active_ipsi(i) = nanmean(nanmean(temp,1));
        temp = squeeze(these_traces_contra(i,firingField.boundaries_wdw(i,1):firingField.boundaries_wdw(i,2),find(firingField.activeTrials_contra(i,:)==1)));
        firingField.meanAmplitude_blSub_active_contra(i) = nanmean(nanmean(temp,1));
    end
end

% mean norm activity in firing field
firingField.meanNormAmplitude_ipsi = nan(prop.numRois,1);
firingField.meanNormAmplitude_contra = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField.boundaries_wdw(i,1)) && ~isnan(firingField.boundaries_wdw(i,2))
        firingField.meanNormAmplitude_ipsi(i) = nanmean(these_normAvgTraces_ipsi(i,firingField.boundaries_wdw(i,1):firingField.boundaries_wdw(i,2)));
        firingField.meanNormAmplitude_contra(i) = nanmean(these_normAvgTraces_contra(i,firingField.boundaries_wdw(i,1):firingField.boundaries_wdw(i,2)));
    end
end

% selectivity index
firingField.selectivity_ipsi = (firingField.meanNormAmplitude_ipsi - firingField.meanNormAmplitude_contra) ./ (firingField.meanNormAmplitude_ipsi + firingField.meanNormAmplitude_contra);
firingField.selectivity_contra = (firingField.meanNormAmplitude_contra - firingField.meanNormAmplitude_ipsi) ./ (firingField.meanNormAmplitude_ipsi + firingField.meanNormAmplitude_contra);
temp = [firingField.selectivity_ipsi;firingField.selectivity_contra];
if any(temp > 1 | temp < -1)
    warning('Absolute selectivity index is sometimes greater than 1.')
end


%% Trial-wise amplitude statistics

firingField.diffBl_base = nan(prop.numRois,size(these_traces_ipsi,3));
firingField.diffBl_resp = nan(prop.numRois,size(these_traces_ipsi,3));
firingField.diffBl_diff = nan(prop.numRois,size(these_traces_ipsi,3));
firingField.diffBl_p = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField.boundaries_wdw(i,1)) && ~isnan(firingField.boundaries_wdw(i,2))
        firingField.diffBl_base(i,:) = squeeze(nanmean(these_traces_ipsi(i,p.general.bins_baselineWindow,:),2))';
        firingField.diffBl_resp(i,:) = squeeze(nanmax(these_traces_ipsi(i,firingField.boundaries_wdw(i,1):firingField.boundaries_wdw(i,2),:),[],2))'; % max in firing field
        firingField.diffBl_diff(i,:) = firingField.diffBl_resp(i,:) - firingField.diffBl_base(i,:);
        firingField.diffBl_p(i) = signrank(firingField.diffBl_base(i,:)',firingField.diffBl_resp(i,:)');
    end
end
[firingField.diffBl_fdr,firingField.diffBl_q,firingField.diffBl_priori] = mafdr(firingField.diffBl_p); % double-checked that it works correctly with NaNs
if firingField.diffBl_priori>=1
    warning('The priori is greater than or equal to 1.')
end
firingField.diffBl = double(firingField.diffBl_q < p.tng.qThreshold);
firingField.diffBl(~prop.iscell==1) = NaN;
firingField.diffBl_pos = zeros(prop.numRois,1);
firingField.diffBl_pos(intersect(find(nanmean(firingField.diffBl_diff,2)>0),find(firingField.diffBl==1))) = 1;
firingField.diffBl_pos(~prop.iscell==1) = NaN;


%% Return

firingField.binary = orderfields(firingField.binary);
firingField = orderfields(firingField);
end

