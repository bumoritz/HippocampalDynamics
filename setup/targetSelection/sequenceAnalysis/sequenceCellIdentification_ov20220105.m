function sca = sequenceCellIdentification_ov20220105(nft_binned,trials,prop,iscell,p)

%% Preparations

% trial type conjunctions
[trials.trials_AB,trials.trials_AB_inA,~] = intersect(trials.trials_A,trials.trials_B);
[trials.trials_AY,trials.trials_AY_inA,~] = intersect(trials.trials_A,trials.trials_Y);
[trials.trials_XY,trials.trials_XY_inX,~] = intersect(trials.trials_X,trials.trials_Y);
[trials.trials_XB,trials.trials_XB_inX,~] = intersect(trials.trials_X,trials.trials_B);

% NaN non-iscells
if p.sca.nannoniscell
    temp = nft_binned;
    nft_binned = nan(size(nft_binned));
    nft_binned(find(iscell),:,:) = temp(find(iscell),:,:);
end

% traces
traces.traces_A = nft_binned(:,:,trials.trials_A);
traces.traces_X = nft_binned(:,:,trials.trials_X);
traces.traces_B = nft_binned(:,:,trials.trials_B);
traces.traces_Y = nft_binned(:,:,trials.trials_Y);
traces.traces_AB = nft_binned(:,:,trials.trials_AB);
traces.traces_AY = nft_binned(:,:,trials.trials_AY);
traces.traces_XY = nft_binned(:,:,trials.trials_XY);
traces.traces_XB = nft_binned(:,:,trials.trials_XB);

% avgTraces
temp = fields(traces);
for i=1:length(temp)
    temp2 = temp{i};
    avgTraces.(['avgTraces_',temp2(8:end)]) = nanmean(traces.(temp2),3);
end

% normAvgTraces (for first odour selectivity index)
temp = [avgTraces.avgTraces_A(:);avgTraces.avgTraces_X(:)];
this_lower = nanmin(temp(:));
this_upper = nanmax(temp(:));
normAvgTraces.normAvgTraces_A = (avgTraces.avgTraces_A-this_lower) ./ (this_upper-this_lower);
normAvgTraces.normAvgTraces_X = (avgTraces.avgTraces_X-this_lower) ./ (this_upper-this_lower);

% normAvgTraces (for second odour selectivity index)
temp = [avgTraces.avgTraces_AB(:);avgTraces.avgTraces_XB(:);avgTraces.avgTraces_AY(:);avgTraces.avgTraces_XY(:)];
this_lower = nanmin(temp(:));
this_upper = nanmax(temp(:));
normAvgTraces.normAvgTraces_AB = (avgTraces.avgTraces_AB-this_lower) ./ (this_upper-this_lower);
normAvgTraces.normAvgTraces_XB = (avgTraces.avgTraces_XB-this_lower) ./ (this_upper-this_lower);
normAvgTraces.normAvgTraces_XY = (avgTraces.avgTraces_XY-this_lower) ./ (this_upper-this_lower);
normAvgTraces.normAvgTraces_AY = (avgTraces.avgTraces_AY-this_lower) ./ (this_upper-this_lower);


%% Identify tuning

% shuffling criterion
sca.shuffling_A = shufflingCriterion_ov20220105(traces.traces_A(:,prop.analysisWindow,:),iscell,p);
sca.shuffling_X = shufflingCriterion_ov20220105(traces.traces_X(:,prop.analysisWindow,:),iscell,p);

% identify firing field peak
[firingField_A.peakAmplitude,temp] = nanmax(avgTraces.avgTraces_A(:,prop.analysisWindow),[],2);
firingField_A.peakLocation = temp + prop.frames_pre_binned;
firingField_A.peakLocation_s = prop.t_binned(firingField_A.peakLocation)';
[firingField_X.peakAmplitude,temp] = nanmax(avgTraces.avgTraces_X(:,prop.analysisWindow),[],2);
firingField_X.peakLocation = temp + prop.frames_pre_binned;
firingField_X.peakLocation_s = prop.t_binned(firingField_X.peakLocation)';

% early vs. late sequence cell
firingField_A.early = firingField_A.peakLocation_s <= p.sca.earlyVsLateThreshold;
firingField_A.late = firingField_A.peakLocation_s > p.sca.earlyVsLateThreshold;
firingField_X.early = firingField_X.peakLocation_s <= p.sca.earlyVsLateThreshold;
firingField_X.late = firingField_X.peakLocation_s > p.sca.earlyVsLateThreshold;

% calculate peak-time variability
sca.peakTimeSd_s_A = nan(prop.numRois,1);
sca.peakTimeSd_s_X = nan(prop.numRois,1);
for i=1:prop.numRois
    [~,temp] = nanmax(squeeze(traces.traces_A(i,prop.analysisWindow,:)),[],1);
    sca.peakTimeSd_s_A(i) = nanstd(prop.t_binned(temp+prop.frames_pre_binned));
    [~,temp] = nanmax(squeeze(traces.traces_X(i,prop.analysisWindow,:)),[],1);
    sca.peakTimeSd_s_X(i) = nanstd(prop.t_binned(temp+prop.frames_pre_binned));
end


%% Firing field extent

% calculate baseline and threshold activity levels (new version)
temp = squeeze(nanmean(traces.traces_A(:,prop.baselineWindow,:),2));
firingField_A.baseline = nanmean(temp,2);
firingField_A.activeInFiringFieldThreshold = firingField_A.baseline + p.sca.activeInFiringFieldSd * nanstd(temp,[],2);
temp = squeeze(nanmean(traces.traces_X(:,prop.baselineWindow,:),2));
firingField_X.baseline = nanmean(temp,2);
firingField_X.activeInFiringFieldThreshold = firingField_X.baseline + p.sca.activeInFiringFieldSd * nanstd(temp,[],2);

% correct peak amplitude by baseline subtraction
firingField_A.peakAmplitude_blSub = firingField_A.peakAmplitude - firingField_A.baseline;
firingField_X.peakAmplitude_blSub = firingField_X.peakAmplitude - firingField_X.baseline;

% identify firing field boundaries
firingField_A.boundaries = nan(prop.numRois,2);
firingField_A.width_s = nan(prop.numRois,1);
temp = avgTraces.avgTraces_A > firingField_A.baseline + p.sca.firingFieldBoundaries*firingField_A.peakAmplitude_blSub;
for i=1:prop.numRois
    n=0;
    try
        while temp(i,firingField_A.peakLocation(i)-n)
            firingField_A.boundaries(i,1) = firingField_A.peakLocation(i)-n;
            n=n+1;
        end
    catch
    end
    n=0;
    try
        while temp(i,firingField_A.peakLocation(i)+n)
            firingField_A.boundaries(i,2) = firingField_A.peakLocation(i)+n;
            n=n+1;
        end
    catch
    end
    if ~isnan(firingField_A.boundaries(i,1)) && ~isnan(firingField_A.boundaries(i,2))
        firingField_A.width_s(i) = prop.t_binned(firingField_A.boundaries(i,2))-prop.t_binned(firingField_A.boundaries(i,1));
    end
end
firingField_A.boundaries_AW = firingField_A.boundaries;
firingField_A.boundaries_AW(firingField_A.boundaries_AW>max(prop.analysisWindow)) = max(prop.analysisWindow);
firingField_A.boundaries_AW(firingField_A.boundaries_AW<min(prop.analysisWindow)) = min(prop.analysisWindow);
firingField_X.boundaries = nan(prop.numRois,2);
firingField_X.width_s = nan(prop.numRois,1);
temp = avgTraces.avgTraces_X > firingField_X.baseline + p.sca.firingFieldBoundaries*firingField_X.peakAmplitude_blSub;
for i=1:prop.numRois
    n=0;
    try
        while temp(i,firingField_X.peakLocation(i)-n)
            firingField_X.boundaries(i,1) = firingField_X.peakLocation(i)-n;
            n=n+1;
        end
    catch
    end
    n=0;
    try
        while temp(i,firingField_X.peakLocation(i)+n)
            firingField_X.boundaries(i,2) = firingField_X.peakLocation(i)+n;
            n=n+1;
        end
    catch
    end
    if ~isnan(firingField_X.boundaries(i,1)) && ~isnan(firingField_X.boundaries(i,2))
        firingField_X.width_s(i) = prop.t_binned(firingField_X.boundaries(i,2))-prop.t_binned(firingField_X.boundaries(i,1));
    end
end
firingField_X.boundaries_AW = firingField_X.boundaries;
firingField_X.boundaries_AW(firingField_X.boundaries_AW>max(prop.analysisWindow)) = max(prop.analysisWindow);
firingField_X.boundaries_AW(firingField_X.boundaries_AW<min(prop.analysisWindow)) = min(prop.analysisWindow);

% peak in analysis window
firingField_A.peakInAW = false(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField_A.boundaries_AW(i,1)) && ~isnan(firingField_A.boundaries_AW(i,2))
        firingField_A.peakInAW(i) = nanmax(avgTraces.avgTraces_A(i,firingField_A.boundaries_AW(i,1):firingField_A.boundaries_AW(i,2))) >= ...
            nanmax(avgTraces.avgTraces_A(i,firingField_A.boundaries(i,1):firingField_A.boundaries(i,2)));
    end
end
firingField_X.peakInAW = false(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField_X.boundaries_AW(i,1)) && ~isnan(firingField_X.boundaries_AW(i,2))
        firingField_X.peakInAW(i) = nanmax(avgTraces.avgTraces_X(i,firingField_X.boundaries_AW(i,1):firingField_X.boundaries_AW(i,2))) >= ...
            nanmax(avgTraces.avgTraces_X(i,firingField_X.boundaries(i,1):firingField_X.boundaries(i,2)));
    end
end


%% Firing field properties

% reliability index
firingField_A.activeTrials_A = false(prop.numRois,length(trials.trials_A));
firingField_A.activeTrials_X = false(prop.numRois,length(trials.trials_X));
firingField_A.trialwiseMeanActivity_A = nan(prop.numRois,1);
firingField_A.trialwiseMeanActivity_X = nan(prop.numRois,1);
firingField_A.trialwiseSdActivity_A = nan(prop.numRois,1);
firingField_A.trialwiseSdActivity_X = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField_A.boundaries_AW(i,1)) && ~isnan(firingField_A.boundaries_AW(i,2))
        firingField_A.activeTrials_A(i,:) = nanmean(squeeze(traces.traces_A(i,firingField_A.boundaries_AW(i,1):firingField_A.boundaries_AW(i,2),:)),1) > firingField_A.activeInFiringFieldThreshold(i);
        firingField_A.activeTrials_X(i,:) = nanmean(squeeze(traces.traces_X(i,firingField_A.boundaries_AW(i,1):firingField_A.boundaries_AW(i,2),:)),1) > firingField_A.activeInFiringFieldThreshold(i);
        firingField_A.trialwiseMeanActivity_A(i) = nanmean(nanmean(squeeze(traces.traces_A(i,firingField_A.boundaries_AW(i,1):firingField_A.boundaries_AW(i,2),:)),1));
        firingField_A.trialwiseMeanActivity_X(i) = nanmean(nanmean(squeeze(traces.traces_X(i,firingField_A.boundaries_AW(i,1):firingField_A.boundaries_AW(i,2),:)),1));
        firingField_A.trialwiseSdActivity_A(i) = nanstd(nanmean(squeeze(traces.traces_A(i,firingField_A.boundaries_AW(i,1):firingField_A.boundaries_AW(i,2),:)),1));
        firingField_A.trialwiseSdActivity_X(i) = nanstd(nanmean(squeeze(traces.traces_X(i,firingField_A.boundaries_AW(i,1):firingField_A.boundaries_AW(i,2),:)),1));
    end
end
firingField_A.activationProbability_A = nanmean(firingField_A.activeTrials_A,2);
firingField_A.activationProbability_X = nanmean(firingField_A.activeTrials_X,2);
firingField_A.reliability_A = firingField_A.activationProbability_A > p.sca.reliabilityCriterion;
firingField_A.reliability_X = firingField_A.activationProbability_X > p.sca.reliabilityCriterion;
firingField_X.activeTrials_A = false(prop.numRois,length(trials.trials_A));
firingField_X.activeTrials_X = false(prop.numRois,length(trials.trials_X));
firingField_X.trialwiseMeanActivity_A = nan(prop.numRois,1);
firingField_X.trialwiseMeanActivity_X = nan(prop.numRois,1);
firingField_X.trialwiseSdActivity_A = nan(prop.numRois,1);
firingField_X.trialwiseSdActivity_X = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField_X.boundaries_AW(i,1)) && ~isnan(firingField_X.boundaries_AW(i,2))
        firingField_X.activeTrials_A(i,:) = nanmean(squeeze(traces.traces_A(i,firingField_X.boundaries_AW(i,1):firingField_X.boundaries_AW(i,2),:)),1) > firingField_X.activeInFiringFieldThreshold(i);
        firingField_X.activeTrials_X(i,:) = nanmean(squeeze(traces.traces_X(i,firingField_X.boundaries_AW(i,1):firingField_X.boundaries_AW(i,2),:)),1) > firingField_X.activeInFiringFieldThreshold(i);
        firingField_X.trialwiseMeanActivity_A(i) = nanmean(nanmean(squeeze(traces.traces_A(i,firingField_X.boundaries_AW(i,1):firingField_X.boundaries_AW(i,2),:)),1));
        firingField_X.trialwiseMeanActivity_X(i) = nanmean(nanmean(squeeze(traces.traces_X(i,firingField_X.boundaries_AW(i,1):firingField_X.boundaries_AW(i,2),:)),1));
        firingField_X.trialwiseSdActivity_A(i) = nanstd(nanmean(squeeze(traces.traces_A(i,firingField_X.boundaries_AW(i,1):firingField_X.boundaries_AW(i,2),:)),1));
        firingField_X.trialwiseSdActivity_X(i) = nanstd(nanmean(squeeze(traces.traces_X(i,firingField_X.boundaries_AW(i,1):firingField_X.boundaries_AW(i,2),:)),1));
    end
end
firingField_X.activationProbability_A = nanmean(firingField_X.activeTrials_A,2);
firingField_X.activationProbability_X = nanmean(firingField_X.activeTrials_X,2);
firingField_X.reliability_A = firingField_X.activationProbability_A > p.sca.reliabilityCriterion;
firingField_X.reliability_X = firingField_X.activationProbability_X > p.sca.reliabilityCriterion;

% mean activity in firing field (incl. bl-sub)
firingField_A.meanAmplitude_A = nan(prop.numRois,1);
firingField_A.meanAmplitude_X = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField_A.boundaries_AW(i,1)) && ~isnan(firingField_A.boundaries_AW(i,2))
        firingField_A.meanAmplitude_A(i) = nanmean(avgTraces.avgTraces_A(i,firingField_A.boundaries_AW(i,1):firingField_A.boundaries_AW(i,2)));
        firingField_A.meanAmplitude_X(i) = nanmean(avgTraces.avgTraces_X(i,firingField_A.boundaries_AW(i,1):firingField_A.boundaries_AW(i,2)));
    end
end
firingField_X.meanAmplitude_A = nan(prop.numRois,1);
firingField_X.meanAmplitude_X = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField_X.boundaries_AW(i,1)) && ~isnan(firingField_X.boundaries_AW(i,2))
        firingField_X.meanAmplitude_A(i) = nanmean(avgTraces.avgTraces_A(i,firingField_X.boundaries_AW(i,1):firingField_X.boundaries_AW(i,2)));
        firingField_X.meanAmplitude_X(i) = nanmean(avgTraces.avgTraces_X(i,firingField_X.boundaries_AW(i,1):firingField_X.boundaries_AW(i,2)));
    end
end
firingField_A.meanAmplitude_blSub_A = firingField_A.meanAmplitude_A - firingField_A.baseline;
firingField_A.meanAmplitude_blSub_X = firingField_A.meanAmplitude_X - firingField_A.baseline;
firingField_X.meanAmplitude_blSub_A = firingField_X.meanAmplitude_A - firingField_X.baseline;
firingField_X.meanAmplitude_blSub_X = firingField_X.meanAmplitude_X - firingField_X.baseline;

% mean bl-sub activity in firing field in active trials only
firingField_A.meanAmplitude_blSub_act_A = nan(prop.numRois,1);
firingField_A.meanAmplitude_blSub_act_X = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField_A.boundaries_AW(i,1)) && ~isnan(firingField_A.boundaries_AW(i,2))
        temp = squeeze(traces.traces_A(i,firingField_A.boundaries_AW(i,1):firingField_A.boundaries_AW(i,2),find(firingField_A.activeTrials_A(i,:))));
        firingField_A.meanAmplitude_blSub_act_A(i) = nanmean(nanmean(temp,1));
        temp = squeeze(traces.traces_X(i,firingField_A.boundaries_AW(i,1):firingField_A.boundaries_AW(i,2),find(firingField_A.activeTrials_X(i,:))));
        firingField_A.meanAmplitude_blSub_act_X(i) = nanmean(nanmean(temp,1));
    end
end
firingField_X.meanAmplitude_blSub_act_A = nan(prop.numRois,1);
firingField_X.meanAmplitude_blSub_act_X = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField_X.boundaries_AW(i,1)) && ~isnan(firingField_X.boundaries_AW(i,2))
        temp = squeeze(traces.traces_A(i,firingField_X.boundaries_AW(i,1):firingField_X.boundaries_AW(i,2),find(firingField_X.activeTrials_A(i,:))));
        firingField_X.meanAmplitude_blSub_act_A(i) = nanmean(nanmean(temp,1));
        temp = squeeze(traces.traces_X(i,firingField_X.boundaries_AW(i,1):firingField_X.boundaries_AW(i,2),find(firingField_X.activeTrials_X(i,:))));
        firingField_X.meanAmplitude_blSub_act_X(i) = nanmean(nanmean(temp,1));
    end
end

% mean norm activity in firing field
firingField_A.meanNormAmplitude_A = nan(prop.numRois,1);
firingField_A.meanNormAmplitude_X = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField_A.boundaries_AW(i,1)) && ~isnan(firingField_A.boundaries_AW(i,2))
        firingField_A.meanNormAmplitude_A(i) = nanmean(normAvgTraces.normAvgTraces_A(i,firingField_A.boundaries_AW(i,1):firingField_A.boundaries_AW(i,2)));
        firingField_A.meanNormAmplitude_X(i) = nanmean(normAvgTraces.normAvgTraces_X(i,firingField_A.boundaries_AW(i,1):firingField_A.boundaries_AW(i,2)));
    end
end
firingField_X.meanNormAmplitude_A = nan(prop.numRois,1);
firingField_X.meanNormAmplitude_X = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField_X.boundaries_AW(i,1)) && ~isnan(firingField_X.boundaries_AW(i,2))
        firingField_X.meanNormAmplitude_A(i) = nanmean(normAvgTraces.normAvgTraces_A(i,firingField_X.boundaries_AW(i,1):firingField_X.boundaries_AW(i,2)));
        firingField_X.meanNormAmplitude_X(i) = nanmean(normAvgTraces.normAvgTraces_X(i,firingField_X.boundaries_AW(i,1):firingField_X.boundaries_AW(i,2)));
    end
end

% selectivity index
firingField_A.selectivity_A = (firingField_A.meanNormAmplitude_A - firingField_A.meanNormAmplitude_X) ./ (firingField_A.meanNormAmplitude_A + firingField_A.meanNormAmplitude_X);
firingField_A.selectivity_X = (firingField_A.meanNormAmplitude_X - firingField_A.meanNormAmplitude_A) ./ (firingField_A.meanNormAmplitude_A + firingField_A.meanNormAmplitude_X);
firingField_X.selectivity_A = (firingField_X.meanNormAmplitude_A - firingField_X.meanNormAmplitude_X) ./ (firingField_X.meanNormAmplitude_A + firingField_X.meanNormAmplitude_X);
firingField_X.selectivity_X = (firingField_X.meanNormAmplitude_X - firingField_X.meanNormAmplitude_A) ./ (firingField_X.meanNormAmplitude_A + firingField_X.meanNormAmplitude_X);
temp = [firingField_A.selectivity_A;firingField_A.selectivity_X;firingField_X.selectivity_A;firingField_X.selectivity_X];
if any(temp > 1 | temp < -1)
    warning('Absolute selectivity index is sometimes greater than 1.')
end


%% Amplitude statistics (calculate p-values, fdr and q-values)

% A
firingField_A.diffBl_base = nan(prop.numRois,length(trials.trials_A));
firingField_A.diffBl_resp = nan(prop.numRois,length(trials.trials_A));
firingField_A.diffBl_diff = nan(prop.numRois,length(trials.trials_A));
firingField_A.diffBl_p = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField_A.boundaries_AW(i,1)) && ~isnan(firingField_A.boundaries_AW(i,2))
        firingField_A.diffBl_base(i,:) = squeeze(nanmean(traces.traces_A(i,prop.baselineWindow,:),2))';
        %these_base = squeeze(nanmax(traces.traces_A(i,prop.baselineWindow,:),[],2));
        % these_resp = squeeze(nanmean(traces.traces_A(i,firingField_A.boundaries_AW(i,1):firingField_A.boundaries_AW(i,2),:),2)); % mean in firing field
        firingField_A.diffBl_resp(i,:) = squeeze(nanmax(traces.traces_A(i,firingField_A.boundaries_AW(i,1):firingField_A.boundaries_AW(i,2),:),[],2))'; % max in firing field
        firingField_A.diffBl_diff(i,:) = firingField_A.diffBl_resp(i,:) - firingField_A.diffBl_base(i,:);
        firingField_A.diffBl_p(i) = signrank(firingField_A.diffBl_base(i,:)',firingField_A.diffBl_resp(i,:)');
    end
end
[firingField_A.diffBl_fdr,firingField_A.diffBl_q,firingField_A.diffBl_priori] = mafdr(firingField_A.diffBl_p); % double-checked that it works correctly with NaNs
if firingField_A.diffBl_priori>=1
    warning('The priori for at least one cluster is greater than or equal to 1.')
end
firingField_A.diffBl = firingField_A.diffBl_q < p.sca.qThreshold;
firingField_A.diffBl_pos = zeros(prop.numRois,1);
firingField_A.diffBl_pos(intersect(find(nanmean(firingField_A.diffBl_diff,2)>0),find(firingField_A.diffBl))) = 1;

% X
firingField_X.diffBl_base = nan(prop.numRois,length(trials.trials_X));
firingField_X.diffBl_resp = nan(prop.numRois,length(trials.trials_X));
firingField_X.diffBl_diff = nan(prop.numRois,length(trials.trials_X));
firingField_X.diffBl_p = nan(prop.numRois,1);
for i=1:prop.numRois
    if ~isnan(firingField_X.boundaries_AW(i,1)) && ~isnan(firingField_X.boundaries_AW(i,2))
        firingField_X.diffBl_base(i,:) = squeeze(nanmean(traces.traces_X(i,prop.baselineWindow,:),2))';
        % these_base = squeeze(nanmax(traces.traces_X(i,prop.baselineWindow,:),[],2));
        % these_resp = squeeze(nanmean(traces.traces_X(i,firingField_X.boundaries_AW(i,1):firingField_X.boundaries_AW(i,2),:),2));  % mean in firing field
        firingField_X.diffBl_resp(i,:) = squeeze(nanmax(traces.traces_X(i,firingField_X.boundaries_AW(i,1):firingField_X.boundaries_AW(i,2),:),[],2))'; % max in firing field
        firingField_X.diffBl_diff(i,:) = firingField_X.diffBl_resp(i,:) - firingField_X.diffBl_base(i,:);
        firingField_X.diffBl_p(i) = signrank(firingField_X.diffBl_base(i,:)',firingField_X.diffBl_resp(i,:)');
    end
end
[firingField_X.diffBl_fdr,firingField_X.diffBl_q,firingField_X.diffBl_priori] = mafdr(firingField_X.diffBl_p); % double-checked that it works correctly with NaNs
if firingField_X.diffBl_priori>=1
    warning('The priori for at least one cluster is greater than or equal to 1.')
end
firingField_X.diffBl = firingField_X.diffBl_q < p.sca.qThreshold;
firingField_X.diffBl_pos = zeros(prop.numRois,1);
firingField_X.diffBl_pos(intersect(find(nanmean(firingField_X.diffBl_diff,2)>0),find(firingField_X.diffBl))) = 1;


%% Apply seqeuence cell criteria

% apply sequence cell criteria
sca.passed.A = floor((sca.shuffling_A.significant + firingField_A.reliability_A + firingField_A.peakInAW + firingField_A.diffBl_pos)/4);
sca.passed.X = floor((sca.shuffling_X.significant + firingField_X.reliability_X + firingField_X.peakInAW + firingField_X.diffBl_pos)/4);

% sequence cell results
sca.passed.A_early = floor((sca.passed.A + firingField_A.early)/2);
sca.passed.A_late = floor((sca.passed.A + firingField_A.late)/2);
sca.passed.X_early = floor((sca.passed.X + firingField_X.early)/2);
sca.passed.X_late = floor((sca.passed.X + firingField_X.late)/2);
sca.passed.Aonly = floor((sca.passed.A + 1-sca.passed.X)/2);
sca.passed.Xonly = floor((1-sca.passed.A + sca.passed.X)/2);
sca.passed.AandX = floor((sca.passed.A + sca.passed.X)/2);
sca.passed.AorX = round((sca.passed.A + sca.passed.X)/2);
sca.passed.early = round((sca.passed.A_early + sca.passed.X_early)/2);
sca.passed.late = round((sca.passed.A_late + sca.passed.X_late)/2);
sca.passed = orderfields(sca.passed);

% sequence cell stats
sca.passed_stats.A_num = nansum(sca.passed.A);
sca.passed_stats.X_num = nansum(sca.passed.X);
sca.passed_stats.A_early_num = nansum(sca.passed.A_early);
sca.passed_stats.A_late_num = nansum(sca.passed.A_late);
sca.passed_stats.X_early_num = nansum(sca.passed.X_early);
sca.passed_stats.X_late_num = nansum(sca.passed.X_late);
sca.passed_stats.A_fractionOfCells = sca.passed_stats.A_num/nansum(iscell);
sca.passed_stats.X_fractionOfCells = sca.passed_stats.X_num/nansum(iscell);
sca.passed_stats.A_early_fractionOfCells = sca.passed_stats.A_early_num/nansum(iscell);
sca.passed_stats.A_late_fractionOfCells = sca.passed_stats.A_late_num/nansum(iscell);
sca.passed_stats.X_early_fractionOfCells = sca.passed_stats.X_early_num/nansum(iscell);
sca.passed_stats.X_late_fractionOfCells = sca.passed_stats.X_late_num/nansum(iscell);
sca.passed_stats.A_early_fractionOfACells = sca.passed_stats.A_early_num/sca.passed_stats.A_num;
sca.passed_stats.A_late_fractionOfACells = sca.passed_stats.A_late_num/sca.passed_stats.A_num;
sca.passed_stats.X_early_fractionOfXCells = sca.passed_stats.X_early_num/sca.passed_stats.X_num;
sca.passed_stats.X_late_fractionOfXCells = sca.passed_stats.X_late_num/sca.passed_stats.X_num;
sca.passed_stats.AorX_num = nansum(sca.passed.AorX);
sca.passed_stats.AandX_num = nansum(sca.passed.AandX);
sca.passed_stats.AorX_fractionOfCells = sca.passed_stats.AorX_num/nansum(iscell);
sca.passed_stats.AandX_fractionOfCells = sca.passed_stats.AandX_num/nansum(iscell);
sca.passed_stats.AandX_fractionOfAorXCells = sca.passed_stats.AandX_num/sca.passed_stats.AorX_num;
sca.passed_stats.early_num = nansum(sca.passed.early);
sca.passed_stats.late_num = nansum(sca.passed.late);
sca.passed_stats.early_fractionOfCells = sca.passed_stats.early_num/nansum(iscell);
sca.passed_stats.late_fractionOfCells = sca.passed_stats.late_num/nansum(iscell);
sca.passed_stats.early_fractionOfAorXCells = sca.passed_stats.early_num/sca.passed_stats.AorX_num;
sca.passed_stats.late_fractionOfAorXCells = sca.passed_stats.late_num/sca.passed_stats.AorX_num;
sca.passed_stats = orderfields(sca.passed_stats);


%% Second odour properties

% selectivity index for go vs no-go
sca.odour2.selectivity_AB_XB = (nanmean(normAvgTraces.normAvgTraces_AB(:,prop.odour2Window),2) - nanmean(normAvgTraces.normAvgTraces_XB(:,prop.odour2Window),2)) ./ (nanmean(normAvgTraces.normAvgTraces_AB(:,prop.odour2Window),2) + nanmean(normAvgTraces.normAvgTraces_XB(:,prop.odour2Window),2));
sca.odour2.selectivity_XY_AY = (nanmean(normAvgTraces.normAvgTraces_XY(:,prop.odour2Window),2) - nanmean(normAvgTraces.normAvgTraces_AY(:,prop.odour2Window),2)) ./ (nanmean(normAvgTraces.normAvgTraces_XY(:,prop.odour2Window),2) + nanmean(normAvgTraces.normAvgTraces_AY(:,prop.odour2Window),2));
temp = [sca.odour2.selectivity_AB_XB;sca.odour2.selectivity_XY_AY];
if any(temp > 1 | temp < -1)
    warning('Absolute second odour selectivity index is sometimes greater than 1.')
end


%% Sort output

sca.trials = orderfields(trials);
sca.traces = orderfields(traces);
sca.avgTraces = orderfields(avgTraces);
sca.normAvgTraces = orderfields(normAvgTraces);
sca.firingField_A = orderfields(firingField_A);
sca.firingField_X = orderfields(firingField_X);
prop.iscells = iscell;
sca.prop = orderfields(prop);
sca.p = p.sca;
sca = orderfields(sca);

end