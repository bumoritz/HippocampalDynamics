function [tngn,nft_binned] = tuningAnalysisDff(info,iscell,ops,p,path,task,s2p_meta,sync_beh,act)
% p = get_p; sync_beh = paq_beh; 

p.tngn.zscore = false


%% Preparations

disp('--- Preparations')

% start parallel pool (if not started yet)
try
    parpool; % can be closed with: delete(gcp('nocreate'));
catch
end

% pre-process activity measure
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.tngn,p,sync_beh,iscell);
bl_trial_frames_binned = [prop.trial_frames_binned(1,:)-[1:5]';prop.trial_frames_binned(1:10,:)];
nft_bl_binned = reshape(nf_binned(:,bl_trial_frames_binned),[],size(bl_trial_frames_binned,1),size(bl_trial_frames_binned,2));
nft_bl_avg_binned = nanmean(nft_bl_binned,3);

% create trials, traces, avgTraces and normAvgTraces structs
if ops.tngn.do_allTrials
    trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
    [traces_all,avgTraces_all,normAvgTraces_all] = createTracesStructs(nft_binned,trials_all);
end
if ops.tngn.do_60t || ops.tngn.do_60t_onlyFirst
    trials_60t = {}; traces_60t = {}; avgTraces_60t = {}; normAvgTraces_60t = {};
    if ~ops.tngn.do_60t
        temp = 1;
    else
        temp = floor(prop.numTrials/60);
    end
    for i=1:temp
        trials_60t{i} = createTrialsStruct(task,(i-1)*60+1:i*60);
        [traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i}] = createTracesStructs(nft_binned,trials_60t{i});
    end
end
if ops.tngn.do_100t || ops.tngn.do_100t_onlyFirst
    trials_100t = {}; traces_100t = {}; avgTraces_100t = {}; normAvgTraces_100t = {};
    if ~ops.tngn.do_100t
        temp = 1;
    else
        temp = floor(prop.numTrials/100);
    end
    for i=1:temp
        trials_100t{i} = createTrialsStruct(task,(i-1)*100+1:i*100);
        [traces_100t{i},avgTraces_100t{i},normAvgTraces_100t{i}] = createTracesStructs(nft_binned,trials_100t{i});
    end
end

% stimVersion
if ops.tngn.do_allTrials_stimVersion
    trials_all_stimVersion = createTrialsStruct(task,1:prop.numTrials_incl,true);
    [traces_all_stimVersion,avgTraces_all_stimVersion,normAvgTraces_all_stimVersion] = createTracesStructs(nft_binned,trials_all_stimVersion);
end
if ops.tngn.do_100t_stimVersion
    trials_100t_stimVersion = {}; traces_100t_stimVersion = {}; avgTraces_100t_stimVersion = {}; normAvgTraces_100t_stimVersion = {};
    for i=1:floor(prop.numTrials/100)
        trials_100t_stimVersion{i} = createTrialsStruct(task,(i-1)*100+1:i*100,true);
        [traces_100t_stimVersion{i},avgTraces_100t_stimVersion{i},normAvgTraces_100t_stimVersion{i}] = createTracesStructs(nft_binned,trials_100t_stimVersion{i});
    end
end


%% Shuffling

disp('--- Shuffling')

if ops.tngn.do_allTrials
    if p.tngn.analysisBlockRestrictedNf
        this_binRange = [prop.sync_all_binned(1)-length(p.general.bins_pre),prop.sync_all_binned(end)+length(p.general.t_binned)-length(p.general.bins_pre)];
    else
        this_binRange = [1,size(nf_binned,2)];
    end
    shuffling_all.A_AW = shufflingCriterion(prop,nf_binned,trials_all.stimuli.A,traces_all.A,p.general.bins_analysisWindow,p,this_binRange);
    shuffling_all.X_AW = shufflingCriterion(prop,nf_binned,trials_all.stimuli.X,traces_all.X,p.general.bins_analysisWindow,p,this_binRange);
    if ops.tngn.do_eventWiseAnalysis
        shuffling_all.A_odour1 = shufflingCriterion(prop,nf_binned,trials_all.stimuli.A,traces_all.A,p.general.bins_odour1Window,p,this_binRange);
        shuffling_all.X_odour1 = shufflingCriterion(prop,nf_binned,trials_all.stimuli.X,traces_all.X,p.general.bins_odour1Window,p,this_binRange);
        shuffling_all.B_odour2 = shufflingCriterion(prop,nf_binned,trials_all.stimuli.B,traces_all.B,p.general.bins_odour2Window,p,this_binRange);
        shuffling_all.Y_odour2 = shufflingCriterion(prop,nf_binned,trials_all.stimuli.Y,traces_all.Y,p.general.bins_odour2Window,p,this_binRange);
        shuffling_all.rew_reward = shufflingCriterion(prop,nf_binned,trials_all.outcome.rew,traces_all.rew,p.general.bins_rewardWindow,p,this_binRange);
        shuffling_all.norew_reward = shufflingCriterion(prop,nf_binned,trials_all.outcome.norew,traces_all.norew,p.general.bins_rewardWindow,p,this_binRange);
    end
end
if ops.tngn.do_60t || ops.tngn.do_60t_onlyFirst
    shuffling_60t = {};
    if ~ops.tngn.do_60t
        temp = 1;
    else
        temp = floor(prop.numTrials/60);
    end
    for i=1:temp
        if p.tngn.analysisBlockRestrictedNf
            this_binRange = [prop.sync_all_binned((i-1)*60+1)-length(p.general.bins_pre),prop.sync_all_binned(i*60)+length(p.general.t_binned)-length(p.general.bins_pre)];
        else
            this_binRange = [1,size(nf_binned,2)];
        end
        shuffling_60t{i}.A_AW = shufflingCriterion(prop,nf_binned,trials_60t{i}.stimuli.A,traces_60t{i}.A,p.general.bins_analysisWindow,p,this_binRange);
        shuffling_60t{i}.X_AW = shufflingCriterion(prop,nf_binned,trials_60t{i}.stimuli.X,traces_60t{i}.X,p.general.bins_analysisWindow,p,this_binRange);
        if ops.tngn.do_eventWiseAnalysis
            shuffling_60t{i}.A_odour1 = shufflingCriterion(prop,nf_binned,trials_60t{i}.stimuli.A,traces_60t{i}.A,p.general.bins_odour1Window,p,this_binRange);
            shuffling_60t{i}.X_odour1 = shufflingCriterion(prop,nf_binned,trials_60t{i}.stimuli.X,traces_60t{i}.X,p.general.bins_odour1Window,p,this_binRange);
            shuffling_60t{i}.B_odour2 = shufflingCriterion(prop,nf_binned,trials_60t{i}.stimuli.B,traces_60t{i}.B,p.general.bins_odour2Window,p,this_binRange);
            shuffling_60t{i}.Y_odour2 = shufflingCriterion(prop,nf_binned,trials_60t{i}.stimuli.Y,traces_60t{i}.Y,p.general.bins_odour2Window,p,this_binRange);
            shuffling_60t{i}.rew_reward = shufflingCriterion(prop,nf_binned,trials_60t{i}.outcome.rew,traces_60t{i}.rew,p.general.bins_rewardWindow,p,this_binRange);
            shuffling_60t{i}.norew_reward = shufflingCriterion(prop,nf_binned,trials_60t{i}.outcome.norew,traces_60t{i}.norew,p.general.bins_rewardWindow,p,this_binRange);
        end
    end
end
if ops.tngn.do_100t || ops.tngn.do_100t_onlyFirst
    shuffling_100t = {};
    if ~ops.tngn.do_100t
        temp = 1;
    else
        temp = floor(prop.numTrials/100);
    end
    for i=1:temp
        if p.tngn.analysisBlockRestrictedNf
            try
                this_binRange = [prop.sync_all_binned((i-1)*100+1)-length(p.general.bins_pre),prop.sync_all_binned(i*100)+length(p.general.t_binned)-length(p.general.bins_pre)];
            catch
                this_binRange = [prop.sync_all_binned((i-1)*100+1)-length(p.general.bins_pre),prop.sync_all_binned(end-2)+length(p.general.t_binned)-length(p.general.bins_pre)];
            end
        else
            this_binRange = [1,size(nf_binned,2)];
        end
        shuffling_100t{i}.A_AW = shufflingCriterion(prop,nf_binned,trials_100t{i}.stimuli.A,traces_100t{i}.A,p.general.bins_analysisWindow,p,this_binRange); % HOW TO DEAL WITH nf_binned trial-wise?
        shuffling_100t{i}.X_AW = shufflingCriterion(prop,nf_binned,trials_100t{i}.stimuli.X,traces_100t{i}.X,p.general.bins_analysisWindow,p,this_binRange);
        if ops.tngn.do_eventWiseAnalysis
            shuffling_100t{i}.A_odour1 = shufflingCriterion(prop,nf_binned,trials_100t{i}.stimuli.A,traces_100t{i}.A,p.general.bins_odour1Window,p,this_binRange);
            shuffling_100t{i}.X_odour1 = shufflingCriterion(prop,nf_binned,trials_100t{i}.stimuli.X,traces_100t{i}.X,p.general.bins_odour1Window,p,this_binRange);
            shuffling_100t{i}.B_odour2 = shufflingCriterion(prop,nf_binned,trials_100t{i}.stimuli.B,traces_100t{i}.B,p.general.bins_odour2Window,p,this_binRange);
            shuffling_100t{i}.Y_odour2 = shufflingCriterion(prop,nf_binned,trials_100t{i}.stimuli.Y,traces_100t{i}.Y,p.general.bins_odour2Window,p,this_binRange);
            shuffling_100t{i}.rew_reward = shufflingCriterion(prop,nf_binned,trials_100t{i}.outcome.rew,traces_100t{i}.rew,p.general.bins_rewardWindow,p,this_binRange);
            shuffling_100t{i}.norew_reward = shufflingCriterion(prop,nf_binned,trials_100t{i}.outcome.norew,traces_100t{i}.norew,p.general.bins_rewardWindow,p,this_binRange);
        end
    end
end

% stimVersion
if ops.tngn.do_allTrials_stimVersion
    if p.tngn.analysisBlockRestrictedNf
        this_binRange = [prop.sync_all_binned(1)-length(p.general.bins_pre),prop.sync_all_binned(end)+length(p.general.t_binned)-length(p.general.bins_pre)];
    else
        this_binRange = [1,size(nf_binned,2)];
    end
    shuffling_all_stimVersion.Acatch_AW = shufflingCriterion(prop,nf_binned,trials_all_stimVersion.stimuli_var0.A,traces_all_stimVersion.A_var0,p.general.bins_analysisWindow,p,this_binRange);
    shuffling_all_stimVersion.Astim_AW = shufflingCriterion(prop,nf_binned,trials_all_stimVersion.stimuli_var1.A,traces_all_stimVersion.A_var1,p.general.bins_analysisWindow,p,this_binRange);
    shuffling_all_stimVersion.Xcatch_AW = shufflingCriterion(prop,nf_binned,trials_all_stimVersion.stimuli_var0.X,traces_all_stimVersion.X_var0,p.general.bins_analysisWindow,p,this_binRange);
    shuffling_all_stimVersion.Xstim_AW = shufflingCriterion(prop,nf_binned,trials_all_stimVersion.stimuli_var1.X,traces_all_stimVersion.X_var1,p.general.bins_analysisWindow,p,this_binRange);
end
if ops.tngn.do_100t_stimVersion
    shuffling_100t_stimVersion = {};
    for i=1:floor(prop.numTrials/100)
        if p.tngn.analysisBlockRestrictedNf
            try
                this_binRange = [prop.sync_all_binned((i-1)*100+1)-length(p.general.bins_pre),prop.sync_all_binned(i*100)+length(p.general.t_binned)-length(p.general.bins_pre)];
            catch
                this_binRange = [prop.sync_all_binned((i-1)*100+1)-length(p.general.bins_pre),prop.sync_all_binned(end-2)+length(p.general.t_binned)-length(p.general.bins_pre)];
            end
        else
            this_binRange = [1,size(nf_binned,2)];
        end
        shuffling_100t_stimVersion{i}.Acatch_AW = shufflingCriterion(prop,nf_binned,trials_100t_stimVersion{i}.stimuli_var0.A,traces_100t_stimVersion{i}.A_var0,p.general.bins_analysisWindow,p,this_binRange);
        shuffling_100t_stimVersion{i}.Astim_AW = shufflingCriterion(prop,nf_binned,trials_100t_stimVersion{i}.stimuli_var1.A,traces_100t_stimVersion{i}.A_var1,p.general.bins_analysisWindow,p,this_binRange);
        shuffling_100t_stimVersion{i}.Xcatch_AW = shufflingCriterion(prop,nf_binned,trials_100t_stimVersion{i}.stimuli_var0.X,traces_100t_stimVersion{i}.X_var0,p.general.bins_analysisWindow,p,this_binRange);        
        shuffling_100t_stimVersion{i}.Xstim_AW = shufflingCriterion(prop,nf_binned,trials_100t_stimVersion{i}.stimuli_var1.X,traces_100t_stimVersion{i}.X_var1,p.general.bins_analysisWindow,p,this_binRange);        
    end
end


%% Firing field properties

disp('--- Calculating firing field properties')

if ops.tngn.do_allTrials
    firingField_all.A_AW = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'AW','A','X',prop,p);
    firingField_all.X_AW = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'AW','X','A',prop,p);
    if ops.tngn.do_eventWiseAnalysis
        firingField_all.A_odour1 = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'odour1','A','X',prop,p);
        firingField_all.X_odour1 = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'odour1','X','A',prop,p);
        firingField_all.B_odour2 = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'odour2','B','Y',prop,p);
        firingField_all.Y_odour2 = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'odour2','Y','B',prop,p);
        firingField_all.rew_reward = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'reward','rew','norew',prop,p);
        firingField_all.norew_reward = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'reward','norew','rew',prop,p);
    end
end
if ops.tngn.do_60t || ops.tngn.do_60t_onlyFirst
    firingField_60t = {};
    if ~ops.tngn.do_60t
        temp = 1;
    else
        temp = floor(prop.numTrials/60);
    end
    for i=1:temp
        firingField_60t{i}.A_AW = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'AW','A','X',prop,p);
        firingField_60t{i}.X_AW = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'AW','X','A',prop,p);
        if ops.tngn.do_eventWiseAnalysis
            firingField_60t{i}.A_odour1 = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'odour1','A','X',prop,p);
            firingField_60t{i}.X_odour1 = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'odour1','X','A',prop,p);
            firingField_60t{i}.B_odour2 = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'odour2','B','Y',prop,p);
            firingField_60t{i}.Y_odour2 = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'odour2','Y','B',prop,p);
            firingField_60t{i}.rew_reward = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'reward','rew','norew',prop,p);
            firingField_60t{i}.norew_reward = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'reward','norew','rew',prop,p);
        end
    end
end
if ops.tngn.do_100t || ops.tngn.do_100t_onlyFirst
    firingField_100t = {};
    if ~ops.tngn.do_100t
        temp = 1;
    else
        temp = floor(prop.numTrials/100);
    end
    for i=1:temp
        firingField_100t{i}.A_AW = firingFieldProperties(traces_100t{i},avgTraces_100t{i},normAvgTraces_100t{i},'AW','A','X',prop,p);
        firingField_100t{i}.X_AW = firingFieldProperties(traces_100t{i},avgTraces_100t{i},normAvgTraces_100t{i},'AW','X','A',prop,p);
        if ops.tngn.do_eventWiseAnalysis
            firingField_100t{i}.A_odour1 = firingFieldProperties(traces_100t{i},avgTraces_100t{i},normAvgTraces_100t{i},'odour1','A','X',prop,p);
            firingField_100t{i}.X_odour1 = firingFieldProperties(traces_100t{i},avgTraces_100t{i},normAvgTraces_100t{i},'odour1','X','A',prop,p);
            firingField_100t{i}.B_odour2 = firingFieldProperties(traces_100t{i},avgTraces_100t{i},normAvgTraces_100t{i},'odour2','B','Y',prop,p);
            firingField_100t{i}.Y_odour2 = firingFieldProperties(traces_100t{i},avgTraces_100t{i},normAvgTraces_100t{i},'odour2','Y','B',prop,p);
            firingField_100t{i}.rew_reward = firingFieldProperties(traces_100t{i},avgTraces_100t{i},normAvgTraces_100t{i},'reward','rew','norew',prop,p);
            firingField_100t{i}.norew_reward = firingFieldProperties(traces_100t{i},avgTraces_100t{i},normAvgTraces_100t{i},'reward','norew','rew',prop,p);
        end
    end
end

% stimVersion
if ops.tngn.do_allTrials_stimVersion
    firingField_all_stimVersion.Acatch_AW = firingFieldProperties(traces_all_stimVersion,avgTraces_all_stimVersion,normAvgTraces_all_stimVersion,'AW','A_var0','X_var0',prop,p);
    firingField_all_stimVersion.Astim_AW = firingFieldProperties(traces_all_stimVersion,avgTraces_all_stimVersion,normAvgTraces_all_stimVersion,'AW','A_var1','X_var1',prop,p);
    firingField_all_stimVersion.Xcatch_AW = firingFieldProperties(traces_all_stimVersion,avgTraces_all_stimVersion,normAvgTraces_all_stimVersion,'AW','X_var0','A_var0',prop,p);
    firingField_all_stimVersion.Xstim_AW = firingFieldProperties(traces_all_stimVersion,avgTraces_all_stimVersion,normAvgTraces_all_stimVersion,'AW','X_var1','A_var1',prop,p);
end
if ops.tngn.do_100t_stimVersion
    firingField_100t_stimVersion = {};
    for i=1:floor(prop.numTrials/100)
        firingField_100t_stimVersion{i}.Acatch_AW = firingFieldProperties(traces_100t_stimVersion{i},avgTraces_100t_stimVersion{i},normAvgTraces_100t_stimVersion{i},'AW','A_var0','X_var0',prop,p);
        firingField_100t_stimVersion{i}.Astim_AW = firingFieldProperties(traces_100t_stimVersion{i},avgTraces_100t_stimVersion{i},normAvgTraces_100t_stimVersion{i},'AW','A_var1','X_var1',prop,p);
        firingField_100t_stimVersion{i}.Xcatch_AW = firingFieldProperties(traces_100t_stimVersion{i},avgTraces_100t_stimVersion{i},normAvgTraces_100t_stimVersion{i},'AW','X_var0','A_var0',prop,p);
        firingField_100t_stimVersion{i}.Xstim_AW = firingFieldProperties(traces_100t_stimVersion{i},avgTraces_100t_stimVersion{i},normAvgTraces_100t_stimVersion{i},'AW','X_var1','A_var1',prop,p);
    end
end


%% Inhibition/artefact evasion

% if ops.tngn.do_allTrials
%     firingField_all.A_1s.negativeGoing = double(nanmean(avgTraces_all.A(:,p.general.bins_1st1s),2) - nanmean(avgTraces_all.A(:,p.general.bins_base1s),2) <0);
%     firingField_all.A_1s.negativeGoing(~prop.iscell==1) = NaN;
%     firingField_all.X_1s.negativeGoing = double(nanmean(avgTraces_all.X(:,p.general.bins_1st1s),2) - nanmean(avgTraces_all.X(:,p.general.bins_base1s),2) <0);
%     firingField_all.X_1s.negativeGoing(~prop.iscell==1) = NaN;
% end
% if ops.tngn.do_100t || ops.tngn.do_100t_onlyFirst
%     if ~ops.tngn.do_100t
%         temp = 1;
%     else
%         temp = floor(prop.numTrials/100);
%     end
%     for i=1:temp
%         firingField_100t{i}.A_1s.negativeGoing = double(nanmean(avgTraces_100t{i}.A(:,p.general.bins_1st1s),2) - nanmean(avgTraces_100t{i}.A(:,p.general.bins_base1s),2) <0);
%         firingField_100t{i}.A_1s.negativeGoing(~prop.iscell==1) = NaN;
%         firingField_100t{i}.X_1s.negativeGoing = double(nanmean(avgTraces_100t{i}.X(:,p.general.bins_1st1s),2) - nanmean(avgTraces_100t{i}.X(:,p.general.bins_base1s),2) <0);
%         firingField_100t{i}.X_1s.negativeGoing(~prop.iscell==1) = NaN;
%     end
% end

%% Applying tuning criteria

disp('--- Applying tuning criteria')

if ops.tngn.do_allTrials
    [passed_all.AW,passed_stats_all.AW] = tuningCriteria(shuffling_all,firingField_all,'AW','A','X',p,prop);
%     [passed_nneg_all.AW,passed_stats_nneg_all.AW] = tuningCriteria_nneg(shuffling_all,firingField_all,'AW','A','X',p,prop);
    if ops.tngn.do_eventWiseAnalysis
        [passed_all.odour1,passed_stats_all.odour1] = tuningCriteria(shuffling_all,firingField_all,'odour1','A','X',p,prop);
        [passed_all.odour2,passed_stats_all.odour2] = tuningCriteria(shuffling_all,firingField_all,'odour2','B','Y',p,prop);
        [passed_all.reward,passed_stats_all.reward] = tuningCriteria(shuffling_all,firingField_all,'reward','rew','norew',p,prop);
    end
end
if ops.tngn.do_60t || ops.tngn.do_60t_onlyFirst
    passed_60t = {}; passed_stats_60t = {};
    if ~ops.tngn.do_60t
        temp = 1;
    else
        temp = floor(prop.numTrials/60);
    end
    for i=1:temp
        [passed_60t{i}.AW,passed_stats_60t{i}.AW] = tuningCriteria(shuffling_60t{i},firingField_60t{i},'AW','A','X',p,prop);
        [passed_nneg_60t{i}.AW,passed_stats_nneg_60t{i}.AW] = tuningCriteria_nneg(shuffling_60t{i},firingField_60t{i},'AW','A','X',p,prop);
        if ops.tngn.do_eventWiseAnalysis
            [passed_60t{i}.odour1,passed_stats_60t{i}.odour1] = tuningCriteria(shuffling_60t{i},firingField_60t{i},'odour1','A','X',p,prop);
            [passed_60t{i}.odour2,passed_stats_60t{i}.odour2] = tuningCriteria(shuffling_60t{i},firingField_60t{i},'odour2','B','Y',p,prop);
            [passed_60t{i}.reward,passed_stats_60t{i}.reward] = tuningCriteria(shuffling_60t{i},firingField_60t{i},'reward','rew','norew',p,prop);
        end
    end
end
if ops.tngn.do_100t || ops.tngn.do_100t_onlyFirst
    passed_100t = {}; passed_stats_100t = {};
    if ~ops.tngn.do_100t
        temp = 1;
    else
        temp = floor(prop.numTrials/100);
    end
    for i=1:temp
        [passed_100t{i}.AW,passed_stats_100t{i}.AW] = tuningCriteria(shuffling_100t{i},firingField_100t{i},'AW','A','X',p,prop);
%         [passed_nneg_100t{i}.AW,passed_stats_nneg_100t{i}.AW] = tuningCriteria_nneg(shuffling_100t{i},firingField_100t{i},'AW','A','X',p,prop);
        if ops.tngn.do_eventWiseAnalysis
            [passed_100t{i}.odour1,passed_stats_100t{i}.odour1] = tuningCriteria(shuffling_100t{i},firingField_100t{i},'odour1','A','X',p,prop);
            [passed_100t{i}.odour2,passed_stats_100t{i}.odour2] = tuningCriteria(shuffling_100t{i},firingField_100t{i},'odour2','B','Y',p,prop);
            [passed_100t{i}.reward,passed_stats_100t{i}.reward] = tuningCriteria(shuffling_100t{i},firingField_100t{i},'reward','rew','norew',p,prop);
        end
    end
end

% stimVersion
if ops.tngn.do_allTrials_stimVersion
    [passed_all_stimVersion_catch.AW,passed_stats_all_stimVersion_catch.AW] = tuningCriteria(shuffling_all_stimVersion,firingField_all_stimVersion,'AW','Acatch','Xcatch',p,prop);
    [passed_all_stimVersion_stim.AW,passed_stats_all_stimVersion_stim.AW] = tuningCriteria(shuffling_all_stimVersion,firingField_all_stimVersion,'AW','Astim','Xstim',p,prop);
end
if ops.tngn.do_100t_stimVersion
    passed_100t_stimVersion_catch = {}; passed_stats_100t_stimVersion_catch = {};
    passed_100t_stimVersion_stim = {}; passed_stats_100t_stimVersion_stim = {};
    for i=1:floor(prop.numTrials/100)
        [passed_100t_stimVersion_catch{i}.AW,passed_stats_100t_stimVersion_catch{i}.AW] = tuningCriteria(shuffling_100t_stimVersion{i},firingField_100t_stimVersion{i},'AW','Acatch','Xcatch',p,prop);
        [passed_100t_stimVersion_stim{i}.AW,passed_stats_100t_stimVersion_stim{i}.AW] = tuningCriteria(shuffling_100t_stimVersion{i},firingField_100t_stimVersion{i},'AW','Astim','Xstim',p,prop);
    end
end


%% Snake plots

% preparations
if ~exist([path.filepart_outX,'plots/SnakePlots/'],'dir')
    mkdir([path.filepart_outX,'plots/SnakePlots/']);
end
% these_files = dir(fullfile([path.filepart_outX,'plots/SnakePlots/'],'*'));
% for k=3:length(these_files)
%   delete(fullfile([path.filepart_outX,'plots/SnakePlots/'],these_files(k).name));
% end
if strcmp(p.tngn.activityMeasure,"casc_50c_beh")
    this_activity_measure = 'casc_50c';
elseif strcmp(p.tngn.activityMeasure,"casc_50g_beh")
    this_activity_measure = 'casc_50g';
elseif strcmp(p.tngn.activityMeasure,"casc_100c_beh")
    this_activity_measure = 'casc_100c';
elseif strcmp(p.tngn.activityMeasure,"casc_100g_beh")
    this_activity_measure = 'casc_100g';
elseif strcmp(p.tngn.activityMeasure,"casc_200g_beh")
    this_activity_measure = 'casc_200g';
elseif strcmp(p.tngn.activityMeasure,"spks_beh")
    this_activity_measure = 'spks';
elseif strcmp(p.tngn.activityMeasure,"dFF_beh")
    this_activity_measure = 'dFF';
elseif strcmp(p.tngn.activityMeasure,"spksn_beh")
    this_activity_measure = 'spksn';
elseif strcmp(p.tngn.activityMeasure,"spksnn_beh")
    this_activity_measure = 'spksnn';
elseif strcmp(p.tngn.activityMeasure,"dFFn_beh")
    this_activity_measure = 'dFFn';
elseif strcmp(p.tngn.activityMeasure,"dFFnn_beh")
    this_activity_measure = 'dFFnn';
elseif strcmp(p.tngn.activityMeasure,"F_beh")
    this_activity_measure = 'F';
end

if ops.tngn.do_allTrials

%     % AX, all trials, 8tile, oddeven, 13
%     plt = struct(); plt.split = 'odd_even';
%     F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[1;3],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=AX, all trials'])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_oddeven_13.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_oddeven_13.png']);
%     

    % AXonly, all trials, 8tile, oddeven, 13
%     this_lower_A = nanmean(nft_bl_avg_binned(find(passed_all.AW.Aonly==1),:),2);
%     this_lower_X = nanmean(nft_bl_avg_binned(find(passed_all.AW.Xonly==1),:),2);
    this_lower = nanmean(nft_bl_avg_binned(:,2));
    plt = struct(); plt.split = 'odd_even'; plt.lower = this_lower; %{this_lower_A,this_lower_X};
    F = snakePlots(traces_all,{find(passed_all.AW.Aonly==1),find(passed_all.AW.Xonly==1)},[],[1;3],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=AXonly, all trials'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AXonly_allTrials_8tile_oddeven_13.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AXonly_allTrials_8tile_oddeven_13.png']);
     
%    % AX, all trials, 8tile, oddeven, 13, two-sided
%     plt = struct(); plt.split = 'odd_even'; plt.normalisationType = 'two-sided'; plt.normalisationWindow = 'all';
%     F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[1;3],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=AX, all trials'])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_twosided_oddeven_13.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_twosided_oddeven_13.png']);
% 
%     % AXonly, all trials, 8tile, oddeven, 13, two-sided
%     plt = struct(); plt.split = 'odd_even'; plt.normalisationType = 'two-sided'; plt.normalisationWindow = 'all';
%     F = snakePlots(traces_all,{find(passed_all.AW.Aonly==1),find(passed_all.AW.Xonly==1)},[],[1;3],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=AXonly, all trials'])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AXonly_allTrials_8tile_twosided_oddeven_13.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AXonly_allTrials_8tile_twosided_oddeven_13.png']);

    if ~ops.tngn.skip_boringSnakePlots
        
        % AX, all trials, 12tile, corrOcorrMincorr, same, 14
        plt = struct(); plt.split = 'corrO_corrM_incorr'; plt.split3_numTrials = 'same';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[1;4],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_same_14.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_same_14.png']);

        % AX, all trials, 4tile, 12
        plt = struct(); plt.split = 'none';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[1;2],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_4tile_12.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_4tile_12.png']);

        % AX, all trials, 4tile, 21
        plt = struct(); plt.split = 'none';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[2;1],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_4tile_21.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_4tile_21.png']);

        % AX, all trials, 8tile, oddeven, 31
        plt = struct(); plt.split = 'odd_even';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[3;1],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_oddeven_31.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_oddeven_31.png']);

        % AX, all trials, 8tile, corrincorr, 13
        plt = struct(); plt.split = 'corr_incorr';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[1;3],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_corrincorr_13.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_corrincorr_13.png']);

        % AX, all trials, 8tile, corrincorr, 24
        plt = struct(); plt.split = 'corr_incorr';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[2;4],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_corrincorr_24.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_corrincorr_24.png']);

        % AX, all trials, 8tile, corrincorr, 31
        plt = struct(); plt.split = 'corr_incorr';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[3;1],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_corrincorr_31.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_corrincorr_31.png']);
    
        % AX, all trials, 12tile, corrOcorrMincorr, max, 14
        plt = struct(); plt.split = 'corrO_corrM_incorr'; plt.split3_numTrials = 'max';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[1;4],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_max_14.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_max_14.png']);

        % AX, all trials, 12tile, corrOcorrMincorr, max, 41
        plt = struct(); plt.split = 'corrO_corrM_incorr'; plt.split3_numTrials = 'max';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[4;1],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_max_41.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_max_41.png']);

        % AX, all trials, 12tile, corrOcorrMincorr, same, 41
        plt = struct(); plt.split = 'corrO_corrM_incorr'; plt.split3_numTrials = 'same';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[4;1],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_same_41.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_same_41.png']);
    end
end


if ops.tngn.do_60t_onlyFirst
    % AX, f60t, 8tile, oddeven, 13
    plt = struct(); plt.split = 'odd_even';
    F = snakePlots(traces_60t{1},{find(passed_60t{1}.AW.A==1),find(passed_60t{1}.AW.X==1)},[],[1;3],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=AX, first 60t'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f60t_8tile_oddeven_13.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f60t_8tile_oddeven_13.png']);
    
    % AX, f60t, 8tile, oddeven, 13
    plt = struct(); plt.split = 'odd_even';
    F = snakePlots(traces_60t{1},{find(passed_60t{1}.AW.A==1),find(passed_60t{1}.AW.X==1)},[],[1;3],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=AX, first 60t'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f60t_8tile_oddeven_13.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f60t_8tile_oddeven_13.png']);

    % AX, f60t, 12tile, corrOcorrMincorr, same, 14
    plt = struct(); plt.split = 'corrO_corrM_incorr'; plt.split3_numTrials = 'same';
    F = snakePlots(traces_60t{1},{find(passed_60t{1}.AW.A==1),find(passed_60t{1}.AW.X==1)},[],[1;4],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=AX, first 60t'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f60t_12tile_corrOcorrMincorr_same_14.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f60t_12tile_corrOcorrMincorr_same_14.png']);
end

if ops.tngn.do_100t
    
    % Aonly, 100t, ktile
    these_idcs = {};
    for i=1:length(passed_100t)
        these_idcs{i} = find(passed_100t{i}.AW.Aonly==1);
    end
    plt = struct(); plt.split = 'none'; plt.normalisationRef = 'plotwise'; plt.sequences = "A"; plt.trialTypes = "A";    
    [F,these_normAvgTraces_A,these_idcs_sorted_A,plt] = snakePlots(traces_100t,these_idcs,[],[1:length(passed_100t)],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=Aonly, 100t, block-wise'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Aonly_100t_ktile.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Aonly_100t_ktile.png']);
    
    % Xonly, 100t, ktile
    these_idcs = {};
    for i=1:length(passed_100t)
        these_idcs{i} = find(passed_100t{i}.AW.Xonly==1);
    end
    plt = struct(); plt.split = 'none'; plt.normalisationRef = 'plotwise'; plt.sequences = "X"; plt.trialTypes = "X";    
    [F,these_normAvgTraces_X,these_idcs_sorted_X,plt] = snakePlots(traces_100t,these_idcs,[],[1:length(passed_100t)],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=Xonly, 100t, block-wise'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Xonly_100t_ktile.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Xonly_100t_ktile.png']);
    
    %Aonlx and Xonly, ktile
    F = default_figure([-20,0.5,20,5]); nrows = 1; ncols = length(passed_100t); numBlocks = length(passed_100t); these_refs_sort = [1:length(passed_100t)];
    for r=1:nrows
        for c=1:ncols
            subplot(nrows,ncols,(r-1)*ncols+c)
            this_data_A = these_normAvgTraces_A{r,c}(these_idcs_sorted_A{c},:);
            this_data_X = these_normAvgTraces_X{r,c}(these_idcs_sorted_X{c},:);
            this_data = [this_data_A;this_data_X];
            
            
            [~,temp] = nanmax(this_data(:,p.general.bins_analysisWindow),[],2);
            [~,temp2] = sort(temp);

            heatMap_task(this_data(temp2,:),NaN,[],p,info,plt);
        end
    end
        n = 0;
        temp = split(plt.trialTypes,'_');
        temp1 = split(plt.sequences,'_');
        temp2 = split(plt.split,'_');
        for r=1:nrows
            for c=1:ncols
                n = n+1;
                subplot(nrows,ncols,(r-1)*ncols+c)

                if numBlocks==1
                    ylabel(['Sequence ',temp{r},' cell'])
                    if strcmp(plt.split,'none')
                        title([temp{c},' trials (',num2str(numTrials(c)),'t)'])
                    else
                        title([temp{ceil(c/numSplits)},' trials, ',temp2{n},' (',num2str(numTrials(c)),'t)'])
                    end
                else
                    if r==1 && c==1
                        ylabel(['Sequence ',temp1{1},' cell'])
                    end
                    title(['Block ',num2str(c)])
                end
            end
        end
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=AonlyXonly, 100t, block-wise'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AonlyXonly_100t_ktile.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AonlyXonly_100t_ktile.png']);
    

       
%     % AX, f100t, 8tile, oddeven, 13
%     this_block = 5;
%     plt = struct(); plt.split = 'odd_even';
%     F = snakePlots(traces_100t{this_block},{find(passed_100t{this_block}.AW.A==1),find(passed_100t{this_block}.AW.X==1)},[],[1;3],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=AX, 100t, block=',num2str(this_block)])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f100t_8tile_oddeven_13.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f100t_8tile_oddeven_13.png']);
%     
%     % AXonly, f100t, 8tile, oddeven, 13
%     this_block = 1;
%     plt = struct(); plt.split = 'odd_even';
%     F = snakePlots(traces_100t{this_block},{find(passed_100t{this_block}.AW.Aonly==1),find(passed_100t{this_block}.AW.Xonly==1)},[],[1;3],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=AXonly, 100t, block=',num2str(this_block)])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AXonly_f100t_8tile_oddeven_13.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AXonly_f100t_8tile_oddeven_13.png']);
%     
%     % A, 100t, wholeAinA
%     plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
%     temp = find(passed_all.AW.A==1); [~,temp1] = nanmax(avgTraces_all.A(temp,p.general.bins_analysisWindow),[],2); [~,temp2] = sort(temp1);
%     F = snakePlots(traces_100t,[],{temp(temp2)},[],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=A, 100t, sorted by whole-session sequence A in A trials'])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_wholeAinA.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_wholeAinA.png']);
%     
%     % A, 100t, wholeAinA (Aonly) %MB20221011 improvised
%     plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
%     temp = find(passed_all.AW.Aonly==1); [~,temp1] = nanmax(avgTraces_all.A(temp,p.general.bins_analysisWindow),[],2); [~,temp2] = sort(temp1);
%     F = snakePlots(traces_100t,[],{temp(temp2)},[],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=Aonly, 100t, sorted by whole-session sequence A in A trials'])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_wholeAinA_Aonly.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_wholeAinA_Aonly.png']);
% 
%     % A, 100t, 4n
%     this_block = 4;
%     plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
%     F = snakePlots(traces_100t,{find(passed_all.AW.A==1)},[],[this_block,NaN],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=A, 100t'])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_AinA_',num2str(this_block),'n.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_AinA_',num2str(this_block),'n.png']);
% 
%     % A, 100t, 4n (Aonly) %MB20221011 improvised
%     this_block = 4;
%     plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
%     F = snakePlots(traces_100t,{find(passed_all.AW.Aonly==1)},[],[this_block,NaN],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=Aonly, 100t'])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_AinA_',num2str(this_block),'n_Aonly.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_AinA_',num2str(this_block),'n_Aonly.png']);
%      
end

if ops.tngn.do_60t && ops.tngn.do_allTrials
    
    % A, 60t, wholeAinA
    plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
    temp = find(passed_all.AW.A==1); [~,temp1] = nanmax(avgTraces_all.A(temp,p.general.bins_analysisWindow),[],2); [~,temp2] = sort(temp1);
    F = snakePlots(traces_60t,[],{temp(temp2)},[],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=A, 60t, sorted by whole-session sequence A in A trials'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_wholeAinA.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_wholeAinA.png']);
    
    % X, 60t, wholeXinX
    plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "X"; 
    temp = find(passed_all.AW.X==1); [~,temp1] = nanmax(avgTraces_all.X(temp,p.general.bins_analysisWindow),[],2); [~,temp2] = sort(temp1);
    F = snakePlots(traces_60t,[],{temp(temp2)},[],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=X, 60t, sorted by whole-session sequence X in X trials'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_wholeXinX.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_wholeXinX.png']);

    % A, 60t, 4n
    this_block = 4;
    plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
    F = snakePlots(traces_60t,{find(passed_all.AW.A==1)},[],[this_block,NaN],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=A, 60t'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_AinA_',num2str(this_block),'n.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_AinA_',num2str(this_block),'n.png']);
    
    % X, 60t, n4
    this_block = 4;
    plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "X"; 
    F = snakePlots(traces_60t,{find(passed_all.AW.X==1)},[],[NaN,this_block],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=X, 60t'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_XinX_n',num2str(this_block),'.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_XinX_n',num2str(this_block),'.png']);
        
    if ~ops.tngn.skip_boringSnakePlots

        % A, 60t, wholeAinX
        plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
        temp = find(passed_all.AW.A==1); [~,temp1] = nanmax(avgTraces_all.X(temp,p.general.bins_analysisWindow),[],2); [~,temp2] = sort(temp1);
        F = snakePlots(traces_60t,[],{temp(temp2)},[],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=A, 60t, sorted by whole-session sequence A in X trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_wholeAinX.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_wholeAinX.png']);

        % X, 60t, wholeXinA
        plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "X"; 
        temp = find(passed_all.AW.X==1); [~,temp1] = nanmax(avgTraces_all.A(temp,p.general.bins_analysisWindow),[],2); [~,temp2] = sort(temp1);
        F = snakePlots(traces_60t,[],{temp(temp2)},[],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=X, 60t, sorted by whole-session sequence X in A trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_wholeXinA.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_wholeXinA.png']);

        % A, 60t, n4
        this_block = 4;
        plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
        F = snakePlots(traces_60t,{find(passed_all.AW.A==1)},[],[NaN,this_block],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=A, 60t'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_AinX_n',num2str(this_block),'.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_AinX_n',num2str(this_block),'.png']);

        % X, 60t, 4n
        this_block = 4;
        plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "X"; 
        F = snakePlots(traces_60t,{find(passed_all.AW.X==1)},[],[this_block,NaN],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tngn.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=X, 60t'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_XinA_',num2str(this_block),'n.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_XinA_',num2str(this_block),'n.png']);
    end
end

disp(['--- Saved snake plots to ',path.filepart_outX,'plots/SnakePlots.'])


%% Save and return

if ops.tngn.do_allTrials
    tngn_all.shuffling = shuffling_all;
    tngn_all.firingField = firingField_all;
    tngn_all.passed = passed_all;
    tngn_all.passed_stats = passed_stats_all;
    tngn_all.prop = prop;
    tngn_all.p = p;
    tngn_all = orderfields(tngn_all);
    save([path.filepart_out,'tngn_all.mat'],'tngn_all','-v7.3');
    tngn_all_cmpr.passed = passed_all.AW;
    tngn_all_cmpr.passed_stats = passed_stats_all.AW;
    save([path.filepart_out,'tngn_all_cmpr.mat'],'tngn_all_cmpr','-v7.3');
    disp(['--- Saved tngn_all file as ',[path.filepart_out,'tngn_all.mat'],'.'])
    
%     tngn_all.shuffling = shuffling_all;
%     tngn_all.firingField = firingField_all;
%     tngn_all.passed = passed_all;
%     tngn_all.passed_stats = passed_stats_all;
%     tngn_all.prop = prop;
%     tngn_all.p = p;
%     tngn_all = orderfields(tngn_all);
%     save([path.filepart_out,'tngn_all.mat'],'tngn_all','-v7.3');
%     tngn_all_cmpr.passed = passed_all.AW;
%     tngn_all_cmpr.passed_stats = passed_stats_all.AW;
%     save([path.filepart_out,'tngn_all_cmpr.mat'],'tngn_all_cmpr','-v7.3');
%     disp(['--- Saved tngn_all file as ',[path.filepart_out,'tngn_all.mat'],'.'])
end
if ops.tngn.do_60t || ops.tngn.do_60t_onlyFirst
    tngn_60t = {};
    if ~ops.tngn.do_60t
        temp = 1;
    else
        temp = floor(prop.numTrials/60);
    end
    for i=1:temp
%         tngn_60t{i}.trials = trials_60t{i};
%         tngn_60t{i}.traces = traces_60t{i};
%         tngn_60t{i}.avgTraces = avgTraces_60t{i};
%         tngn_60t{i}.normAvgTraces = normAvgTraces_60t{i};
        tngn_60t{i}.shuffling = shuffling_60t{i};
        tngn_60t{i}.firingField = firingField_60t{i};
        tngn_60t{i}.passed = passed_60t{i};
        tngn_60t{i}.passed_stats = passed_stats_60t{i};
%         tngn_60t{i}.passed_nneg = passed_nneg_60t{i};
%         tngn_60t{i}.passed_stats_nneg = passed_stats_nneg_60t{i};
        tngn_60t{i}.prop = prop;
        tngn_60t{i}.p = p;
        tngn_60t{i} = orderfields(tngn_60t{i});
    end
    save([path.filepart_outX,'tngn_60t.mat'],'tngn_60t','-v7.3');
    disp(['--- Saved tngn_60t file as ',[path.filepart_outX,'tngn_60t.mat'],'.'])
end
if ops.tngn.do_100t || ops.tngn.do_100t_onlyFirst
    tngn_100t = {};
    if ~ops.tngn.do_100t
        temp = 1;
    else
        temp = floor(prop.numTrials/100);
    end
    for i=1:temp
%         tngn_100t{i}.trials = trials_100t{i};
%         tngn_100t{i}.traces = traces_100t{i};
%         tngn_100t{i}.avgTraces = avgTraces_100t{i};
%         tngn_100t{i}.normAvgTraces = normAvgTraces_100t{i};
        tngn_100t{i}.shuffling = shuffling_100t{i};
        tngn_100t{i}.firingField = firingField_100t{i};
        tngn_100t{i}.passed = passed_100t{i};
        tngn_100t{i}.passed_stats = passed_stats_100t{i};
%         tngn_100t{i}.passed_nneg = passed_nneg_100t{i};
%         tngn_100t{i}.passed_stats_nneg = passed_stats_nneg_100t{i};
        tngn_100t{i}.prop = prop;
        tngn_100t{i}.p = p;
        tngn_100t{i} = orderfields(tngn_100t{i});
        tngn_100t_cmpr{i}.passed = passed_100t{i}.AW;
        tngn_100t_cmpr{i}.passed_stats = passed_stats_100t{i}.AW;
%         tngn_100t_cmpr{i}.passed_nneg = passed_nneg_100t{i}.AW;
%         tngn_100t_cmpr{i}.passed_stats_nneg = passed_stats_nneg_100t{i}.AW;
    end
    save([path.filepart_outX,'tngn_100t.mat'],'tngn_100t','-v7.3');
    save([path.filepart_out,'tngn_100t_cmpr.mat'],'tngn_100t_cmpr','-v7.3');
    disp(['--- Saved tngn_100t file as ',[path.filepart_outX,'tngn_100t.mat'],'.'])
end
if ops.tngn.do_allTrials
    tngn = tngn_all;
elseif ops.tngn.do_60t || ops.tngn.do_60t_onlyFirst
    tngn = tngn_60t{1};
elseif ops.tngn.do_100t || ops.tngn.do_100t_onlyFirst
    tngn = tngn_100t{1};
end

% stimVersion
if ops.tngn.do_allTrials_stimVersion
    tngn_all_stimVersion.shuffling = shuffling_all_stimVersion;
    tngn_all_stimVersion.firingField = firingField_all_stimVersion;
    tngn_all_stimVersion.passed_catch = passed_all_stimVersion_catch;
    tngn_all_stimVersion.passed_stats_catch = passed_stats_all_stimVersion_catch;
    tngn_all_stimVersion.passed_stim = passed_all_stimVersion_stim;
    tngn_all_stimVersion.passed_stats_stim = passed_stats_all_stimVersion_stim;
    tngn_all_stimVersion.prop = prop;
    tngn_all_stimVersion.p = p;
    tngn_all_stimVersion = orderfields(tngn_all_stimVersion);
    save([path.filepart_outX,'tngn_all_stimVersion.mat'],'tngn_all_stimVersion','-v7.3');
    disp(['--- Saved tngn_all_stimVersion file as ',[path.filepart_outX,'tngn_all_stimVersion.mat'],'.'])    
end
if ops.tngn.do_100t_stimVersion
    tngn_100t_stimVersion = {};
    for i=1:floor(prop.numTrials/100)
        tngn_100t_stimVersion{i}.shuffling = shuffling_100t_stimVersion{i};
        tngn_100t_stimVersion{i}.firingField = firingField_100t_stimVersion{i};
        tngn_100t_stimVersion{i}.passed_catch = passed_100t_stimVersion_catch{i};
        tngn_100t_stimVersion{i}.passed_stats_catch = passed_stats_100t_stimVersion_catch{i};
        tngn_100t_stimVersion{i}.passed_stim = passed_100t_stimVersion_stim{i};
        tngn_100t_stimVersion{i}.passed_stats_stim = passed_stats_100t_stimVersion_stim{i};
        tngn_100t_stimVersion{i}.prop = prop;
        tngn_100t_stimVersion{i}.p = p;
        tngn_100t_stimVersion{i} = orderfields(tngn_100t_stimVersion{i});
    end
    save([path.filepart_outX,'tngn_100t_stimVersion.mat'],'tngn_100t_stimVersion','-v7.3');
    disp(['--- Saved tngn_100t_stimVersion file as ',[path.filepart_outX,'tngn_100t_stimVersion.mat'],'.'])
end
if ops.tngn.do_allTrials_stimVersion
    tngn_stimVersion = tngn_all_stimVersion;
elseif ops.tngn.do_100t_stimVersion
    tngn_stimVersion = tngn_100t_stimVersion{1};
end

% try
%     idx = 1;
%     population = find(tngn.passed.AW.A==1);
%     F = singleCellFigure(info,p,tngn,population(idx),sync_beh);
%     if ~exist([path.filepart_outX,'plots/SingleCellFigure/'],'dir')
%         mkdir([path.filepart_outX,'plots/SingleCellFigure/']);
%     end
%     saveas(F,[path.filepart_outX,'plots/SingleCellFigure/',info.animal,'_',info.date,'_1.png']);
% catch
% end
if ops.close_figures
    close all;
end
%end


