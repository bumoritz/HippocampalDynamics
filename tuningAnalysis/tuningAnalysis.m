function [tng,nft_binned] = tuningAnalysis(info,iscell,ops,p,path,task,s2p_meta,sync_beh,act)
%sync_beh = paq_beh;

%% Preparations

disp('--- Preparations')

% pre-process activity measure
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.tng,p,sync_beh,iscell);
bl_trial_frames_binned = [prop.trial_frames_binned(1,:)-[1:5]';prop.trial_frames_binned(1:10,:)];
nft_bl_binned = reshape(nf_binned(:,bl_trial_frames_binned),[],size(bl_trial_frames_binned,1),size(bl_trial_frames_binned,2));
nft_bl_avg_binned = nanmean(nft_bl_binned,3);

% create trials, traces, avgTraces and normAvgTraces structs
if ops.tng.do_allTrials
    trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
    [traces_all,avgTraces_all,normAvgTraces_all] = createTracesStructs(nft_binned,trials_all);
end
if ops.tng.do_60t || ops.tng.do_60t_onlyFirst
    trials_60t = {}; traces_60t = {}; avgTraces_60t = {}; normAvgTraces_60t = {};
    if ~ops.tng.do_60t
        temp = 1;
    else
        temp = floor(prop.numTrials/60);
    end
    for i=1:temp
        trials_60t{i} = createTrialsStruct(task,(i-1)*60+1:i*60);
        [traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i}] = createTracesStructs(nft_binned,trials_60t{i});
    end
end
if ops.tng.do_100t || ops.tng.do_100t_onlyFirst
    trials_100t = {}; traces_100t = {}; avgTraces_100t = {}; normAvgTraces_100t = {};
    if ~ops.tng.do_100t
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
if ops.tng.do_allTrials_stimVersion
    trials_all_stimVersion = createTrialsStruct(task,1:prop.numTrials_incl,true);
    [traces_all_stimVersion,avgTraces_all_stimVersion,normAvgTraces_all_stimVersion] = createTracesStructs(nft_binned,trials_all_stimVersion);
end
if ops.tng.do_100t_stimVersion
    trials_100t_stimVersion = {}; traces_100t_stimVersion = {}; avgTraces_100t_stimVersion = {}; normAvgTraces_100t_stimVersion = {};
    for i=1:floor(prop.numTrials/100)
        trials_100t_stimVersion{i} = createTrialsStruct(task,(i-1)*100+1:i*100,true);
        [traces_100t_stimVersion{i},avgTraces_100t_stimVersion{i},normAvgTraces_100t_stimVersion{i}] = createTracesStructs(nft_binned,trials_100t_stimVersion{i});
    end
end


%% Shuffling

disp('--- Shuffling')

if ops.tng.do_allTrials
    if p.tng.analysisBlockRestrictedNf
        this_binRange = [prop.sync_all_binned(1)-length(p.general.bins_pre),prop.sync_all_binned(end)+length(p.general.t_binned)-length(p.general.bins_pre)];
    else
        this_binRange = [1,size(nf_binned,2)];
    end
    shuffling_all.A_AW = shufflingCriterion(prop,nf_binned,trials_all.stimuli.A,traces_all.A,p.general.bins_analysisWindow,p,this_binRange);
    shuffling_all.X_AW = shufflingCriterion(prop,nf_binned,trials_all.stimuli.X,traces_all.X,p.general.bins_analysisWindow,p,this_binRange);
    if ops.tng.do_eventWiseAnalysis
        shuffling_all.A_odour1 = shufflingCriterion(prop,nf_binned,trials_all.stimuli.A,traces_all.A,p.general.bins_odour1Window,p,this_binRange);
        shuffling_all.X_odour1 = shufflingCriterion(prop,nf_binned,trials_all.stimuli.X,traces_all.X,p.general.bins_odour1Window,p,this_binRange);
        shuffling_all.B_odour2 = shufflingCriterion(prop,nf_binned,trials_all.stimuli.B,traces_all.B,p.general.bins_odour2Window,p,this_binRange);
        shuffling_all.Y_odour2 = shufflingCriterion(prop,nf_binned,trials_all.stimuli.Y,traces_all.Y,p.general.bins_odour2Window,p,this_binRange);
        shuffling_all.rew_reward = shufflingCriterion(prop,nf_binned,trials_all.outcome.rew,traces_all.rew,p.general.bins_rewardWindow,p,this_binRange);
        shuffling_all.norew_reward = shufflingCriterion(prop,nf_binned,trials_all.outcome.norew,traces_all.norew,p.general.bins_rewardWindow,p,this_binRange);
    end
end
if ops.tng.do_60t || ops.tng.do_60t_onlyFirst
    shuffling_60t = {};
    if ~ops.tng.do_60t
        temp = 1;
    else
        temp = floor(prop.numTrials/60);
    end
    for i=1:temp
        if p.tng.analysisBlockRestrictedNf
            this_binRange = [prop.sync_all_binned((i-1)*60+1)-length(p.general.bins_pre),prop.sync_all_binned(i*60)+length(p.general.t_binned)-length(p.general.bins_pre)];
        else
            this_binRange = [1,size(nf_binned,2)];
        end
        shuffling_60t{i}.A_AW = shufflingCriterion(prop,nf_binned,trials_60t{i}.stimuli.A,traces_60t{i}.A,p.general.bins_analysisWindow,p,this_binRange);
        shuffling_60t{i}.X_AW = shufflingCriterion(prop,nf_binned,trials_60t{i}.stimuli.X,traces_60t{i}.X,p.general.bins_analysisWindow,p,this_binRange);
        if ops.tng.do_eventWiseAnalysis
            shuffling_60t{i}.A_odour1 = shufflingCriterion(prop,nf_binned,trials_60t{i}.stimuli.A,traces_60t{i}.A,p.general.bins_odour1Window,p,this_binRange);
            shuffling_60t{i}.X_odour1 = shufflingCriterion(prop,nf_binned,trials_60t{i}.stimuli.X,traces_60t{i}.X,p.general.bins_odour1Window,p,this_binRange);
            shuffling_60t{i}.B_odour2 = shufflingCriterion(prop,nf_binned,trials_60t{i}.stimuli.B,traces_60t{i}.B,p.general.bins_odour2Window,p,this_binRange);
            shuffling_60t{i}.Y_odour2 = shufflingCriterion(prop,nf_binned,trials_60t{i}.stimuli.Y,traces_60t{i}.Y,p.general.bins_odour2Window,p,this_binRange);
            shuffling_60t{i}.rew_reward = shufflingCriterion(prop,nf_binned,trials_60t{i}.outcome.rew,traces_60t{i}.rew,p.general.bins_rewardWindow,p,this_binRange);
            shuffling_60t{i}.norew_reward = shufflingCriterion(prop,nf_binned,trials_60t{i}.outcome.norew,traces_60t{i}.norew,p.general.bins_rewardWindow,p,this_binRange);
        end
    end
end
if ops.tng.do_100t || ops.tng.do_100t_onlyFirst
    shuffling_100t = {};
    if ~ops.tng.do_100t
        temp = 1;
    else
        temp = floor(prop.numTrials/100);
    end
    for i=1:temp
        if p.tng.analysisBlockRestrictedNf
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
        if ops.tng.do_eventWiseAnalysis
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
if (ops.tng.do_allTrials_stimVersion || ops.tng.do_100t_stimVersion)
    nf_binned_stimVersion = nf_binned;
    for i=1:size(prop.trial_frames_binned,2)
        if task.var(i)==1
            nf_binned_stimVersion(:,prop.trial_frames_binned(11:56,i)) = NaN;
        end
    end  
end
if ops.tng.do_allTrials_stimVersion
    if p.tng.analysisBlockRestrictedNf
        this_binRange = [prop.sync_all_binned(1)-length(p.general.bins_pre),prop.sync_all_binned(end)+length(p.general.t_binned)-length(p.general.bins_pre)];
    else
        this_binRange = [1,size(nf_binned_stimVersion,2)];
    end
    shuffling_all_stimVersion.Acatch_AW = shufflingCriterion(prop,nf_binned_stimVersion,trials_all_stimVersion.stimuli_var0.A,traces_all_stimVersion.A_var0,p.general.bins_analysisWindow,p,this_binRange);
    shuffling_all_stimVersion.Astim_AW = shufflingCriterion(prop,nf_binned_stimVersion,trials_all_stimVersion.stimuli_var1.A,traces_all_stimVersion.A_var1,p.general.bins_analysisWindow,p,this_binRange);
    shuffling_all_stimVersion.Xcatch_AW = shufflingCriterion(prop,nf_binned_stimVersion,trials_all_stimVersion.stimuli_var0.X,traces_all_stimVersion.X_var0,p.general.bins_analysisWindow,p,this_binRange);
    shuffling_all_stimVersion.Xstim_AW = shufflingCriterion(prop,nf_binned_stimVersion,trials_all_stimVersion.stimuli_var1.X,traces_all_stimVersion.X_var1,p.general.bins_analysisWindow,p,this_binRange);
end
if ops.tng.do_100t_stimVersion
    shuffling_100t_stimVersion = {};
    for i=1:floor(prop.numTrials/100)
        if p.tng.analysisBlockRestrictedNf
            try
                this_binRange = [prop.sync_all_binned((i-1)*100+1)-length(p.general.bins_pre),prop.sync_all_binned(i*100)+length(p.general.t_binned)-length(p.general.bins_pre)];
            catch
                this_binRange = [prop.sync_all_binned((i-1)*100+1)-length(p.general.bins_pre),prop.sync_all_binned(end-2)+length(p.general.t_binned)-length(p.general.bins_pre)];
            end
        else
            this_binRange = [1,size(nf_binned_stimVersion,2)];
        end
        shuffling_100t_stimVersion{i}.Acatch_AW = shufflingCriterion(prop,nf_binned_stimVersion,trials_100t_stimVersion{i}.stimuli_var0.A,traces_100t_stimVersion{i}.A_var0,p.general.bins_analysisWindow,p,this_binRange);
        shuffling_100t_stimVersion{i}.Astim_AW = shufflingCriterion(prop,nf_binned_stimVersion,trials_100t_stimVersion{i}.stimuli_var1.A,traces_100t_stimVersion{i}.A_var1,p.general.bins_analysisWindow,p,this_binRange);
        shuffling_100t_stimVersion{i}.Xcatch_AW = shufflingCriterion(prop,nf_binned_stimVersion,trials_100t_stimVersion{i}.stimuli_var0.X,traces_100t_stimVersion{i}.X_var0,p.general.bins_analysisWindow,p,this_binRange);        
        shuffling_100t_stimVersion{i}.Xstim_AW = shufflingCriterion(prop,nf_binned_stimVersion,trials_100t_stimVersion{i}.stimuli_var1.X,traces_100t_stimVersion{i}.X_var1,p.general.bins_analysisWindow,p,this_binRange);        
    end
end


%% Firing field properties

disp('--- Calculating firing field properties')

if ops.tng.do_allTrials
    firingField_all.A_AW = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'AW','A','X',prop,p);
    firingField_all.X_AW = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'AW','X','A',prop,p);
    if ops.tng.do_eventWiseAnalysis
        firingField_all.A_odour1 = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'odour1','A','X',prop,p);
        firingField_all.X_odour1 = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'odour1','X','A',prop,p);
        firingField_all.B_odour2 = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'odour2','B','Y',prop,p);
        firingField_all.Y_odour2 = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'odour2','Y','B',prop,p);
        firingField_all.rew_reward = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'reward','rew','norew',prop,p);
        firingField_all.norew_reward = firingFieldProperties(traces_all,avgTraces_all,normAvgTraces_all,'reward','norew','rew',prop,p);
    end
end
if ops.tng.do_60t || ops.tng.do_60t_onlyFirst
    firingField_60t = {};
    if ~ops.tng.do_60t
        temp = 1;
    else
        temp = floor(prop.numTrials/60);
    end
    for i=1:temp
        firingField_60t{i}.A_AW = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'AW','A','X',prop,p);
        firingField_60t{i}.X_AW = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'AW','X','A',prop,p);
        if ops.tng.do_eventWiseAnalysis
            firingField_60t{i}.A_odour1 = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'odour1','A','X',prop,p);
            firingField_60t{i}.X_odour1 = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'odour1','X','A',prop,p);
            firingField_60t{i}.B_odour2 = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'odour2','B','Y',prop,p);
            firingField_60t{i}.Y_odour2 = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'odour2','Y','B',prop,p);
            firingField_60t{i}.rew_reward = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'reward','rew','norew',prop,p);
            firingField_60t{i}.norew_reward = firingFieldProperties(traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i},'reward','norew','rew',prop,p);
        end
    end
end
if ops.tng.do_100t || ops.tng.do_100t_onlyFirst
    firingField_100t = {};
    if ~ops.tng.do_100t
        temp = 1;
    else
        temp = floor(prop.numTrials/100);
    end
    for i=1:temp
        firingField_100t{i}.A_AW = firingFieldProperties(traces_100t{i},avgTraces_100t{i},normAvgTraces_100t{i},'AW','A','X',prop,p);
        firingField_100t{i}.X_AW = firingFieldProperties(traces_100t{i},avgTraces_100t{i},normAvgTraces_100t{i},'AW','X','A',prop,p);
        if ops.tng.do_eventWiseAnalysis
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
if ops.tng.do_allTrials_stimVersion
    firingField_all_stimVersion.Acatch_AW = firingFieldProperties(traces_all_stimVersion,avgTraces_all_stimVersion,normAvgTraces_all_stimVersion,'AW','A_var0','X_var0',prop,p);
    firingField_all_stimVersion.Astim_AW = firingFieldProperties(traces_all_stimVersion,avgTraces_all_stimVersion,normAvgTraces_all_stimVersion,'AW','A_var1','X_var1',prop,p);
    firingField_all_stimVersion.Xcatch_AW = firingFieldProperties(traces_all_stimVersion,avgTraces_all_stimVersion,normAvgTraces_all_stimVersion,'AW','X_var0','A_var0',prop,p);
    firingField_all_stimVersion.Xstim_AW = firingFieldProperties(traces_all_stimVersion,avgTraces_all_stimVersion,normAvgTraces_all_stimVersion,'AW','X_var1','A_var1',prop,p);
end
if ops.tng.do_100t_stimVersion
    firingField_100t_stimVersion = {};
    for i=1:floor(prop.numTrials/100)
        firingField_100t_stimVersion{i}.Acatch_AW = firingFieldProperties(traces_100t_stimVersion{i},avgTraces_100t_stimVersion{i},normAvgTraces_100t_stimVersion{i},'AW','A_var0','X_var0',prop,p);
        firingField_100t_stimVersion{i}.Astim_AW = firingFieldProperties(traces_100t_stimVersion{i},avgTraces_100t_stimVersion{i},normAvgTraces_100t_stimVersion{i},'AW','A_var1','X_var1',prop,p);
        firingField_100t_stimVersion{i}.Xcatch_AW = firingFieldProperties(traces_100t_stimVersion{i},avgTraces_100t_stimVersion{i},normAvgTraces_100t_stimVersion{i},'AW','X_var0','A_var0',prop,p);
        firingField_100t_stimVersion{i}.Xstim_AW = firingFieldProperties(traces_100t_stimVersion{i},avgTraces_100t_stimVersion{i},normAvgTraces_100t_stimVersion{i},'AW','X_var1','A_var1',prop,p);
    end
end


%% Inhibition/artefact evasion

% if ops.tng.do_allTrials
%     firingField_all.A_1s.negativeGoing = double(nanmean(avgTraces_all.A(:,p.general.bins_1st1s),2) - nanmean(avgTraces_all.A(:,p.general.bins_base1s),2) <0);
%     firingField_all.A_1s.negativeGoing(~prop.iscell==1) = NaN;
%     firingField_all.X_1s.negativeGoing = double(nanmean(avgTraces_all.X(:,p.general.bins_1st1s),2) - nanmean(avgTraces_all.X(:,p.general.bins_base1s),2) <0);
%     firingField_all.X_1s.negativeGoing(~prop.iscell==1) = NaN;
% end
% if ops.tng.do_100t || ops.tng.do_100t_onlyFirst
%     if ~ops.tng.do_100t
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

if ops.tng.do_allTrials
    [passed_all.AW,passed_stats_all.AW] = tuningCriteria(shuffling_all,firingField_all,'AW','A','X',p,prop);
%     [passed_nneg_all.AW,passed_stats_nneg_all.AW] = tuningCriteria_nneg(shuffling_all,firingField_all,'AW','A','X',p,prop);
    if ops.tng.do_eventWiseAnalysis
        [passed_all.odour1,passed_stats_all.odour1] = tuningCriteria(shuffling_all,firingField_all,'odour1','A','X',p,prop);
        [passed_all.odour2,passed_stats_all.odour2] = tuningCriteria(shuffling_all,firingField_all,'odour2','B','Y',p,prop);
        [passed_all.reward,passed_stats_all.reward] = tuningCriteria(shuffling_all,firingField_all,'reward','rew','norew',p,prop);
    end
end
if ops.tng.do_60t || ops.tng.do_60t_onlyFirst
    passed_60t = {}; passed_stats_60t = {};
    if ~ops.tng.do_60t
        temp = 1;
    else
        temp = floor(prop.numTrials/60);
    end
    for i=1:temp
        [passed_60t{i}.AW,passed_stats_60t{i}.AW] = tuningCriteria(shuffling_60t{i},firingField_60t{i},'AW','A','X',p,prop);
        [passed_nneg_60t{i}.AW,passed_stats_nneg_60t{i}.AW] = tuningCriteria_nneg(shuffling_60t{i},firingField_60t{i},'AW','A','X',p,prop);
        if ops.tng.do_eventWiseAnalysis
            [passed_60t{i}.odour1,passed_stats_60t{i}.odour1] = tuningCriteria(shuffling_60t{i},firingField_60t{i},'odour1','A','X',p,prop);
            [passed_60t{i}.odour2,passed_stats_60t{i}.odour2] = tuningCriteria(shuffling_60t{i},firingField_60t{i},'odour2','B','Y',p,prop);
            [passed_60t{i}.reward,passed_stats_60t{i}.reward] = tuningCriteria(shuffling_60t{i},firingField_60t{i},'reward','rew','norew',p,prop);
        end
    end
end
if ops.tng.do_100t || ops.tng.do_100t_onlyFirst
    passed_100t = {}; passed_stats_100t = {};
    if ~ops.tng.do_100t
        temp = 1;
    else
        temp = floor(prop.numTrials/100);
    end
    for i=1:temp
        [passed_100t{i}.AW,passed_stats_100t{i}.AW] = tuningCriteria(shuffling_100t{i},firingField_100t{i},'AW','A','X',p,prop);
%         [passed_nneg_100t{i}.AW,passed_stats_nneg_100t{i}.AW] = tuningCriteria_nneg(shuffling_100t{i},firingField_100t{i},'AW','A','X',p,prop);
        if ops.tng.do_eventWiseAnalysis
            [passed_100t{i}.odour1,passed_stats_100t{i}.odour1] = tuningCriteria(shuffling_100t{i},firingField_100t{i},'odour1','A','X',p,prop);
            [passed_100t{i}.odour2,passed_stats_100t{i}.odour2] = tuningCriteria(shuffling_100t{i},firingField_100t{i},'odour2','B','Y',p,prop);
            [passed_100t{i}.reward,passed_stats_100t{i}.reward] = tuningCriteria(shuffling_100t{i},firingField_100t{i},'reward','rew','norew',p,prop);
        end
    end
end

% stimVersion
if ops.tng.do_allTrials_stimVersion
    [passed_all_stimVersion_catch.AW,passed_stats_all_stimVersion_catch.AW] = tuningCriteria(shuffling_all_stimVersion,firingField_all_stimVersion,'AW','Acatch','Xcatch',p,prop);
    [passed_all_stimVersion_stim.AW,passed_stats_all_stimVersion_stim.AW] = tuningCriteria(shuffling_all_stimVersion,firingField_all_stimVersion,'AW','Astim','Xstim',p,prop);
end
if ops.tng.do_100t_stimVersion
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
if strcmp(p.tng.activityMeasure,"casc_50c_beh")
    this_activity_measure = 'casc_50c';
elseif strcmp(p.tng.activityMeasure,"casc_50g_beh")
    this_activity_measure = 'casc_50g';
elseif strcmp(p.tng.activityMeasure,"casc_100c_beh")
    this_activity_measure = 'casc_100c';
elseif strcmp(p.tng.activityMeasure,"casc_100g_beh")
    this_activity_measure = 'casc_100g';
elseif strcmp(p.tng.activityMeasure,"casc_200g_beh")
    this_activity_measure = 'casc_200g';
elseif strcmp(p.tng.activityMeasure,"spks_beh")
    this_activity_measure = 'spks';
elseif strcmp(p.tng.activityMeasure,"dFF_beh")
    this_activity_measure = 'dFF';
elseif strcmp(p.tng.activityMeasure,"spksn_beh")
    this_activity_measure = 'spksn';
elseif strcmp(p.tng.activityMeasure,"spksnn_beh")
    this_activity_measure = 'spksnn';
elseif strcmp(p.tng.activityMeasure,"dFFn_beh")
    this_activity_measure = 'dFFn';
elseif strcmp(p.tng.activityMeasure,"dFFnn_beh")
    this_activity_measure = 'dFFnn';
elseif strcmp(p.tng.activityMeasure,"F_beh")
    this_activity_measure = 'F';
end

if ops.tng.do_allTrials

%     % AX, all trials, 8tile, oddeven, 13
%     plt = struct(); plt.split = 'odd_even';
%     F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[1;3],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
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
        this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=AXonly, all trials'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AXonly_allTrials_8tile_oddeven_13.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AXonly_allTrials_8tile_oddeven_13.png']);
     
%    % AX, all trials, 8tile, oddeven, 13, two-sided
%     plt = struct(); plt.split = 'odd_even'; plt.normalisationType = 'two-sided'; plt.normalisationWindow = 'all';
%     F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[1;3],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=AX, all trials'])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_twosided_oddeven_13.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_twosided_oddeven_13.png']);
% 
%     % AXonly, all trials, 8tile, oddeven, 13, two-sided
%     plt = struct(); plt.split = 'odd_even'; plt.normalisationType = 'two-sided'; plt.normalisationWindow = 'all';
%     F = snakePlots(traces_all,{find(passed_all.AW.Aonly==1),find(passed_all.AW.Xonly==1)},[],[1;3],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=AXonly, all trials'])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AXonly_allTrials_8tile_twosided_oddeven_13.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AXonly_allTrials_8tile_twosided_oddeven_13.png']);

    if ~ops.tng.skip_boringSnakePlots
        
        % AX, all trials, 12tile, corrOcorrMincorr, same, 14
        plt = struct(); plt.split = 'corrO_corrM_incorr'; plt.split3_numTrials = 'same';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[1;4],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_same_14.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_same_14.png']);

        % AX, all trials, 4tile, 12
        plt = struct(); plt.split = 'none';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[1;2],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_4tile_12.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_4tile_12.png']);

        % AX, all trials, 4tile, 21
        plt = struct(); plt.split = 'none';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[2;1],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_4tile_21.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_4tile_21.png']);

        % AX, all trials, 8tile, oddeven, 31
        plt = struct(); plt.split = 'odd_even';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[3;1],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_oddeven_31.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_oddeven_31.png']);

        % AX, all trials, 8tile, corrincorr, 13
        plt = struct(); plt.split = 'corr_incorr';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[1;3],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_corrincorr_13.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_corrincorr_13.png']);

        % AX, all trials, 8tile, corrincorr, 24
        plt = struct(); plt.split = 'corr_incorr';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[2;4],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_corrincorr_24.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_corrincorr_24.png']);

        % AX, all trials, 8tile, corrincorr, 31
        plt = struct(); plt.split = 'corr_incorr';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[3;1],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_corrincorr_31.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_8tile_corrincorr_31.png']);
    
        % AX, all trials, 12tile, corrOcorrMincorr, max, 14
        plt = struct(); plt.split = 'corrO_corrM_incorr'; plt.split3_numTrials = 'max';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[1;4],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_max_14.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_max_14.png']);

        % AX, all trials, 12tile, corrOcorrMincorr, max, 41
        plt = struct(); plt.split = 'corrO_corrM_incorr'; plt.split3_numTrials = 'max';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[4;1],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_max_41.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_max_41.png']);

        % AX, all trials, 12tile, corrOcorrMincorr, same, 41
        plt = struct(); plt.split = 'corrO_corrM_incorr'; plt.split3_numTrials = 'same';
        F = snakePlots(traces_all,{find(passed_all.AW.A==1),find(passed_all.AW.X==1)},[],[4;1],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=AX, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_same_41.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_allTrials_12tile_corrOcorrMincorr_same_41.png']);
    end
end


if ops.tng.do_60t_onlyFirst
    % AX, f60t, 8tile, oddeven, 13
    plt = struct(); plt.split = 'odd_even';
    F = snakePlots(traces_60t{1},{find(passed_60t{1}.AW.A==1),find(passed_60t{1}.AW.X==1)},[],[1;3],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=AX, first 60t'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f60t_8tile_oddeven_13.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f60t_8tile_oddeven_13.png']);
    
    % AX, f60t, 8tile, oddeven, 13
    plt = struct(); plt.split = 'odd_even';
    F = snakePlots(traces_60t{1},{find(passed_60t{1}.AW.A==1),find(passed_60t{1}.AW.X==1)},[],[1;3],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=AX, first 60t'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f60t_8tile_oddeven_13.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f60t_8tile_oddeven_13.png']);

    % AX, f60t, 12tile, corrOcorrMincorr, same, 14
    plt = struct(); plt.split = 'corrO_corrM_incorr'; plt.split3_numTrials = 'same';
    F = snakePlots(traces_60t{1},{find(passed_60t{1}.AW.A==1),find(passed_60t{1}.AW.X==1)},[],[1;4],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=AX, first 60t'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f60t_12tile_corrOcorrMincorr_same_14.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f60t_12tile_corrOcorrMincorr_same_14.png']);
end

if ops.tng.do_100t
    
    % Aonly, 100t, ktile
    these_idcs = {};
    for i=1:length(passed_100t)
        these_idcs{i} = find(passed_100t{i}.AW.Aonly==1);
    end
    plt = struct(); plt.split = 'none'; plt.normalisationRef = 'plotwise'; plt.sequences = "A"; plt.trialTypes = "A";    
    [F,these_normAvgTraces_A,these_idcs_sorted_A,plt] = snakePlots(traces_100t,these_idcs,[],[1:length(passed_100t)],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
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
        this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=Xonly, 100t, block-wise'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Xonly_100t_ktile.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Xonly_100t_ktile.png']);
    
    %Aonlx and Xonly, ktile
    F = paper_figure([0,0.5,mm2inch(10*34),mm2inch(3*34)]); hold on;
    %F = default_figure([-20,0.5,20,5]); 
    nrows = 1; ncols = length(passed_100t); numBlocks = length(passed_100t); these_refs_sort = [1:length(passed_100t)];
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
                    %title(['Block ',num2str(c)])
                end
                
                % paper version
%                 xlim([6,53])
%                 set(gca,'xtick',[]);
%                 set(gca,'xticklabel',[]);
%                 set(gca,'ytick',[]);
%                 set(gca,'yticklabel',[]);
%                 set(gca,'box','off');
%                 set(gca,'xcolor','none');
%                 set(gca,'ycolor','none');
%                 box off
            end
        end
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=AonlyXonly, 100t, block-wise'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AonlyXonly_100t_ktile.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AonlyXonly_100t_ktile.png']);
    
    
    % save plot
%     save_root_fig = [path.root_summary,'figures\Fig3_fig\'];
%     save_root_png = [path.root_summary,'figures\Fig3_png\'];
%     save_root_pdf = [path.root_summary,'figures\Fig3_pdf\'];
%     save_root_txt = [path.root_summary,'figures\Fig3_txt\'];
%     savefig(F,[save_root_fig,'\Fig3_ExampleSnakePlots_Elrond_d3.fig']);
%     saveas(F,[save_root_png,'\Fig3_ExampleSnakePlots_Elrond_d3.png']);
%     set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_ExampleSnakePlots_Elrond_d3.pdf']); set(gcf,'Color',[1,1,1])

%     % AX, f100t, 8tile, oddeven, 13
%     this_block = 5;
%     plt = struct(); plt.split = 'odd_even';
%     F = snakePlots(traces_100t{this_block},{find(passed_100t{this_block}.AW.A==1),find(passed_100t{this_block}.AW.X==1)},[],[1;3],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=AX, 100t, block=',num2str(this_block)])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f100t_8tile_oddeven_13.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AX_f100t_8tile_oddeven_13.png']);
%     
%     % AXonly, f100t, 8tile, oddeven, 13
%     this_block = 1;
%     plt = struct(); plt.split = 'odd_even';
%     F = snakePlots(traces_100t{this_block},{find(passed_100t{this_block}.AW.Aonly==1),find(passed_100t{this_block}.AW.Xonly==1)},[],[1;3],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=AXonly, 100t, block=',num2str(this_block)])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AXonly_f100t_8tile_oddeven_13.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_AXonly_f100t_8tile_oddeven_13.png']);
%     
%     % A, 100t, wholeAinA
%     plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
%     temp = find(passed_all.AW.A==1); [~,temp1] = nanmax(avgTraces_all.A(temp,p.general.bins_analysisWindow),[],2); [~,temp2] = sort(temp1);
%     F = snakePlots(traces_100t,[],{temp(temp2)},[],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=A, 100t, sorted by whole-session sequence A in A trials'])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_wholeAinA.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_wholeAinA.png']);
%     
%     % A, 100t, wholeAinA (Aonly) %MB20221011 improvised
%     plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
%     temp = find(passed_all.AW.Aonly==1); [~,temp1] = nanmax(avgTraces_all.A(temp,p.general.bins_analysisWindow),[],2); [~,temp2] = sort(temp1);
%     F = snakePlots(traces_100t,[],{temp(temp2)},[],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=Aonly, 100t, sorted by whole-session sequence A in A trials'])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_wholeAinA_Aonly.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_wholeAinA_Aonly.png']);
% 
%     % A, 100t, 4n
%     this_block = 4;
%     plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
%     F = snakePlots(traces_100t,{find(passed_all.AW.A==1)},[],[this_block,NaN],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=A, 100t'])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_AinA_',num2str(this_block),'n.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_AinA_',num2str(this_block),'n.png']);
% 
%     % A, 100t, 4n (Aonly) %MB20221011 improvised
%     this_block = 4;
%     plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
%     F = snakePlots(traces_100t,{find(passed_all.AW.Aonly==1)},[],[this_block,NaN],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%         'passed=Aonly, 100t'])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_AinA_',num2str(this_block),'n_Aonly.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_100t_AinA_',num2str(this_block),'n_Aonly.png']);
%      
end

if ops.tng.do_60t && ops.tng.do_allTrials
    
    % A, 60t, wholeAinA
    plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
    temp = find(passed_all.AW.A==1); [~,temp1] = nanmax(avgTraces_all.A(temp,p.general.bins_analysisWindow),[],2); [~,temp2] = sort(temp1);
    F = snakePlots(traces_60t,[],{temp(temp2)},[],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=A, 60t, sorted by whole-session sequence A in A trials'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_wholeAinA.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_wholeAinA.png']);
    
    % X, 60t, wholeXinX
    plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "X"; 
    temp = find(passed_all.AW.X==1); [~,temp1] = nanmax(avgTraces_all.X(temp,p.general.bins_analysisWindow),[],2); [~,temp2] = sort(temp1);
    F = snakePlots(traces_60t,[],{temp(temp2)},[],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=X, 60t, sorted by whole-session sequence X in X trials'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_wholeXinX.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_wholeXinX.png']);

    % A, 60t, 4n
    this_block = 4;
    plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
    F = snakePlots(traces_60t,{find(passed_all.AW.A==1)},[],[this_block,NaN],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=A, 60t'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_AinA_',num2str(this_block),'n.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_AinA_',num2str(this_block),'n.png']);
    
    % X, 60t, n4
    this_block = 4;
    plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "X"; 
    F = snakePlots(traces_60t,{find(passed_all.AW.X==1)},[],[NaN,this_block],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=X, 60t'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_XinX_n',num2str(this_block),'.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_XinX_n',num2str(this_block),'.png']);
        
    if ~ops.tng.skip_boringSnakePlots

        % A, 60t, wholeAinX
        plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
        temp = find(passed_all.AW.A==1); [~,temp1] = nanmax(avgTraces_all.X(temp,p.general.bins_analysisWindow),[],2); [~,temp2] = sort(temp1);
        F = snakePlots(traces_60t,[],{temp(temp2)},[],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=A, 60t, sorted by whole-session sequence A in X trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_wholeAinX.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_wholeAinX.png']);

        % X, 60t, wholeXinA
        plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "X"; 
        temp = find(passed_all.AW.X==1); [~,temp1] = nanmax(avgTraces_all.A(temp,p.general.bins_analysisWindow),[],2); [~,temp2] = sort(temp1);
        F = snakePlots(traces_60t,[],{temp(temp2)},[],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=X, 60t, sorted by whole-session sequence X in A trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_wholeXinA.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_wholeXinA.png']);

        % A, 60t, n4
        this_block = 4;
        plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
        F = snakePlots(traces_60t,{find(passed_all.AW.A==1)},[],[NaN,this_block],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=A, 60t'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_AinX_n',num2str(this_block),'.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_A_60t_AinX_n',num2str(this_block),'.png']);

        % X, 60t, 4n
        this_block = 4;
        plt = struct(); plt.split = 'none'; plt.normalisationRef = 'uniform'; plt.sequences = "X"; 
        F = snakePlots(traces_60t,{find(passed_all.AW.X==1)},[],[this_block,NaN],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'passed=X, 60t'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_XinA_',num2str(this_block),'n.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_X_60t_XinA_',num2str(this_block),'n.png']);
    end
end

disp(['--- Saved snake plots to ',path.filepart_outX,'plots/SnakePlots.'])


%% Save and return

if ops.tng.do_allTrials
    tng_all.shuffling = shuffling_all;
    tng_all.firingField = firingField_all;
    tng_all.passed = passed_all;
    tng_all.passed_stats = passed_stats_all;
    tng_all.prop = prop;
    tng_all.p = p;
    tng_all = orderfields(tng_all);
    save([path.filepart_out,'tng_all.mat'],'tng_all','-v7.3');
    tng_all_cmpr.passed = passed_all.AW;
    tng_all_cmpr.passed_stats = passed_stats_all.AW;
    save([path.filepart_out,'tng_all_cmpr.mat'],'tng_all_cmpr','-v7.3');
    disp(['--- Saved tng_all file as ',[path.filepart_out,'tng_all.mat'],'.'])
    
%     tng_all.shuffling = shuffling_all;
%     tng_all.firingField = firingField_all;
%     tng_all.passed = passed_all;
%     tng_all.passed_stats = passed_stats_all;
%     tng_all.prop = prop;
%     tng_all.p = p;
%     tng_all = orderfields(tng_all);
%     save([path.filepart_out,'tng_all.mat'],'tng_all','-v7.3');
%     tng_all_cmpr.passed = passed_all.AW;
%     tng_all_cmpr.passed_stats = passed_stats_all.AW;
%     save([path.filepart_out,'tng_all_cmpr.mat'],'tng_all_cmpr','-v7.3');
%     disp(['--- Saved tng_all file as ',[path.filepart_out,'tng_all.mat'],'.'])
end
if ops.tng.do_60t || ops.tng.do_60t_onlyFirst
    tng_60t = {};
    if ~ops.tng.do_60t
        temp = 1;
    else
        temp = floor(prop.numTrials/60);
    end
    for i=1:temp
%         tng_60t{i}.trials = trials_60t{i};
%         tng_60t{i}.traces = traces_60t{i};
%         tng_60t{i}.avgTraces = avgTraces_60t{i};
%         tng_60t{i}.normAvgTraces = normAvgTraces_60t{i};
        tng_60t{i}.shuffling = shuffling_60t{i};
        tng_60t{i}.firingField = firingField_60t{i};
        tng_60t{i}.passed = passed_60t{i};
        tng_60t{i}.passed_stats = passed_stats_60t{i};
%         tng_60t{i}.passed_nneg = passed_nneg_60t{i};
%         tng_60t{i}.passed_stats_nneg = passed_stats_nneg_60t{i};
        tng_60t{i}.prop = prop;
        tng_60t{i}.p = p;
        tng_60t{i} = orderfields(tng_60t{i});
    end
    save([path.filepart_outX,'tng_60t.mat'],'tng_60t','-v7.3');
    disp(['--- Saved tng_60t file as ',[path.filepart_outX,'tng_60t.mat'],'.'])
end
if ops.tng.do_100t || ops.tng.do_100t_onlyFirst
    tng_100t = {};
    if ~ops.tng.do_100t
        temp = 1;
    else
        temp = floor(prop.numTrials/100);
    end
    for i=1:temp
%         tng_100t{i}.trials = trials_100t{i};
%         tng_100t{i}.traces = traces_100t{i};
%         tng_100t{i}.avgTraces = avgTraces_100t{i};
%         tng_100t{i}.normAvgTraces = normAvgTraces_100t{i};
        tng_100t{i}.shuffling = shuffling_100t{i};
        tng_100t{i}.firingField = firingField_100t{i};
        tng_100t{i}.passed = passed_100t{i};
        tng_100t{i}.passed_stats = passed_stats_100t{i};
%         tng_100t{i}.passed_nneg = passed_nneg_100t{i};
%         tng_100t{i}.passed_stats_nneg = passed_stats_nneg_100t{i};
        tng_100t{i}.prop = prop;
        tng_100t{i}.p = p;
        tng_100t{i} = orderfields(tng_100t{i});
        tng_100t_cmpr{i}.passed = passed_100t{i}.AW;
        tng_100t_cmpr{i}.passed_stats = passed_stats_100t{i}.AW;
%         tng_100t_cmpr{i}.passed_nneg = passed_nneg_100t{i}.AW;
%         tng_100t_cmpr{i}.passed_stats_nneg = passed_stats_nneg_100t{i}.AW;
    end
    save([path.filepart_outX,'tng_100t.mat'],'tng_100t','-v7.3');
    save([path.filepart_out,'tng_100t_cmpr.mat'],'tng_100t_cmpr','-v7.3');
    disp(['--- Saved tng_100t file as ',[path.filepart_outX,'tng_100t.mat'],'.'])
end
if ops.tng.do_allTrials
    tng = tng_all;
elseif ops.tng.do_60t || ops.tng.do_60t_onlyFirst
    tng = tng_60t{1};
elseif ops.tng.do_100t || ops.tng.do_100t_onlyFirst
    tng = tng_100t{1};
end

% stimVersion
if ops.tng.do_allTrials_stimVersion
    tng_all_stimVersion.shuffling = shuffling_all_stimVersion;
    tng_all_stimVersion.firingField = firingField_all_stimVersion;
    tng_all_stimVersion.passed_catch = passed_all_stimVersion_catch;
    tng_all_stimVersion.passed_stats_catch = passed_stats_all_stimVersion_catch;
    tng_all_stimVersion.passed_stim = passed_all_stimVersion_stim;
    tng_all_stimVersion.passed_stats_stim = passed_stats_all_stimVersion_stim;
    tng_all_stimVersion.prop = prop;
    tng_all_stimVersion.p = p;
    tng_all_stimVersion = orderfields(tng_all_stimVersion);
    save([path.filepart_outX,'tng_all_stimVersion.mat'],'tng_all_stimVersion','-v7.3');
    disp(['--- Saved tng_all_stimVersion file as ',[path.filepart_outX,'tng_all_stimVersion.mat'],'.'])    
end
if ops.tng.do_100t_stimVersion
    tng_100t_stimVersion = {};
    for i=1:floor(prop.numTrials/100)
        tng_100t_stimVersion{i}.shuffling = shuffling_100t_stimVersion{i};
        tng_100t_stimVersion{i}.firingField = firingField_100t_stimVersion{i};
        tng_100t_stimVersion{i}.passed_catch = passed_100t_stimVersion_catch{i};
        tng_100t_stimVersion{i}.passed_stats_catch = passed_stats_100t_stimVersion_catch{i};
        tng_100t_stimVersion{i}.passed_stim = passed_100t_stimVersion_stim{i};
        tng_100t_stimVersion{i}.passed_stats_stim = passed_stats_100t_stimVersion_stim{i};
        tng_100t_stimVersion{i}.prop = prop;
        tng_100t_stimVersion{i}.p = p;
        tng_100t_stimVersion{i} = orderfields(tng_100t_stimVersion{i});
    end
    save([path.filepart_outX,'tng_100t_stimVersion.mat'],'tng_100t_stimVersion','-v7.3');
    disp(['--- Saved tng_100t_stimVersion file as ',[path.filepart_outX,'tng_100t_stimVersion.mat'],'.'])
end
if ops.tng.do_allTrials_stimVersion
    tng_stimVersion = tng_all_stimVersion;
elseif ops.tng.do_100t_stimVersion
    tng_stimVersion = tng_100t_stimVersion{1};
end

% try
%     idx = 1;
%     population = find(tng.passed.AW.A==1);
%     F = singleCellFigure(info,p,tng,population(idx),sync_beh);
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


%% Make single cell figure

% load('F:\SniffinHippo\AnalysisX\Turing\Turing_20220916\Turing_20220916_tng_all_stimVersion.mat');

%%

% idx = 91; % 425; 587;
% 
% tng_all.prop = prop;
% tng_all.trials = trials_all_stimVersion;
% tng_all.traces = traces_all_stimVersion;
% tng_all.avgTraces = avgTraces_all_stimVersion;
% 
% numBlocks = 10;
% tng_all.trials.stimuli_inA_var0.AB = tng_all.trials.stimuli_inA_var0.AB(1:1*numBlocks);
% tng_all.trials.stimuli_inA_var1.AB = tng_all.trials.stimuli_inA_var1.AB(1:4*numBlocks);
% tng_all.trials.stimuli_inA_var0.AY = tng_all.trials.stimuli_inA_var0.AY(1:1*numBlocks);
% tng_all.trials.stimuli_inA_var1.AY = tng_all.trials.stimuli_inA_var1.AY(1:4*numBlocks);
% tng_all.trials.stimuli_inX_var0.XY = tng_all.trials.stimuli_inX_var0.XY(1:1*numBlocks);
% tng_all.trials.stimuli_inX_var1.XY = tng_all.trials.stimuli_inX_var1.XY(1:4*numBlocks);
% tng_all.trials.stimuli_inX_var0.XB = tng_all.trials.stimuli_inX_var0.XB(1:1*numBlocks);
% tng_all.trials.stimuli_inX_var1.XB = tng_all.trials.stimuli_inX_var1.XB(1:4*numBlocks);
% 
% F = singleCellFigure(info,p,tng_all,idx,paq_beh,'_var1'); % '_var0' or '_var1'
% title([info.animal,'-',info.date,'-idx',num2str(idx),'-stim'])
% 
% F = singleCellFigure(info,p,tng_all,idx,paq_beh,'_var0'); % '_var0' or '_var1'
% title([info.animal,'-',info.date,'-idx',num2str(idx),'-nostim'])

%%


