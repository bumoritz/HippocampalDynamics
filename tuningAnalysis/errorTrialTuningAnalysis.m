%function [ett] = errorTrialTuningAnalysis(info,iscell,ops,p,path,task,s2p_meta,sync_beh,act)
act = spks_beh; p = get_p; sync_beh = paq_beh; %spks_beh; %Fneu_beh - values_beh'; %casc_100c_beh;

%% Preparations

disp('--- Preparations')

% start parallel pool (if not started yet)
try
    parpool; % can be closed with: delete(gcp('nocreate'));
catch
end

% pre-process activity measure
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.ett,p,sync_beh,iscell);

% create trials struct
trials = createTrialsStruct(task,1:prop.numTrials);

% prepare cross-validation
rng(p.general.rgnSeed);
if strcmp(p.ett.numIncorrectTrialsAX,'same')
    these_numTestTrials_A = nanmin([length(trials.outcome.A_incorrect),length(trials.outcome.X_incorrect)]);
    these_numTestTrials_X = nanmin([length(trials.outcome.A_incorrect),length(trials.outcome.X_incorrect)]);
elseif strcmp(p.ett.numIncorrectTrialsAX,'max')
    these_numTestTrials_A = length(trials.outcome.A_incorrect);
    these_numTestTrials_X = length(trials.outcome.X_incorrect);
end
for i=1:p.ett.cvFolds
    trials_cv{i}.outcome.A_incorrect =  datasample(trials.outcome.A_incorrect,these_numTestTrials_A,'Replace',false);
    trials_cv{i}.outcome.A_correct_test = datasample(trials.outcome.A_correct,these_numTestTrials_A,'Replace',false);
    trials_cv{i}.outcome.A_correct_train = setdiff(trials.outcome.A_correct,trials_cv{i}.outcome.A_correct_test);
    trials_cv{i}.outcome.X_incorrect = datasample(trials.outcome.X_incorrect,these_numTestTrials_X,'Replace',false);
    trials_cv{i}.outcome.X_correct_test = datasample(trials.outcome.X_correct,these_numTestTrials_X,'Replace',false);
    trials_cv{i}.outcome.X_correct_train = setdiff(trials.outcome.X_correct,trials_cv{i}.outcome.X_correct_test);
end

% create traces, avgTraces and normAvgTraces structs
for i=1:p.ett.cvFolds
    [traces_cv{i},avgTraces_cv{i},normAvgTraces_cv{i}] = createTracesStructs(nft_binned,trials_cv{i});
end


%% Shuffling

disp('--- Shuffling')

if p.ett.analysisBlockRestrictedNf
    this_binRange = [prop.sync_all_binned(1)-length(p.general.bins_pre),prop.sync_all_binned(prop.numTrials)+length(p.general.t_binned)-length(p.general.bins_pre)];
else
    this_binRange = [1,size(nf_binned,2)];
end

for i=1:p.ett.cvFolds
    shuffling_cv{i}.A_AW = shufflingCriterion(prop,nf_binned,trials_cv{i}.outcome.A_correct_train,traces_cv{i}.A_correct_train,p.general.bins_analysisWindow,p,this_binRange);
    shuffling_cv{i}.X_AW = shufflingCriterion(prop,nf_binned,trials_cv{i}.outcome.X_correct_train,traces_cv{i}.X_correct_train,p.general.bins_analysisWindow,p,this_binRange);
    if ops.ett.do_eventWiseAnalysis
    end
end


%% Firing field properties

disp('--- Calculating firing field properties')

for i=1:p.ett.cvFolds
    firingField_cv{i}.A_AW = firingFieldProperties(traces_cv{i},avgTraces_cv{i},normAvgTraces_cv{i},'AW','A_correct_train','X_correct_train',prop,p);
    firingField_cv{i}.X_AW = firingFieldProperties(traces_cv{i},avgTraces_cv{i},normAvgTraces_cv{i},'AW','X_correct_train','A_correct_train',prop,p);
    if ops.ett.do_eventWiseAnalysis
    end
end


%% Applying tuning criteria

disp('--- Applying tuning criteria')

for i=1:p.ett.cvFolds
    [passed_cv{i}.AW,passed_stats_cv{i}.AW] = tuningCriteria(shuffling_cv{i},firingField_cv{i},'AW','A','X',p,prop);
    if ops.ett.do_eventWiseAnalysis
    end
end


%% Snake plots

% preparations
if ~exist([path.filepart_out,'plots/SnakePlots/'],'dir')
    mkdir([path.filepart_out,'plots/SnakePlots/']);
end
if strcmp(p.ett.activityMeasure,"casc_50c_beh")
    this_activity_measure = 'casc_50c';
elseif strcmp(p.ett.activityMeasure,"casc_50g_beh")
    this_activity_measure = 'casc_50g';
elseif strcmp(p.ett.activityMeasure,"casc_100c_beh")
    this_activity_measure = 'casc_100c';
elseif strcmp(p.ett.activityMeasure,"casc_100g_beh")
    this_activity_measure = 'casc_100g';
elseif strcmp(p.ett.activityMeasure,"casc_200g_beh")
    this_activity_measure = 'casc_200g';
elseif strcmp(p.ett.activityMeasure,"spks_beh")
    this_activity_measure = 'spks';
elseif strcmp(p.ett.activityMeasure,"dFF_beh")
    this_activity_measure = 'dFF';
end

for i=1:p.ett.cvFolds
    
    % AX, all trials, 12tile, corrTraincorrTestincorrTest, k, 14
    plt = struct(); plt.split = 'corrTrain_corrTest_incorrTest';
    F = snakePlots(traces_cv{i},{find(passed_cv{i}.AW.A==1),find(passed_cv{i}.AW.X==1)},[],[1;4],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.ett.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'passed=AX, all trials, ',num2str(p.ett.cvFolds),'-fold CV, k=',num2str(i)])
    savefig(F,[path.filepart_out,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AX_allTrials_12tile_corrTraincorrTestincorrTest_k',num2str(i),'_14.fig']);
    saveas(F,[path.filepart_out,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AX_allTrials_12tile_corrTraincorrTestincorrTest_k',num2str(i),'_14.png']);
    
end

disp(['--- Saved snake plots to ',path.filepart_out,'plots/SnakePlots.'])


%% Sequence cell activity by outcome

for i=1:p.ett.cvFolds
    sequenceCellActivity_cv{i}.A_AW = sequenceCellActivityByOutcome(traces_cv{i},avgTraces_cv{i},normAvgTraces_cv{i},firingField_cv{i}.A_AW,find(passed_cv{i}.AW.A==1),'AW','A','X',prop,p);
    sequenceCellActivity_cv{i}.A_FF = sequenceCellActivityByOutcome(traces_cv{i},avgTraces_cv{i},normAvgTraces_cv{i},firingField_cv{i}.A_AW,find(passed_cv{i}.AW.A==1),'FF','A','X',prop,p);
    sequenceCellActivity_cv{i}.X_AW = sequenceCellActivityByOutcome(traces_cv{i},avgTraces_cv{i},normAvgTraces_cv{i},firingField_cv{i}.X_AW,find(passed_cv{i}.AW.X==1),'AW','X','A',prop,p);
    sequenceCellActivity_cv{i}.X_FF = sequenceCellActivityByOutcome(traces_cv{i},avgTraces_cv{i},normAvgTraces_cv{i},firingField_cv{i}.X_AW,find(passed_cv{i}.AW.X==1),'FF','X','A',prop,p);
end


%% Figure





