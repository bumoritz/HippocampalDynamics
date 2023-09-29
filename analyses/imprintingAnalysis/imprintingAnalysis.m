function [impr] = imprintingAnalysis(info,iscell,ops,p,path,task,trg,sync_beh,act)
% p = get_p; sync_beh = paq_beh;

%% Preparations

disp('--- Preparations')

% pre-process activity measure
[prop,~,nft_binned] = preprocessActivityMeasure(act,p.impr,p,sync_beh,iscell);

% create trials, traces, avgTraces and normAvgTraces structs
if ops.impr.do_allTrials
    trials_all = createTrialsStruct(task,1:prop.numTrials_incl,true);
    [traces_all,~,~] = createTracesStructs(nft_binned,trials_all);
end
if ops.impr.do_100t
    trials_100t = {};
    if ~ops.tng.do_100t
        temp = 1;
    else
        temp = floor(prop.numTrials/100);
    end
    for i=1:temp
        trials_100t{i} = createTrialsStruct(task,(i-1)*100+1:i*100,true);
        [traces_100t{i},~,~] = createTracesStructs(nft_binned,trials_100t{i});
    end
end


%% Store trg information

impr.trg.sequenceClusters = trg.sequenceClusters;
impr.trg.idcs_targetedCells = trg.idcs_targetedCells;
impr.trg.idcs_targetedCells_A_raw = trg.idcs_targetedCells(:,trg.sequenceClusters(:,1));
impr.trg.idcs_targetedCells_X_raw = trg.idcs_targetedCells(:,trg.sequenceClusters(:,2));
impr.trg.idcs_targetedCells_A = rmmissing(impr.trg.idcs_targetedCells_A_raw(:));
impr.trg.idcs_targetedCells_X = rmmissing(impr.trg.idcs_targetedCells_X_raw(:));
impr.trg.idcs_nonTargetedCells = setdiff(setdiff(find(prop.iscell==1),impr.trg.idcs_targetedCells_A),impr.trg.idcs_targetedCells_X);
temp = impr.trg.idcs_targetedCells_A_raw;
temp2 = repmat(1:size(impr.trg.idcs_targetedCells_A_raw,2),size(impr.trg.idcs_targetedCells_A_raw,1),1);
impr.trg.clusters_targetedCells_A = temp2(~isnan(temp));
temp = impr.trg.idcs_targetedCells_X_raw;
temp2 = repmat(1:size(impr.trg.idcs_targetedCells_X_raw,2),size(impr.trg.idcs_targetedCells_X_raw,1),1);
impr.trg.clusters_targetedCells_X = temp2(~isnan(temp));


%% Core _all

% so far this is all during analysis window. What about during window of specific stim time?

if ops.impr.do_allTrials
    
    % get relevant traces
    these_traces = {};
    these_traces.Aseq_Atrials_catch = traces_all.A_var0(impr.trg.idcs_targetedCells_A,p.general.bins_analysisWindow,:);
    these_traces.Aseq_Atrials_stim = traces_all.A_var1(impr.trg.idcs_targetedCells_A,p.general.bins_analysisWindow,:);
    these_traces.Aseq_Xtrials_catch = traces_all.X_var0(impr.trg.idcs_targetedCells_A,p.general.bins_analysisWindow,:);
    these_traces.Aseq_Xtrials_stim = traces_all.X_var1(impr.trg.idcs_targetedCells_A,p.general.bins_analysisWindow,:);
    these_traces.Xseq_Atrials_catch = traces_all.A_var0(impr.trg.idcs_targetedCells_X,p.general.bins_analysisWindow,:);
    these_traces.Xseq_Atrials_stim = traces_all.A_var1(impr.trg.idcs_targetedCells_X,p.general.bins_analysisWindow,:);
    these_traces.Xseq_Xtrials_catch = traces_all.X_var0(impr.trg.idcs_targetedCells_X,p.general.bins_analysisWindow,:);
    these_traces.Xseq_Xtrials_stim = traces_all.X_var1(impr.trg.idcs_targetedCells_X,p.general.bins_analysisWindow,:);
    these_traces.none_Atrials_catch = traces_all.A_var0(impr.trg.idcs_nonTargetedCells,p.general.bins_analysisWindow,:);
    these_traces.none_Atrials_stim = traces_all.A_var1(impr.trg.idcs_nonTargetedCells,p.general.bins_analysisWindow,:);
    these_traces.none_Xtrials_catch = traces_all.X_var0(impr.trg.idcs_nonTargetedCells,p.general.bins_analysisWindow,:);
    these_traces.none_Xtrials_stim = traces_all.X_var1(impr.trg.idcs_nonTargetedCells,p.general.bins_analysisWindow,:);

    % get activity
    temp = fields(these_traces);
    for i=1:length(temp)
        %impr.act.(temp{i}) = squeeze(nanmean(these_traces.(temp{i}),2));
        %impr.avgAct.(temp{i}) = nanmean(impr.act.(temp{i}),2);
        impr.avgAct.(temp{i}) = nanmean(squeeze(nanmean(these_traces.(temp{i}),2)),2);
    end

    % get pair-wise correlations
    temp = fields(these_traces);
    for i=1:length(temp)
        %impr.corr.(temp{i}) = getPairwiseCorrelations(these_traces.(temp{i}));
        %impr.avgCorr.(temp{i}) = nanmean(impr.corr.(temp{i}),3);
        impr.avgCorr.(temp{i}) = nanmean(getPairwiseCorrelations(these_traces.(temp{i})),3);
    end
    
    % get average pairwise correlations during stims
    temp = fields(these_traces);
    for i=1:length(temp)
        impr.avgCorr_pw.(temp{i}) = nan(1,size(impr.trg.sequenceClusters,1));
        for j=1:size(impr.trg.sequenceClusters,1)
            
            % average pairwise correlations between neurons that were stimulated in the same cluster
            if strcmp(temp{i}(1:4),'Aseq')
                these_idcs = find(impr.trg.clusters_targetedCells_A==j);
            elseif strcmp(temp{i}(1:4),'Xseq')
                these_idcs = find(impr.trg.clusters_targetedCells_X==j);
            elseif strcmp(temp{i}(1:4),'none')
                these_idcs = 1:length(impr.trg.idcs_nonTargetedCells);
            end
            impr.avgCorr_pw_raw.(temp{i}){j} = impr.avgCorr.(temp{i})(these_idcs,these_idcs);
            temp2 = impr.avgCorr_pw_raw.(temp{i}){j};
            impr.avgCorr_pw.(temp{i})(j) = nanmean(temp2(~eye(size(temp2))));
        end
    end
    
end


%% Core _100t

if ops.impr.do_100t
    impr_100t = {};
    if ~ops.tng.do_100t
        niter = 1;
    else
        niter = floor(prop.numTrials/100);
    end
    for k=1:niter
        impr_100t{k}.trg = impr.trg;

        % get relevant traces
        these_traces = {};
        these_traces.Aseq_Atrials_catch = traces_100t{k}.A_var0(impr_100t{k}.trg.idcs_targetedCells_A,p.general.bins_analysisWindow,:);
        these_traces.Aseq_Atrials_stim = traces_100t{k}.A_var1(impr_100t{k}.trg.idcs_targetedCells_A,p.general.bins_analysisWindow,:);
        these_traces.Aseq_Xtrials_catch = traces_100t{k}.X_var0(impr_100t{k}.trg.idcs_targetedCells_A,p.general.bins_analysisWindow,:);
        these_traces.Aseq_Xtrials_stim = traces_100t{k}.X_var1(impr_100t{k}.trg.idcs_targetedCells_A,p.general.bins_analysisWindow,:);
        these_traces.Xseq_Atrials_catch = traces_100t{k}.A_var0(impr_100t{k}.trg.idcs_targetedCells_X,p.general.bins_analysisWindow,:);
        these_traces.Xseq_Atrials_stim = traces_100t{k}.A_var1(impr_100t{k}.trg.idcs_targetedCells_X,p.general.bins_analysisWindow,:);
        these_traces.Xseq_Xtrials_catch = traces_100t{k}.X_var0(impr_100t{k}.trg.idcs_targetedCells_X,p.general.bins_analysisWindow,:);
        these_traces.Xseq_Xtrials_stim = traces_100t{k}.X_var1(impr_100t{k}.trg.idcs_targetedCells_X,p.general.bins_analysisWindow,:);
        these_traces.none_Atrials_catch = traces_100t{k}.A_var0(impr_100t{k}.trg.idcs_nonTargetedCells,p.general.bins_analysisWindow,:);
        these_traces.none_Atrials_stim = traces_100t{k}.A_var1(impr_100t{k}.trg.idcs_nonTargetedCells,p.general.bins_analysisWindow,:);
        these_traces.none_Xtrials_catch = traces_100t{k}.X_var0(impr_100t{k}.trg.idcs_nonTargetedCells,p.general.bins_analysisWindow,:);
        these_traces.none_Xtrials_stim = traces_100t{k}.X_var1(impr_100t{k}.trg.idcs_nonTargetedCells,p.general.bins_analysisWindow,:);

        % get activity
        temp = fields(these_traces);
        for i=1:length(temp)
            %impr_100t{k}.act.(temp{i}) = squeeze(nanmean(these_traces.(temp{i}),2));
            %impr_100t{k}.avgAct.(temp{i}) = nanmean(impr_100t{k}.act.(temp{i}),2);
            impr_100t{k}.avgAct.(temp{i}) = nanmean(squeeze(nanmean(these_traces.(temp{i}),2)),2);
        end
        
        % get pair-wise correlations
        temp = fields(these_traces);
        for i=1:length(temp)
            %impr_100t{k}.corr.(temp{i}) = getPairwiseCorrelations(these_traces.(temp{i}));
            %impr_100t{k}.avgCorr.(temp{i}) = nanmean(impr_100t{k}.corr.(temp{i}),3);
            impr_100t{k}.avgCorr.(temp{i}) = nanmean(getPairwiseCorrelations(these_traces.(temp{i})),3);
        end

        % get average pairwise correlations during stims
        temp = fields(these_traces);
        for i=1:length(temp)
            impr_100t{k}.avgCorr_pw.(temp{i}) = nan(1,size(impr_100t{k}.trg.sequenceClusters,1));
            for j=1:size(impr_100t{k}.trg.sequenceClusters,1)

                % average pairwise correlations between neurons that were stimulated in the same cluster
                if strcmp(temp{i}(1:4),'Aseq')
                    these_idcs = find(impr_100t{k}.trg.clusters_targetedCells_A==j);
                elseif strcmp(temp{i}(1:4),'Xseq')
                    these_idcs = find(impr_100t{k}.trg.clusters_targetedCells_X==j);
                elseif strcmp(temp{i}(1:4),'none')
                    these_idcs = 1:length(impr_100t{k}.trg.idcs_nonTargetedCells);
                end
                impr_100t{k}.avgCorr_pw_raw.(temp{i}){j} = impr_100t{k}.avgCorr.(temp{i})(these_idcs,these_idcs);
                temp2 = impr_100t{k}.avgCorr_pw_raw.(temp{i}){j};
                impr_100t{k}.avgCorr_pw.(temp{i})(j) = nanmean(temp2(~eye(size(temp2))));
            end
        end
    end
end


%% Snake plots

% preparations
if ~exist([path.filepart_outX,'plots/SnakePlots/'],'dir')
    mkdir([path.filepart_outX,'plots/SnakePlots/']);
end
if strcmp(p.impr.activityMeasure,"casc_50c_beh")
    this_activity_measure = 'casc_50c';
elseif strcmp(p.impr.activityMeasure,"casc_50g_beh")
    this_activity_measure = 'casc_50g';
elseif strcmp(p.impr.activityMeasure,"casc_100c_beh")
    this_activity_measure = 'casc_100c';
elseif strcmp(p.impr.activityMeasure,"casc_100g_beh")
    this_activity_measure = 'casc_100g';
elseif strcmp(p.impr.activityMeasure,"casc_200g_beh")
    this_activity_measure = 'casc_200g';
elseif strcmp(p.impr.activityMeasure,"spks_beh")
    this_activity_measure = 'spks';
elseif strcmp(p.impr.activityMeasure,"spksn_beh")
    this_activity_measure = 'spksn';
elseif strcmp(p.impr.activityMeasure,"spksnn_beh")
    this_activity_measure = 'spksnn';
elseif strcmp(p.impr.activityMeasure,"dFF_beh")
    this_activity_measure = 'dFF';
elseif strcmp(p.impr.activityMeasure,"dFFn_beh")
    this_activity_measure = 'dFFn';
elseif strcmp(p.impr.activityMeasure,"dFFnn_beh")
    this_activity_measure = 'dFFnn';
elseif strcmp(p.impr.activityMeasure,"F_beh")
    this_activity_measure = 'F';
end

if ops.impr.do_allTrials

    % AXtargets, all trials, 8tile, stim_catch
    plt = struct(); plt.split = 'stim_catch';
    temp1 = trg.idcs_targetedCells(:,trg.sequenceClusters(:,1)); 
    temp2 = trg.idcs_targetedCells(:,trg.sequenceClusters(:,2)); 
    F = snakePlots(traces_all,[],{rmmissing(temp1(:)),rmmissing(temp2(:))},[],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.impr.smoothingSd_preBinning),'-bin',num2str(p.impr.binSize),', ',...
        'AXtargets, all trials'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AXtargets_allTrials_8tile_stimcatch.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AXtargets_allTrials_8tile_stimcatch.png']);
% 
%     % AXtargets, all trials, 12tile, stimOstimMcatch, same, 14
%     plt = struct(); plt.split = 'stimO_stimM_catch'; plt.split3_numTrials = 'same';
%     F = snakePlots(traces_all,{rmmissing(temp1(:)),rmmissing(temp2(:))},[],[1;4],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.impr.smoothingSd_preBinning),'-bin',num2str(p.impr.binSize),', ',...
%         'AXtargets, all trials'])
%     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AXtargets_allTrials_12tile_stimOstimMcatch_same_14.fig']);
%     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AXtargets_allTrials_12tile_stimOstimMcatch_same_14.png']);
% 
%     % AXtargets, all trials, 12tile, stimOstimMcatch, same, 14
%     plt = struct(); plt.split = 'stimO_stimM_catch'; plt.split3_numTrials = 'same'; plt.normalisationRef = 'plotwise';
%     F = snakePlots(traces_all,{rmmissing(temp1(:)),rmmissing(temp2(:))},[],[1;4],[],p,info,plt);
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%         this_activity_measure,'-smo',num2str(p.impr.smoothingSd_preBinning),'-bin',num2str(p.impr.binSize),', ',...
%         'AXtargets, all trials, normalised within'])
% %     savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AXtargets_allTrials_12tile_stimOstimMcatch_same_14.fig']);
% %     saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AXtargets_allTrials_12tile_stimOstimMcatch_same_14.png']);

    if ~ops.impr.skip_boringSnakePlots
        
        % AXtargets, all trials, 12tile, stimOstimMcatch, same, 41
        plt = struct(); plt.split = 'stimO_stimM_catch'; plt.split3_numTrials = 'same';
        F = snakePlots(traces_all,{rmmissing(temp1(:)),rmmissing(temp2(:))},[],[4;1],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.impr.smoothingSd_preBinning),'-bin',num2str(p.impr.binSize),', ',...
            'AXtargets, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AXtargets_allTrials_12tile_stimOstimMcatch_same_41.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AXtargets_allTrials_12tile_stimOstimMcatch_same_41.png']);

        % AXtargets, all trials, 12tile, stimOstimMcatch, same, 36
        plt = struct(); plt.split = 'stimO_stimM_catch'; plt.split3_numTrials = 'same';
        F = snakePlots(traces_all,{rmmissing(temp1(:)),rmmissing(temp2(:))},[],[3;6],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.impr.smoothingSd_preBinning),'-bin',num2str(p.impr.binSize),', ',...
            'AXtargets, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AXtargets_allTrials_12tile_stimOstimMcatch_same_36.fig']);
            saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AXtargets_allTrials_12tile_stimOstimMcatch_same_36.png']);

        % AXtargets, all trials, 12tile, stimOstimMcatch, same, 63
        plt = struct(); plt.split = 'stimO_stimM_catch'; plt.split3_numTrials = 'same';
        F = snakePlots(traces_all,{rmmissing(temp1(:)),rmmissing(temp2(:))},[],[6;3],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.impr.smoothingSd_preBinning),'-bin',num2str(p.impr.binSize),', ',...
            'AXtargets, all trials'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AXtargets_allTrials_12tile_stimOstimMcatch_same_63.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AXtargets_allTrials_12tile_stimOstimMcatch_same_63.png']);
    end
end

if ops.impr.do_60t && ops.impr.do_allTrials
    
    % Atargets, 60t, AinA
    plt = struct(); plt.split = 'none'; plt.trialTypes = "A"; plt.split = 'stim_catch'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
    temp1 = trg.idcs_targetedCells(:,trg.sequenceClusters(:,1));
    F = snakePlots(traces_60t,[],{rmmissing(temp1(:))},[],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.impr.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'Atargets, 60t, AinA'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Atargets_60t_AinA.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Atargets_60t_AinA.png']);

    % Xtargets, 60t, XinX
    plt = struct(); plt.split = 'none'; plt.trialTypes = "X"; plt.split = 'stim_catch'; plt.normalisationRef = 'uniform'; plt.sequences = "X"; 
    temp1 = trg.idcs_targetedCells(:,trg.sequenceClusters(:,2));
    F = snakePlots(traces_60t,[],{rmmissing(temp1(:))},[],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.impr.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
        'Xtargets, 60t, XinX'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Xtargets_60t_XinX.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Xtargets_60t_XinX.png']);

    if ~ops.impr.skip_boringSnakePlots
    
        % Atargets, 60t, AinX
        plt = struct(); plt.split = 'none'; plt.trialTypes = "X"; plt.split = 'stim_catch'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
        temp1 = trg.idcs_targetedCells(:,trg.sequenceClusters(:,1));
        F = snakePlots(traces_60t,[],{rmmissing(temp1(:))},[],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.impr.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'Atargets, 60t, AinX'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Atargets_60t_AinX.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Atargets_60t_AinX.png']);
        
        % Xtargets, 60t, XinA
        plt = struct(); plt.split = 'none'; plt.trialTypes = "A"; plt.split = 'stim_catch'; plt.normalisationRef = 'uniform'; plt.sequences = "X"; 
        temp1 = trg.idcs_targetedCells(:,trg.sequenceClusters(:,2));
        F = snakePlots(traces_60t,[],{rmmissing(temp1(:))},[],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.impr.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'Xtargets, 60t, XinA'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Xtargets_60t_XinA.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Xtargets_60t_XinA.png']);
    
        % Atargets, 60t, AinA, n4
        this_block = 4;
        plt = struct(); plt.split = 'none'; plt.trialTypes = "A"; plt.split = 'stim_catch'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
        temp1 = trg.idcs_targetedCells(:,trg.sequenceClusters(:,1));
        F = snakePlots(traces_60t,{rmmissing(temp1(:))},[],[NaN,this_block],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.impr.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'Atargets, 60t, AinA'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Atargets_60t_AinA_n',num2str(this_block),'.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Atargets_60t_AinA_n',num2str(this_block),'.png']);

        % Atargets, 60t, AinX, n4
        this_block = 4;
        plt = struct(); plt.split = 'none'; plt.trialTypes = "X"; plt.split = 'stim_catch'; plt.normalisationRef = 'uniform'; plt.sequences = "A"; 
        temp1 = trg.idcs_targetedCells(:,trg.sequenceClusters(:,1));
        F = snakePlots(traces_60t,{rmmissing(temp1(:))},[],[NaN,this_block],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.impr.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'Atargets, 60t, AinX'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Atargets_60t_AinX_n',num2str(this_block),'.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Atargets_60t_AinX_n',num2str(this_block),'.png']);

        % Xtargets, 60t, XinX, n4
        this_block = 4;
        plt = struct(); plt.split = 'none'; plt.trialTypes = "X"; plt.split = 'stim_catch'; plt.normalisationRef = 'uniform'; plt.sequences = "X"; 
        temp1 = trg.idcs_targetedCells(:,trg.sequenceClusters(:,2));
        F = snakePlots(traces_60t,{rmmissing(temp1(:))},[],[NaN,this_block],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.impr.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'Xtargets, 60t, XinX'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Xtargets_60t_XinX_n',num2str(this_block),'.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Xtargets_60t_XinX_n',num2str(this_block),'.png']);

        % Xtargets, 60t, XinA, n4
        this_block = 4;
        plt = struct(); plt.split = 'none'; plt.trialTypes = "A"; plt.split = 'stim_catch'; plt.normalisationRef = 'uniform'; plt.sequences = "X"; 
        temp1 = trg.idcs_targetedCells(:,trg.sequenceClusters(:,2));
        F = snakePlots(traces_60t,{rmmissing(temp1(:))},[],[NaN,this_block],[],p,info,plt);
        suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
            this_activity_measure,'-smo',num2str(p.impr.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
            'Xtargets, 60t, XinA'])
        savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Xtargets_60t_XinA_n',num2str(this_block),'.fig']);
        saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_Xtargets_60t_XinA_n',num2str(this_block),'.png']);
    end
end

disp(['--- Saved snake plots to ',path.filepart_outX,'plots/SnakePlots.'])


%% Imprinting analysis figure

% F = figure;
% savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','imprintingAnalysis.fig']);
% saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','imprintingAnalysis.png']);
% disp(['--- Saved imprinting analysis figure to ',path.filepart_outX,'plots.'])


%% Save and return

if ops.impr.do_allTrials
    impr_all = impr;
    save([path.filepart_outX,'impr_all.mat'],'impr_all','-v7.3');
    disp(['--- Saved impr_all file as ',[path.filepart_outX,'impr_all.mat'],'.'])
end
if ops.impr.do_100t
    save([path.filepart_outX,'impr_100t.mat'],'impr_100t','-v7.3');
    disp(['--- Saved impr_100t file as ',[path.filepart_outX,'impr_100t.mat'],'.'])
end
if ops.impr.do_allTrials
    impr = impr_all;
elseif ops.impr.do_60t || ops.impr.do_60t_onlyFirst
    impr = impr_60t{1};
else
    impr = impr_100t{1};
end


if ops.close_figures
    close all;
end
%end