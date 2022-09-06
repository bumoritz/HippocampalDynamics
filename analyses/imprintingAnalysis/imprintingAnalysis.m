%function [impr] = imprintingAnalysis(info,iscell,ops,p,path,task,s2p_meta,sync_beh,act)
act = spks_beh; p = get_p; sync_beh = thor_beh;

%% Preparations

disp('--- Preparations')

% pre-process activity measure
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.impr,p,sync_beh,iscell);

% create trials, traces, avgTraces and normAvgTraces structs
if ops.impr.do_allTrials
    trials_all = createTrialsStruct(task,1:prop.numTrials,true);
    [traces_all,avgTraces_all,normAvgTraces_all] = createTracesStructs(nft_binned,trials_all);
end
if ops.impr.do_60t || ops.impr.do_60t_onlyFirst
    trials_60t = {}; traces_60t = {}; avgTraces_60t = {}; normAvgTraces_60t = {};
    if ~ops.impr.do_60t
        temp = 1;
    else
        temp = floor(prop.numTrials/60);
    end
    for i=1:temp
        trials_60t{i} = createTrialsStruct(task,(i-1)*60+1:i*60,true);
        [traces_60t{i},avgTraces_60t{i},normAvgTraces_60t{i}] = createTracesStructs(nft_binned,trials_60t{i});
    end
end
if ops.impr.do_100t || ops.impr.do_100t_onlyFirst
    trials_100t = {}; traces_100t = {}; avgTraces_100t = {}; normAvgTraces_100t = {};
    if ~ops.impr.do_100t
        temp = 1;
    else
        temp = floor(prop.numTrials/100);
    end
    for i=1:temp
        trials_100t{i} = createTrialsStruct(task,(i-1)*100+1:i*100,true);
        [traces_100t{i},avgTraces_100t{i},normAvgTraces_100t{i}] = createTracesStructs(nft_binned,trials_100t{i});
    end
end


%%













%% temp

%%% start with running Preparations section of sequence cell analysis
% trg = trg_nonrigid;
% 
% 
% %% kind of from snake plot utility
% 
% data_A = sca.traces_A;
% data_X = sca.traces_X;
% 
% avgTraces_A_1 = nanmean(data_A(:,:, find(task.var(find(task.odour1=="A"))) ),3);
% avgTraces_A_2 = nanmean(data_A(:,:, find(~task.var(find(task.odour1=="A"))) ),3);
% avgTraces_X_1 = nanmean(data_X(:,:, find(task.var(find(task.odour1=="X"))) ),3);
% avgTraces_X_2 = nanmean(data_X(:,:, find(~task.var(find(task.odour1=="X"))) ),3);
% 
% % 1 is stim, 2 is catch
% normAvgTraces_A_1_Aseq = (avgTraces_A_1-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
% normAvgTraces_A_2_Aseq = (avgTraces_A_2-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
% normAvgTraces_X_1_Aseq = (avgTraces_X_1-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
% normAvgTraces_X_2_Aseq = (avgTraces_X_2-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
% normAvgTraces_A_1_Xseq = (avgTraces_A_1-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2));
% normAvgTraces_A_2_Xseq = (avgTraces_A_2-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2));
% normAvgTraces_X_1_Xseq = (avgTraces_X_1-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2));
% normAvgTraces_X_2_Xseq = (avgTraces_X_2-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2));
% 
% temp = trg.idcs_targetedCells(:,trg.sequenceClusters(:,1));
% idcs_sorted_Aseq = rmmissing(temp(:));
% numSel_A = NaN;
% temp = trg.idcs_targetedCells(:,trg.sequenceClusters(:,2));
% idcs_sorted_Xseq = rmmissing(temp(:));
% numSel_X = NaN;
% 
% AstimAseq = normAvgTraces_A_1_Aseq(idcs_sorted_Aseq,:);
% AcatchAseq = normAvgTraces_A_2_Aseq(idcs_sorted_Aseq,:);
% XstimAseq = normAvgTraces_X_1_Aseq(idcs_sorted_Aseq,:);
% XcatchAseq = normAvgTraces_X_2_Aseq(idcs_sorted_Aseq,:);
% AstimXseq = normAvgTraces_A_1_Xseq(idcs_sorted_Xseq,:);
% AcatchXseq = normAvgTraces_A_2_Xseq(idcs_sorted_Xseq,:);
% XstimXseq = normAvgTraces_X_1_Xseq(idcs_sorted_Xseq,:);
% XcatchXseq = normAvgTraces_X_2_Xseq(idcs_sorted_Xseq,:);
% 
% 
% %% A seq correlation 
% 
% corr_AstimAseq_AcatchAseq = nan(length(idcs_sorted_Aseq),1);
% for i=1:length(idcs_sorted_Aseq)
%     corr_AstimAseq_AcatchAseq(i) = corr(AstimAseq(i,:)',AcatchAseq(i,:)');
% end
% 
% corr_AstimAseq_XcatchAseq = nan(length(idcs_sorted_Aseq),1);
% for i=1:length(idcs_sorted_Aseq)
%     corr_AstimAseq_XcatchAseq(i) = corr(AstimAseq(i,:)',XcatchAseq(i,:)');
% end
% 
% figure;
% v = violinplot([corr_AstimAseq_AcatchAseq,corr_AstimAseq_XcatchAseq],{'corr(A-stim,A-nostim)','corr(A-stim,X-nostim)'});
% title('SeqA cell activity correlation between stim trials and no-stim trials')
% 
% 
% %% X seq correlation 
% 
% corr_XstimXseq_XcatchXseq = nan(length(idcs_sorted_Xseq),1);
% for i=1:length(idcs_sorted_Xseq)
%     corr_XstimXseq_XcatchXseq(i) = corr(XstimXseq(i,:)',XcatchXseq(i,:)');
% end
% 
% corr_XstimXseq_AcatchXseq = nan(length(idcs_sorted_Xseq),1);
% for i=1:length(idcs_sorted_Xseq)
%     corr_XstimXseq_AcatchXseq(i) = corr(XstimXseq(i,:)',AcatchXseq(i,:)');
% end
% 
% figure;
% v = violinplot([corr_XstimXseq_XcatchXseq,corr_XstimXseq_AcatchXseq],{'corr(X-stim,X-nostim)','corr(X-stim,A-nostim)'});
% title('SeqX cell activity correlation between stim trials and no-stim trials')
% 
% 
% %% --- %%% --- %%%
% 
% %% Pop vector correlation - step size of 5
% 
% figure;
% 
% binSize = 3;
% 
% temp = movmean(avgTraces_A_1(idcs_sorted_Aseq,sca.prop.analysisWindow),binSize,2,'omitnan');
% templateA = temp(:,floor(binSize/2)+1:binSize:end,:);
% 
% catch_trials = find(~task.var(find(task.odour1=="A")));
% catchDataA = {};
% for i=1:16
%     temp = movmean(   nanmean(data_A(idcs_sorted_Aseq,sca.prop.analysisWindow, catch_trials((i-1)*5+1:i*5) ),3)   ,binSize,2,'omitnan');
%     catchDataA{i} = temp(:,floor(binSize/2)+1:binSize:end,:);
% end
% 
% corr_templateA_catchDataA = nan(size(templateA,2),16);
% for i=1:16
%     for j=1:size(templateA,2)
%         corr_templateA_catchDataA(j,i) = corr(templateA(:,j),catchDataA{i}(:,j));
%     end
% end
% 
% subplot(1,2,1)
% imagesc(corr_templateA_catchDataA)
% colorbar
% xlabel('Catch trial bin (1 bin = 5 catch trials)')
% ylabel('Time bin (1 bin = 300 ms)')
% title('PV corr for A: catch trial bin with template from all stim trials')
% 
% 
% binSize = 3;
% temp = movmean(avgTraces_X_1(idcs_sorted_Xseq,sca.prop.analysisWindow),binSize,2,'omitnan');
% templateX = temp(:,floor(binSize/2)+1:binSize:end,:);
% 
% catch_trials = find(~task.var(find(task.odour1=="X")));
% catchDataX = {};
% for i=1:16
%     temp = movmean(   nanmean(data_X(idcs_sorted_Xseq,sca.prop.analysisWindow, catch_trials((i-1)*5+1:i*5) ),3)   ,binSize,2,'omitnan');
%     catchDataX{i} = temp(:,floor(binSize/2)+1:binSize:end,:);
% end
% 
% corr_templateX_catchDataX = nan(size(templateX,2),16);
% for i=1:16
%     for j=1:size(templateX,2)
%         corr_templateX_catchDataX(j,i) = corr(templateX(:,j),catchDataX{i}(:,j));
%     end
% end
% 
% subplot(1,2,2)
% imagesc(corr_templateX_catchDataX)
% colorbar
% xlabel('Catch trial bin (1 bin = 5 catch trials)')
% ylabel('Time bin (1 bin = 300 ms)')
% title('PV corr for X: catch trial bin with template from all stim trials')
% 
% 
% %% Pop vector correlation - step size of 1
% 
% figure;
% 
% binSize = 3;
% 
% temp = movmean(avgTraces_A_1(idcs_sorted_Aseq,sca.prop.analysisWindow),binSize,2,'omitnan');
% templateA = temp(:,floor(binSize/2)+1:binSize:end,:);
% 
% catch_trials = find(~task.var(find(task.odour1=="A")));
% catchDataA = {};
% for i=1:size(data_A,3)/5
%     temp = movmean(   nanmean(data_A(idcs_sorted_Aseq,sca.prop.analysisWindow, catch_trials((i-1)*1+1:i*1) ),3)   ,binSize,2,'omitnan');
%     catchDataA{i} = temp(:,floor(binSize/2)+1:binSize:end,:);
% end
% 
% corr_templateA_catchDataA = nan(size(templateA,2),size(data_A,3)/5);
% for i=1:size(data_A,3)/5
%     for j=1:size(templateA,2)
%         corr_templateA_catchDataA(j,i) = corr(templateA(:,j),catchDataA{i}(:,j));
%     end
% end
% 
% subplot(1,2,1)
% imagesc(corr_templateA_catchDataA)
% colorbar
% xlabel('Catch trial')
% ylabel('Time bin (1 bin = 300 ms)')
% title('PV corr for A: catch trial bin with template from all stim trials')
% 
% 
% binSize = 3;
% temp = movmean(avgTraces_X_1(idcs_sorted_Xseq,sca.prop.analysisWindow),binSize,2,'omitnan');
% templateX = temp(:,floor(binSize/2)+1:binSize:end,:);
% 
% catch_trials = find(~task.var(find(task.odour1=="X")));
% catchDataX = {};
% for i=1:size(data_X,3)/5
%     temp = movmean(   nanmean(data_X(idcs_sorted_Xseq,sca.prop.analysisWindow, catch_trials((i-1)*1+1:i*1) ),3)   ,binSize,2,'omitnan');
%     catchDataX{i} = temp(:,floor(binSize/2)+1:binSize:end,:);
% end
% 
% corr_templateX_catchDataX = nan(size(templateX,2),size(data_X,3)/5);
% for i=1:size(data_X,3)/5
%     for j=1:size(templateX,2)
%         corr_templateX_catchDataX(j,i) = corr(templateX(:,j),catchDataX{i}(:,j));
%     end
% end
% 
% subplot(1,2,2)
% imagesc(corr_templateX_catchDataX)
% colorbar
% xlabel('Catch trial')
% ylabel('Time bin (1 bin = 300 ms)')
% title('PV corr for X: catch trial bin with template from all stim trials')
% %%
% 
% figure;
% hold on
% plot(nanmean(corr_templateA_catchDataA,1),'Color',p.col.odourA)
% plot(nanmean(corr_templateX_catchDataX,1),'Color',p.col.odourX)
% xlabel('Catch trial bin (1 bin = 5 catch trials)')
% ylabel('Correlation')
% title('PV corr: catch trial bin with template from corresponding stim trials')
% hold off
% 
% 
% %%
% temp = nchoosek(1:prop.numTrials/100,2)
% stability = ones(prop.numTrials/100,prop.numTrials/100);
% for i=1:length(temp)
%     temp3 = size(sca_100t{1,temp(i,1)}.avgTraces_A(find(s2p_meta.iscell(:,1)),:),2);
%     temp2 = nan(temp3,1);
%     for j=1:temp3
%         temp2(j) = corr(sca_100t{1,temp(i,1)}.avgTraces_A(find(s2p_meta.iscell(:,1)),j),sca_100t{1,temp(i,2)}.avgTraces_A(find(s2p_meta.iscell(:,1)),j),'Rows','complete');
%     end
%     stability(temp(i,1),temp(i,2)) = nanmean(temp2);
%     stability(temp(i,2),temp(i,1)) = nanmean(temp2);
%     i
% end
% 
% 
% 
% 
% 
% 





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
elseif strcmp(p.impr.activityMeasure,"dFF_beh")
    this_activity_measure = 'dFF';
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

    % AXtargets, all trials, 12tile, stimOstimMcatch, same, 14
    plt = struct(); plt.split = 'stimO_stimM_catch'; plt.split3_numTrials = 'same';
    F = snakePlots(traces_all,{rmmissing(temp1(:)),rmmissing(temp2(:))},[],[1;4],[],p,info,plt);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
        this_activity_measure,'-smo',num2str(p.impr.smoothingSd_preBinning),'-bin',num2str(p.impr.binSize),', ',...
        'AXtargets, all trials'])
    savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AXtargets_allTrials_12tile_stimOstimMcatch_same_14.fig']);
    saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_','snakePlot_AXtargets_allTrials_12tile_stimOstimMcatch_same_14.png']);

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

F = figure;
savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','imprintingAnalysis.fig']);
saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','imprintingAnalysis.png']);
disp(['--- Saved imprinting analysis figure to ',path.filepart_outX,'plots.'])


%% Save and return

if ops.impr.do_allTrials
    impr_all.trials = trials_all;
    impr_all.traces = traces_all;
    impr_all.avgTraces = avgTraces_all;
    impr_all.normAvgTraces = normAvgTraces_all;
    impr_all.prop = prop;
    impr_all.p = p;
    impr_all = orderfields(impr_all);
    save([path.filepart_out,'impr_all.mat'],'impr_all','-v7.3');
    disp(['--- Saved impr_all file as ',[path.filepart_out,'impr_all.mat'],'.'])
end
if ops.impr.do_60t || ops.impr.do_60t_onlyFirst
    impr_60t.trials = trials_60t;
    impr_60t.traces = traces_60t;
    impr_60t.avgTraces = avgTraces_60t;
    impr_60t.normAvgTraces = normAvgTraces_60t;
    impr_60t.prop = prop;
    impr_60t.p = p;
    impr_60t = orderfields(impr_60t);
    save([path.filepart_outX,'impr_60t.mat'],'impr_60t','-v7.3');
    disp(['--- Saved impr_60t file as ',[path.filepart_outX,'impr_60t.mat'],'.'])
end
if ops.impr.do_100t || ops.impr.do_100t_onlyFirst
    impr_100t.trials = trials_100t;
    impr_100t.traces = traces_100t;
    impr_100t.avgTraces = avgTraces_100t;
    impr_100t.normAvgTraces = normAvgTraces_100t;
    impr_100t.prop = prop;
    impr_100t.p = p;
    impr_100t = orderfields(impr_100t);
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