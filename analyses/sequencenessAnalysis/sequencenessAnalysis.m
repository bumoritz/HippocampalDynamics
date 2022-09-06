function [sqn] = sequencenessAnalysis(info,iscell,ops,p,path,task,s2p_meta,sync_beh,act,tng_all,tng_100t,nem_all_cmpr)
%act = spks_beh; p = get_p; sync_beh = paq_beh; 

%% Preparations

disp('--- Preparations')

% pre-process activity measure
[prop,~,nft_binned] = preprocessActivityMeasure(act,p.sqn,p,sync_beh,iscell);

% create trials, traces, avgTraces and normAvgTraces structs
if ops.sqn.do_allTrials
    trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
    trials_all = balanceTrialSelelection(trials_all,info,task);
    [traces_all,avgTraces_all,~] = createTracesStructs(nft_binned,trials_all);
end
if ops.sqn.do_100t
    trials_100t = {}; traces_100t = {}; avgTraces_100t = {};
    for i=1:floor(prop.numTrials/100)
        trials_100t{i} = createTrialsStruct(task,(i-1)*100+1:i*100);
        [traces_100t{i},avgTraces_100t{i},~] = createTracesStructs(nft_binned,trials_100t{i});
    end
end


%% Template correlation

disp('--- Template correlation')

if ops.sqn.do_allTrials
    tempCorr_all = templateCorrelation(prop,traces_all,avgTraces_all,tng_all,p);
end
if ops.sqn.do_100t
    tempCorr_100t = {};
    for i=1:floor(prop.numTrials/100)
        tempCorr_100t{i} = templateCorrelation(prop,traces_100t{i},avgTraces_100t{i},tng_100t{i},p);
    end
end


%% Population vector similarity

disp('--- Population vector similarity')

if ops.sqn.do_allTrials
    povSim_all = populationVectorSimilarity(prop,traces_all,avgTraces_all,tng_all,nem_all_cmpr,p);
end
if ops.sqn.do_100t
    povSim_100t = {};
    for i=1:floor(prop.numTrials/100)
        povSim_100t{i} = populationVectorSimilarity(prop,traces_100t{i},avgTraces_100t{i},tng_100t{i},p);
    end
end


%% Make figure

% preparations
if ~exist([path.filepart_outX,'plots/SequencenessPlots/'],'dir')
    mkdir([path.filepart_outX,'plots/SequencenessPlots/']);
end

if ops.sqn.do_allTrials

    F = templateCorrelationFigure(tempCorr_all,'maxBinCorr','Spearman',p);
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType]);
    savefig(F,[path.filepart_outX,'plots/SequencenessPlots/',info.animal,'_',info.date,'_','templateCorrelationFigure.fig']);
    saveas(F,[path.filepart_outX,'plots/SequencenessPlots/',info.animal,'_',info.date,'_','templateCorrelationFigure.png']);
    
    disp(['--- Saved sequenceness plots to ',path.filepart_outX,'plots/SequencenessPlots.'])
end


%% Save and return

if ops.sqn.do_allTrials
    sqn_all.tempCorr = tempCorr_all;
    sqn_all.povSim = povSim_all;
    sqn_all = orderfields(sqn_all);
    save([path.filepart_out,'sqn_all.mat'],'sqn_all','-v7.3');
    disp(['--- Saved sqn_all file as ',[path.filepart_out,'sqn_all.mat'],'.'])
end
if ops.sqn.do_100t
    sqn_100t = {};
    for i=1:floor(prop.numTrials/100)
        sqn_100t{i}.tempCorr = tempCorr_100t{i};
        sqn_100t{i}.povSim = povSim_100t{i};
        sqn_100t{i} = orderfields(sqn_100t{i});
    end
    save([path.filepart_outX,'sqn_100t.mat'],'sqn_100t','-v7.3');
    disp(['--- Saved sqn_100t file as ',[path.filepart_outX,'sqn_100t.mat'],'.'])
end
if ops.sqn.do_allTrials
    sqn = sqn_all;
else
    sqn = sqn_100t{1};
end

if ops.close_figures
    close all;
end
end







