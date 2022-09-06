function decodingAnalysis(info,iscell,ops,p,path,task,paq_beh,act)
%act = spks_beh; p = get_p;

%% Preparations

disp('--- Preparations')

% pre-process activity measure
[prop,~,nft_binned] = preprocessActivityMeasure(act,p.dec,p,paq_beh,iscell);

% create trials, traces, avgTraces and normAvgTraces structs
if ops.dec.do_allTrials
    trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
    trials_all = balanceTrialSelelection(trials_all,info,task);
end


%% Core - all trials

if ops.dec.do_allTrials
    dec_all = prepareDecoding(trials_all,p);
    [dec_all,~] = decodingCore(dec_all,nft_binned,prop,p);
    dec_all.analysis.trainingSet = analyseDecoding(dec_all,'trainingSet',1:dec_all.cv.numTrials_trainingSet,p);

    dec_all.analysis.testSet_corrM = analyseDecoding(dec_all,'completeSet',dec_all.cv.trialsForTestSet_corrM.combined,p);
    dec_all.analysis.testSet_incorrM = analyseDecoding(dec_all,'completeSet',dec_all.cv.trialsForTestSet_incorrM.combined,p);
    dec_all.analysis.completeSet = analyseDecoding(dec_all,'completeSet',1:dec_all.cv.numTrials_completeSet,p);
end


%% Save

if ops.dec.do_allTrials
    save([path.filepart_outX,'dec_all.mat'],'dec_all','-v7.3');
    disp(['--- Saved dec_all file as ',[path.filepart_outX,'dec_all.mat'],'.'])
end


%% Return

if ops.close_figures
    close all;
end
end



