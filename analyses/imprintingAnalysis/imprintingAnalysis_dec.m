function [impr] = imprintingAnalysis_dec(info,iscell,ops,p,path,task,s2p_meta,sync_beh,act)
p = get_p; sync_beh = thor_beh;

%% Preparations

disp('--- Preparations')

% pre-process activity measure
[prop,~,nft_binned] = preprocessActivityMeasure(act,p.impr,p,sync_beh,iscell);

% create trials, traces, avgTraces and normAvgTraces structs
if ops.impr.do_allTrials
    trials_all = createTrialsStruct(task,1:prop.numTrials_incl,true);
    [traces_all,avgTraces_all,~] = createTracesStructs(nft_binned,trials_all);
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
    end
end


%% Core - all trials

% if ops.impr.do_allTrials        
%     dec_all = prepareDecoding_catch_seq(trials_all,p);
%     [dec_all,~] = decodingCore(dec_all,nft_binned,prop,p,true);
%     dec_all.analysis.trainingSet = analyseDecoding(dec_all,'trainingSet',1:dec_all.cv.numTrials_trainingSet,p);
%     dec_all.analysis.testSet = analyseDecoding(dec_all,'completeSet',dec_all.cv.trialsForTestSet.combined,p);
%     dec_all.analysis.completeSet = analyseDecoding(dec_all,'completeSet',1:dec_all.cv.numTrials_completeSet,p);
% end




