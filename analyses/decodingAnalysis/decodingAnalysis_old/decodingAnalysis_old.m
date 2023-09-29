function decodingAnalysis_old(info,iscell,ops,p,path,task,paq_beh,act)
%act = spksn_beh; p = get_p;

%% Preparations

disp('--- Preparations')

% pre-process activity measure
[prop,~,nft_binned] = preprocessActivityMeasure(act,p.dec,p,paq_beh,iscell);

% create trials, traces, avgTraces and normAvgTraces structs
if ops.dec.do_allTrials
    trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
    trials_all = balanceTrialSelelection(trials_all,info,task);
end
if ops.dec.do_100t
    trials_100t = {};
    if ~ops.tng.do_100t
        temp = 1;
    else
        temp = floor(prop.numTrials/100);
    end
    for i=1:temp
        trials_100t{i} = createTrialsStruct(task,(i-1)*100+1:i*100);
    end
end


%% Core - all trials
 
if ops.dec.do_allTrials
    
    p.dec.trainOnCompleteSet = false;
    
    dec_all = prepareDecoding(trials_all,p);
    [dec_all,~] = decodingCore(dec_all,nft_binned,prop,p);
    dec_all.analysis.trainingSet = analyseDecoding(dec_all,'trainingSet',1:dec_all.cv.numTrials_trainingSet,p);

    dec_all.analysis.testSet_corrM = analyseDecoding(dec_all,'completeSet',dec_all.cv.trialsForTestSet_corrM.combined,p);
    dec_all.analysis.testSet_incorrM = analyseDecoding(dec_all,'completeSet',dec_all.cv.trialsForTestSet_incorrM.combined,p);
    dec_all.analysis.completeSet = analyseDecoding(dec_all,'completeSet',1:dec_all.cv.numTrials_completeSet,p);
end


%% Core - block-wise

if ops.dec.do_100t
    dec_100t = {};
    if ~ops.tng.do_100t
        temp = 1;
    else
        temp = floor(prop.numTrials/100);
    end
    for i=1:temp
%         p.dec.trainOnCompleteSet = true;
%         dec_100t{i} = prepareDecoding(trials_100t{i},p);
%         [dec_100t{i},~] = decodingCore(dec_100t{i},nft_binned,prop,p);
%         dec_100t{i}.analysis.trainingSet = analyseDecoding(dec_100t{i},'trainingSet',1:dec_100t{i}.cv.numTrials_trainingSet,p);

        dec_100t{i} = prepareDecoding_100t_seq(trials_100t{i},p);
        [dec_100t{i},~] = decodingCore(dec_100t{i},nft_binned,prop,p);
        dec_100t{i}.analysis.trainingSet = analyseDecoding(dec_100t{i},'trainingSet',1:dec_100t{i}.cv.numTrials_trainingSet,p);
    end
end


%% Save

% if ops.dec.do_allTrials
%     save([path.filepart_outX,'dec_all.mat'],'dec_all','-v7.3');
%     disp(['--- Saved dec_all file as ',[path.filepart_outX,'dec_all.mat'],'.'])
% end
if ops.dec.do_100t
    save([path.filepart_outX,'dec_100t.mat'],'dec_100t','-v7.3');
    disp(['--- Saved dec_100t file as ',[path.filepart_outX,'dec_100t.mat'],'.'])
end


%% Return

if ops.close_figures
    close all;
end
end



