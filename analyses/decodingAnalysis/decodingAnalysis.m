function decodingAnalysis(info,iscell,ops,p,path,task,paq_beh,act,tng_all,tng_100t)
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
    dec_all = prepareDecoding(trials_all,p);
    [dec_all,~] = decodingCore(dec_all,nft_binned,prop,p);
    [dec_all,~] = decodingCore_seq(dec_all,nft_binned,prop,p,false,tng_all);
    dec_all.analysis.trainingSet = analyseDecoding(dec_all,'trainingSet',1:dec_all.cv.numTrials_trainingSet,p);
    dec_all.analysis.trainingSet_seq = analyseDecoding(dec_all,'trainingSet',1:dec_all.cv.numTrials_trainingSet,p,true);
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
        dec_100t{i} = prepareDecoding(trials_100t{i},p);
        [dec_100t{i},~] = decodingCore(dec_100t{i},nft_binned,prop,p);
        [dec_100t{i},~] = decodingCore_seq(dec_100t{i},nft_binned,prop,p,false,tng_100t{i});
        dec_100t{i}.analysis.trainingSet = analyseDecoding(dec_100t{i},'trainingSet',1:dec_100t{i}.cv.numTrials_trainingSet,p);
        dec_100t{i}.analysis.trainingSet_seq = analyseDecoding(dec_100t{i},'trainingSet',1:dec_100t{i}.cv.numTrials_trainingSet,p,true);
    end
end


%% Save

if ops.dec.do_allTrials
    save([path.filepart_outX,'dec_all.mat'],'dec_all','-v7.3');
    disp(['--- Saved dec_all file as ',[path.filepart_outX,'dec_all.mat'],'.'])
end
if ops.dec.do_100t
    save([path.filepart_outX,'dec_100t.mat'],'dec_100t','-v7.3');
    disp(['--- Saved dec_100t file as ',[path.filepart_outX,'dec_100t.mat'],'.'])
end


%% Return

if ops.close_figures
    close all;
end
end
