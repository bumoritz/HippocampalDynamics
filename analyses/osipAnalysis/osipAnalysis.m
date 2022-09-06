%function osipAnalysis(info,iscell,ops,p,path,task,paq_beh,act,tng_nem,nem_all)
act = spks_beh;

%% Preparations

disp('--- Preparations.')

% pre-process activity measure and paq events
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.osip,p,paq_beh,iscell);
[events_binned] = binPaqEvents(paq_beh,task,p,prop);

% create trials, traces and avgTraces structs
trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
[traces_all,avgTraces_all,~] = createTracesStructs(nft_binned,trials_all);



disp('--- Preparations')

% pre-process activity measure
[prop,~,nft_binned] = preprocessActivityMeasure(act,p.sqn,p,sync_beh,iscell);

% create trials, traces, avgTraces and normAvgTraces structs
if ops.sqn.do_allTrials
    trials_all = createTrialsStruct(task,1:prop.numTrials);
    [traces_all,avgTraces_all,normAvgTraces_all] = createTracesStructs(nft_binned,trials_all);
end
if ops.sqn.do_100t
    trials_100t = {}; traces_100t = {}; avgTraces_100t = {};
    for i=1:floor(prop.numTrials/100)
        trials_100t{i} = createTrialsStruct(task,(i-1)*100+1:i*100);
        [traces_100t{i},avgTraces_100t{i},~] = createTracesStructs(nft_binned,trials_100t{i});
    end
end




%%

% odour B cells (testGroup==11)

find(nem_all_cmpr.sigM_sigTG{11}==1)








%% Return

% if ops.close_figures
%     close all;
% end
% end

