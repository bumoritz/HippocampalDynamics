function encodingCellActivityAnalysis(info,iscell,ops,p,path,task,paq_beh,act,nem_all_cmpr,nem_100t_cmpr)
% act = spks_beh; p = get_p;

%% Preparations

disp('--- Preparations')

% pre-process activity measure
[prop,~,nft_binned] = preprocessActivityMeasure(act,p.eca,p,paq_beh,iscell);

% create trials, traces and avgTraces structs
if ops.eca.do_allTrials
    trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
    trials_all = balanceTrialSelelection(trials_all,info,task);
    [traces_all,avgTraces_all,~] = createTracesStructs(nft_binned,trials_all);
end

% repeat for unsmoothed data
p.eca_unsmoothed = p.eca;
p.eca_unsmoothed.smoothingSd_preBinning = 0;
[prop,~,nft_binned_unsmoothed] = preprocessActivityMeasure(act,p.eca_unsmoothed,p,paq_beh,iscell);
if ops.eca.do_allTrials
    [traces_unsmoothed_all,avgTraces_unsmoothed_all,~] = createTracesStructs(nft_binned_unsmoothed,trials_all);
end


%% Go for it!

if ops.sqn.do_allTrials
    
    eca_all.dyn = getDynamicECA(prop,traces_all,avgTraces_all,nem_all_cmpr,p);
    eca_all.dyn_unsmoothed = getDynamicECA(prop,traces_unsmoothed_all,avgTraces_unsmoothed_all,nem_all_cmpr,p);
end


%% Save and return

if ops.eca.do_allTrials
    eca_all = orderfields(eca_all);
    save([path.filepart_out,'eca_all.mat'],'eca_all','-v7.3');
    disp(['--- Saved eca_all file as ',[path.filepart_out,'eca_all.mat'],'.'])
end
if ops.eca.do_allTrials
    eca = eca_all;
end

if ops.close_figures
    close all;
end
%end






