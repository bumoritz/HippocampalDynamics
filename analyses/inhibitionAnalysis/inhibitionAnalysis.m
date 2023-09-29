%function [inh] = inhibitionAnalysis(info,iscell,ops,p,path,task,s2p_meta,sync_beh,act)
%act = spks_beh; p = get_p; sync_beh = paq_beh; 
act = F_beh; p = get_p; sync_beh = paq_beh; 

%% Preparations

disp('--- Preparations')

% pre-process activity measure
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);
[events_binned] = binPaqEvents(paq_beh,task,p,prop);

% create trials, traces, avgTraces and normAvgTraces structs
trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
[traces_all,avgTraces_all,normAvgTraces_all] = createTracesStructs(nft_binned,trials_all);


%% Identifying inhibited neurons

% a) hard cut-off at zero
% b) statistical test: 1 s pre vs 1 s post
% c) GLM: early odour predictors


%% Correlation with running speed

nrows = 2; ncols = 2;
F = default_figure([20,0.5,20,9.9]);

[rho,pval] = corr(nf_binned(idcs_A_pos,:)',events_binned.velocity,'Type','Pearson','Rows','Complete');

figure;
histogram(rho)

[rho,pval] = corr(nf_binned(idcs_A_neg,:)',events_binned.velocity,'Type','Pearson','Rows','Complete');

figure;
histogram(rho)


%%

figure;

[rho,pval] = corr(nf_binned(idcs_A_pos,:)',nf_binned(idcs_A_pos,:)','Type','Pearson','Rows','Complete');
subplot(2,2,1)
histogram(rho)
title('A-pos')
subplot(2,2,2)
imshow(rho,[0,1])
title('A-pos')

[rho,pval] = corr(nf_binned(idcs_A_neg,:)',nf_binned(idcs_A_neg,:)','Type','Pearson','Rows','Complete');
subplot(2,2,3)
histogram(rho)
title('A-neg')
subplot(2,2,4)
imshow(rho,[0,1])
title('A-neg')

%%

[rho,pval] = corr(nf_binned(find(iscell==1),:)',nf_binned(find(iscell==1),:)','Type','Pearson','Rows','Complete');


cm = struct('GroupNumber',{27,31,25,21,24,33,29,18,32},...
    'Annotation',{'woody','almond','waxy','glue','vomit','fruity','fruity','Christmas','piss'},...
    'Color',{rgb('Brown'),rgb('Orange'),rgb('Yellow'),rgb('Purple'),rgb('Blue'),rgb('Red'),rgb('Red'),rgb('DarkOrange'),rgb('DarkBlue')});

cgo = clustergram(rho,...
    'Colormap','pink','Symmetric',false);




