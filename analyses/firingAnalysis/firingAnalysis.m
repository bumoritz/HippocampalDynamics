%%function firingAnalysis
act = spks_beh; p = get_p; sync_beh = paq_beh;


%% Preparations

disp('--- Preparations')

% pre-process activity measure
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.tng,p,sync_beh,iscell);


%%


this_data = nanmean(nf_binned(find(iscell),:),2);



%%

figure;
histogram(this_data)