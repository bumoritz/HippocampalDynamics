function saveTraces(info,iscell,ops,p,path,task,s2p_meta,paq_beh,act)

%% Preparations

disp('--- Processing')

% pre-process
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.tng,p,paq_beh,iscell);
bl_trial_frames_binned = [prop.trial_frames_binned(1,:)-[1:5]';prop.trial_frames_binned(1:10,:)];
nft_bl_binned = reshape(nf_binned(:,bl_trial_frames_binned),[],size(bl_trial_frames_binned,1),prop.numTrials);
nft_bl_avg_binned = nanmean(nft_bl_binned,3);

% create trials and traces structs
trials_all_stimVersion = createTrialsStruct(task,1:prop.numTrials_incl,true);
[~,avgTraces_all_stimVersion,~] = createTracesStructs(nft_binned,trials_all_stimVersion);


%% Create str struct

str.prop = prop;
str.nft_bl_avg_binned = nft_bl_avg_binned;
str.avgTraces_all_stimVersion.A_var0 = avgTraces_all_stimVersion.A_var0;
str.avgTraces_all_stimVersion.A_var1 = avgTraces_all_stimVersion.A_var1;
str.avgTraces_all_stimVersion.X_var0 = avgTraces_all_stimVersion.X_var0;
str.avgTraces_all_stimVersion.X_var1 = avgTraces_all_stimVersion.X_var1;


%% Save str struct

str = orderfields(str);
save([path.filepart_out,'str.mat'],'str','-v7.3');
disp(['--- Saved str file as ',[path.filepart_out,'str.mat'],'.'])
