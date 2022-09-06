%% Preparations

act = spks_beh;

% pre-process activity measure
[prop,~,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);
trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
[traces_all,avgTraces_all,~] = createTracesStructs(nft_binned,trials_all);




%%

q = 100;

[coeff,Data_PCA,latent,tsquared,explained,mu] = pca(spks_beh','NumComponents',q);

disp([num2str(sum(explained(1:q)))])

mdl = rica(Data_PCA, q);
Data_ICA = transform(mdl,Data_PCA);


%%

act = repmat(Data_ICA',20,1);

% pre-process activity measure
[prop,~,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);
trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
[traces_all,avgTraces_all,~] = createTracesStructs(nft_binned,trials_all);


a = nanmean(nft_binned,3);

figure
for i=1:100
    hold on;
    plot(a(i,:))
end




%%

p.inh.smoothingSd_preBinning=0;
p.inh.zscore=0;

idx = 6;

data_in = [F_beh(6,:);Fneu_beh(6,:)]';
mdl = rica(data_in,2);
data_out = (transform(mdl,data_in))';

act = repmat(nanmean(data_out(1,:),1),2000,1);
[~,~,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);
data_out_main_avg = nanmean(nft_binned,3);

act = repmat(nanmean(data_out(2,:),1),2000,1);
[~,~,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);
data_out_avg = nanmean(nft_binned,3);



%% Simple ICA with F and Fneu

data_in = [this_F;this_Fneu]';

mdl = rica(data_in,2);
data_out = (transform(mdl,data_in))';

data_out(2,:)=0;
data_out_tformed = (data_out'.*mdl.TransformWeights(1,:))';

act = repmat(nanmean(data_out_tformed(1,:),1),2000,1);
[~,~,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);
data_out_main_tformed_avg = nanmean(nft_binned,3);




%% Fissa-like ICA with F and 4 Fhalo










%% Example raw traces

idx = 1;

wdw = 1:10000;
left = 1;
right = 10000;

F=default_figure([20,0.5,20,9.9]);

%subplot(4,1,1); hold on;
plot(F_beh(idx,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['F (idx=',num2str(idx),')'])

%subplot(4,1,2); hold on;
plot(Fneu_beh(idx,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['Fneu (idx=',num2str(idx),')'])

%subplot(4,1,3); hold on;
this_data = F_beh(idx,wdw) - nansum(data_out(2:5,:),1);
plot(this_data)
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['F ic sub (idx=',num2str(idx),')'])

% subplot(4,1,3); hold on;
% plot(Fhalo_beh(idx,wdw))
% for i=1:10
%     xline(paq_beh.sync(i));
% end
% xlim([left,right])
% title(['Fhalo (idx=',num2str(idx),')'])





