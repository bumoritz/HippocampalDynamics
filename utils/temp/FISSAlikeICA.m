%% Get data

idx = 1;

Fraw = F_beh(idx,:);
Fhalo = Fhalo_beh(idx,:);
Fhalo_1 = Fhalo_1_beh(idx,:);
Fhalo_2 = Fhalo_2_beh(idx,:);
Fhalo_3 = Fhalo_3_beh(idx,:);
Fhalo_4 = Fhalo_4_beh(idx,:);


%% Make trial-averaged traces

p = get_p; sync_beh = paq_beh;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[~,~,Fraw_avg] = preprocessActivityMeasure(repmat(Fraw,length(iscell),1),p.inh,p,sync_beh,iscell);
Fraw_avg = nanmean(nanmean(Fraw_avg,3),1);
[~,~,Fhalo_avg] = preprocessActivityMeasure(repmat(Fhalo,length(iscell),1),p.inh,p,sync_beh,iscell);
Fhalo_avg = nanmean(nanmean(Fhalo_avg,3),1);
[~,~,Fhalo_1_avg] = preprocessActivityMeasure(repmat(Fhalo_1,length(iscell),1),p.inh,p,sync_beh,iscell);
Fhalo_1_avg = nanmean(nanmean(Fhalo_1_avg,3),1);
[~,~,Fhalo_2_avg] = preprocessActivityMeasure(repmat(Fhalo_2,length(iscell),1),p.inh,p,sync_beh,iscell);
Fhalo_2_avg = nanmean(nanmean(Fhalo_2_avg,3),1);
[~,~,Fhalo_3_avg] = preprocessActivityMeasure(repmat(Fhalo_3,length(iscell),1),p.inh,p,sync_beh,iscell);
Fhalo_3_avg = nanmean(nanmean(Fhalo_3_avg,3),1);
[~,~,Fhalo_4_avg] = preprocessActivityMeasure(repmat(Fhalo_4,length(iscell),1),p.inh,p,sync_beh,iscell);
Fhalo_4_avg = nanmean(nanmean(Fhalo_4_avg,3),1);


%% Plot raw traces

wdw = 1:10000;
left = 1;
right = 10000;

nrows = 6; ncols = 1; n=0;
F=default_figure([20,0.5,20,9.9]);

n=n+1; subplot(nrows,ncols,n); hold on;
plot(Fraw(:,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['F'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(Fhalo(:,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['Fhalo'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(Fhalo_1(:,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['Fhalo-1'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(Fhalo_2(:,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['Fhalo-2'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(Fhalo_3(:,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['Fhalo-3'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(Fhalo_4(:,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['Fhalo-4'])

suptitle(['idx=',num2str(idx)])


%% Plot averaged raw traces

wdw = 1:10000;

nrows = 6; ncols = 1; n=0;
F=default_figure([20,0.5,20,9.9]);

n=n+1; subplot(nrows,ncols,n); hold on;
plot(Fraw_avg)
title(['F'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(Fhalo_avg)
title(['Fhalo'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(Fhalo_1_avg)
title(['Fhalo-1'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(Fhalo_2_avg)
title(['Fhalo-2'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(Fhalo_3_avg)
title(['Fhalo-3'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(Fhalo_4_avg)
title(['Fhalo-4'])

suptitle(['idx=',num2str(idx)])


%% Do ICA

ica_in = [Fraw;Fhalo_1;Fhalo_2;Fhalo_3;Fhalo_4]';
mdl = rica(ica_in,size(ica_in,2));
ica_out = (transform(mdl,ica_in))';


%% Make trial-averaged ICA traces

p = get_p; sync_beh = paq_beh;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[~,~,ica_out_1_avg] = preprocessActivityMeasure(repmat(ica_out(1,:),length(iscell),1),p.inh,p,sync_beh,iscell);
ica_out_1_avg = nanmean(nanmean(ica_out_1_avg,3),1);
[~,~,ica_out_2_avg] = preprocessActivityMeasure(repmat(ica_out(2,:),length(iscell),1),p.inh,p,sync_beh,iscell);
ica_out_2_avg = nanmean(nanmean(ica_out_2_avg,3),1);
[~,~,ica_out_3_avg] = preprocessActivityMeasure(repmat(ica_out(3,:),length(iscell),1),p.inh,p,sync_beh,iscell);
ica_out_3_avg = nanmean(nanmean(ica_out_3_avg,3),1);
[~,~,ica_out_4_avg] = preprocessActivityMeasure(repmat(ica_out(4,:),length(iscell),1),p.inh,p,sync_beh,iscell);
ica_out_4_avg = nanmean(nanmean(ica_out_4_avg,3),1);
[~,~,ica_out_5_avg] = preprocessActivityMeasure(repmat(ica_out(5,:),length(iscell),1),p.inh,p,sync_beh,iscell);
ica_out_5_avg = nanmean(nanmean(ica_out_5_avg,3),1);


%% Plot ICA_out traces

wdw = 1:10000;
left = 1;
right = 10000;

nrows = 5; ncols = 1; n=0;
F=default_figure([20,0.5,20,9.9]);

n=n+1; subplot(nrows,ncols,n); hold on;
plot(ica_out(1,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['IC1'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(ica_out(2,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['IC2'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(ica_out(3,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['IC3'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(ica_out(4,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['IC4'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(ica_out(5,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['IC5'])

suptitle(['idx=',num2str(idx)])


%% Plot averaged raw traces

wdw = 1:10000;

nrows = 5; ncols = 1; n=0;
F=default_figure([20,0.5,20,9.9]);

n=n+1; subplot(nrows,ncols,n); hold on;
plot(ica_out_1_avg)
title(['IC1'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(ica_out_2_avg)
title(['IC2'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(ica_out_3_avg)
title(['IC3'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(ica_out_4_avg)
title(['IC4'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(ica_out_5_avg)
title(['IC5'])

suptitle(['idx=',num2str(idx)])


%% Do NMF

nmf_in = [Fraw;Fhalo_1;Fhalo_2;Fhalo_3;Fhalo_4]';
[W,H] = nnmf(nmf_in,size(nmf_in,2));
nmf_out = W';
%ica_out = (transform(mdl,nmf_in))';


%% Plot nmf_out traces

wdw = 1:10000;
left = 1;
right = 10000;

nrows = 5; ncols = 1; n=0;
F=default_figure([20,0.5,20,9.9]);

n=n+1; subplot(nrows,ncols,n); hold on;
plot(nmf_out(1,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['NMF1'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(nmf_out(2,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['NMF2'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(nmf_out(3,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['NMF3'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(nmf_out(4,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['NMF4'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(nmf_out(5,wdw))
for i=1:10
    xline(paq_beh.sync(i));
end
xlim([left,right])
title(['NMF5'])

suptitle(['idx=',num2str(idx)])




%% --- BACK TO NEUROPIL SUBRACTION ---

figure;
scatter(Fhalo,Fraw)
xlabel('Neuropil signal')
ylabel('Raw signal')


%%

mean_activity = nanmean(F_beh(find(iscell)==1,:),1);
mean_neuropil = nanmean(Fneu_beh(find(iscell)==1,:),1);
mean_halo = nanmean(Fhalo_beh(find(iscell)==1,:),1);

p = get_p; sync_beh = paq_beh;
p.inh.smoothingSd_preBinning = 0;
p.inh.zscore = 0;
[~,~,mean_activity_avg] = preprocessActivityMeasure(repmat(mean_activity,length(iscell),1),p.inh,p,sync_beh,iscell);
mean_activity_avg = nanmean(nanmean(mean_activity_avg,3),1);
[~,~,mean_neuropil_avg] = preprocessActivityMeasure(repmat(mean_neuropil,length(iscell),1),p.inh,p,sync_beh,iscell);
mean_neuropil_avg = nanmean(nanmean(mean_neuropil_avg,3),1);
[~,~,mean_halo_avg] = preprocessActivityMeasure(repmat(mean_halo,length(iscell),1),p.inh,p,sync_beh,iscell);
mean_halo_avg = nanmean(nanmean(mean_halo_avg,3),1);


%% Plot averaged raw traces

wdw = 1:10000;

nrows = 3; ncols = 1; n=0;
F=default_figure([20,0.5,20,9.9]);

n=n+1; subplot(nrows,ncols,n); hold on;
plot(mean_activity_avg)
title(['mean activity'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(mean_neuropil_avg)
title(['mean neuropil'])

n=n+1; subplot(nrows,ncols,n); hold on;
plot(mean_halo_avg)
title(['mean halo'])


%% Approaches left

% - brutal common-average referencing
% - subtracting highly multiplied halo signal
% - 












