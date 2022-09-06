%% Preparations

% start with running all pre-processing like in tuningAnalysis
% load nem_all_cmpr

%%

wdw = p.general.bins_2nd1s; %p.general.bins_2nd1s, p.general.bins_2nd2s, p.general.bins_2nd3s

% B - A vs X
idcs_B = setdiff(find(nem_all_cmpr.sigM_sigTG_pos{11}),find(nem_all_cmpr.sigM_sigTG{16}));
activity_Atrials = nanmean(normAvgTraces_all.trialType.AB(idcs_B,wdw),2);
activity_Xtrials = nanmean(normAvgTraces_all.trialType.XB(idcs_B,wdw),2);
differenceOverSum_B = (activity_Atrials - activity_Xtrials) ./ (activity_Atrials + activity_Xtrials);

% Y - A vs X
idcs_Y = setdiff(find(nem_all_cmpr.sigM_sigTG_pos{12}),find(nem_all_cmpr.sigM_sigTG{16}));
activity_Xtrials = nanmean(normAvgTraces_all.trialType.XY(idcs_Y,wdw),2);
activity_Atrials = nanmean(normAvgTraces_all.trialType.AY(idcs_Y,wdw),2);
differenceOverSum_Y = (activity_Xtrials - activity_Atrials) ./ (activity_Xtrials + activity_Atrials);

figure;
subplot(1,2,1)
histogram(differenceOverSum_B);
xlabel(['Selectivity index',newline,'AB trials vs XB trials'])
ylabel(['Count'])
title(['Odour B cells (that are not L+R modulated)'])
subplot(1,2,2)
histogram(differenceOverSum_Y);
xlabel(['Selectivity index',newline,'XY trials vs AY trials'])
ylabel(['Count'])
title(['Odour Y cells (that are not L+R modulated)'])


%%

idx = idcs_B(76)
F = singleCellFigure(info,p,tng,idx,paq_beh);


%% Trialwise decoding error -> integration

prctile_dec = 20;
prctile_sel = 10;

wdw_decoding = 11:20; % used to be 1:26

% --- decoding errors ---

decodingErrors = nanmean(dec_all.analysis.completeSet.timeDecodingError_withinCat_s(:,21:26),2);

% AB trials
these_trials = trials_all.stimuli.AB;
these_decodingErrors = decodingErrors(these_trials);
decoding_AB_best = these_trials(find(these_decodingErrors <= prctile(these_decodingErrors,prctile_dec)));
decoding_AB_worst = these_trials(find(these_decodingErrors >= prctile(these_decodingErrors,100-prctile_dec)));

% XY trials
these_trials = trials_all.stimuli.XY;
these_decodingErrors = decodingErrors(these_trials);
decoding_XY_best = these_trials(find(these_decodingErrors <= prctile(these_decodingErrors,prctile_dec)));
decoding_XY_worst = these_trials(find(these_decodingErrors >= prctile(these_decodingErrors,100-prctile_dec)));

% AY trials
these_trials = trials_all.stimuli.AY;
these_decodingErrors = decodingErrors(these_trials);
decoding_AY_best = these_trials(find(these_decodingErrors <= prctile(these_decodingErrors,prctile_dec)));
decoding_AY_worst = these_trials(find(these_decodingErrors >= prctile(these_decodingErrors,100-prctile_dec)));

% XB trials
these_trials = trials_all.stimuli.XB;
these_decodingErrors = decodingErrors(these_trials);
decoding_XB_best = these_trials(find(these_decodingErrors <= prctile(these_decodingErrors,prctile_dec)));
decoding_XB_worst = these_trials(find(these_decodingErrors >= prctile(these_decodingErrors,100-prctile_dec)));


% --- selectivity indices ---

selectivity_B_AoverX = idcs_B(find(differenceOverSum_B >= prctile(differenceOverSum_B,100-prctile_sel)));
selectivity_B_XoverA = idcs_B(find(differenceOverSum_B <= prctile(differenceOverSum_B,prctile_sel)));
selectivity_Y_XoverA = idcs_Y(find(differenceOverSum_Y >= prctile(differenceOverSum_Y,100-prctile_sel)));
selectivity_Y_AoverX = idcs_Y(find(differenceOverSum_Y <= prctile(differenceOverSum_Y,prctile_sel)));


%% Figure: Trialwise decoding error -> integration

nrows=2; ncols=4;
F = default_figure([20,0.5,20,9.9]);

r=1; c=1; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_best = squeeze(nanmean(nanmean(nft_binned(selectivity_B_AoverX,wdw,decoding_AB_best),2),1)) - nanmean(nanmean(avgTraces_all.B(selectivity_B_AoverX,wdw),2),1);
this_data_worst = squeeze(nanmean(nanmean(nft_binned(selectivity_B_AoverX,wdw,decoding_AB_worst),2),1)) - nanmean(nanmean(avgTraces_all.B(selectivity_B_AoverX,wdw),2),1);
h=bar(diag([nanmean(this_data_best),nanmean(this_data_worst)]),'stacked','BaseValue',0,'FaceColor',p.col.darkGray,'EdgeColor','none');
scatter([1*ones(length(this_data_best),1);2*ones(length(this_data_worst),1)],[this_data_best;this_data_worst],'k')
ylim([-1,2])
xticks(1:2)
xticklabels({'Good seq','Bad seq'})
ylabel('Mean activity of neurons (zscore, rel to avg of all B trials)')
title('AB trials, B cells with A>X')

r=1; c=2; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_best = squeeze(nanmean(nanmean(nft_binned(selectivity_B_XoverA,wdw,decoding_AB_best),2),1)) - nanmean(nanmean(avgTraces_all.B(selectivity_B_XoverA,wdw),2),1);
this_data_worst = squeeze(nanmean(nanmean(nft_binned(selectivity_B_XoverA,wdw,decoding_AB_worst),2),1)) - nanmean(nanmean(avgTraces_all.B(selectivity_B_XoverA,wdw),2),1);
h=bar(diag([nanmean(this_data_best),nanmean(this_data_worst)]),'stacked','BaseValue',0,'FaceColor',p.col.darkGray,'EdgeColor','none');
scatter([1*ones(length(this_data_best),1);2*ones(length(this_data_worst),1)],[this_data_best;this_data_worst],'k')
ylim([-1,2])
xticks(1:2)
xticklabels({'Good seq','Bad seq'})
ylabel('Mean activity of neurons (zscore, rel to avg of all B trials)')
title('AB trials, B cells with X>A')

r=1; c=3; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_best = squeeze(nanmean(nanmean(nft_binned(selectivity_B_AoverX,wdw,decoding_XB_best),2),1)) - nanmean(nanmean(avgTraces_all.B(selectivity_B_AoverX,wdw),2),1);
this_data_worst = squeeze(nanmean(nanmean(nft_binned(selectivity_B_AoverX,wdw,decoding_XB_worst),2),1)) - nanmean(nanmean(avgTraces_all.B(selectivity_B_AoverX,wdw),2),1);
h=bar(diag([nanmean(this_data_best),nanmean(this_data_worst)]),'stacked','BaseValue',0,'FaceColor',p.col.darkGray,'EdgeColor','none');
scatter([1*ones(length(this_data_best),1);2*ones(length(this_data_worst),1)],[this_data_best;this_data_worst],'k')
ylim([-1,2])
xticks(1:2)
xticklabels({'Good seq','Bad seq'})
ylabel('Mean activity of neurons (zscore, rel to avg of all B trials)')
title('XB trials, B cells with A>X')

r=1; c=4; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_best = squeeze(nanmean(nanmean(nft_binned(selectivity_B_XoverA,wdw,decoding_XB_best),2),1)) - nanmean(nanmean(avgTraces_all.B(selectivity_B_XoverA,wdw),2),1);
this_data_worst = squeeze(nanmean(nanmean(nft_binned(selectivity_B_XoverA,wdw,decoding_XB_worst),2),1)) - nanmean(nanmean(avgTraces_all.B(selectivity_B_XoverA,wdw),2),1);
h=bar(diag([nanmean(this_data_best),nanmean(this_data_worst)]),'stacked','BaseValue',0,'FaceColor',p.col.darkGray,'EdgeColor','none');
scatter([1*ones(length(this_data_best),1);2*ones(length(this_data_worst),1)],[this_data_best;this_data_worst],'k')
ylim([-1,2])
xticks(1:2)
xticklabels({'Good seq','Bad seq'})
ylabel('Mean activity of neurons (zscore, rel to avg of all B trials)')
title('XB trials, B cells with X>A')

r=2; c=1; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_best = squeeze(nanmean(nanmean(nft_binned(selectivity_Y_XoverA,wdw,decoding_AY_best),2),1)) - nanmean(nanmean(avgTraces_all.Y(selectivity_Y_XoverA,wdw),2),1);
this_data_worst = squeeze(nanmean(nanmean(nft_binned(selectivity_Y_XoverA,wdw,decoding_AY_worst),2),1)) - nanmean(nanmean(avgTraces_all.Y(selectivity_Y_XoverA,wdw),2),1);
h=bar(diag([nanmean(this_data_best),nanmean(this_data_worst)]),'stacked','BaseValue',0,'FaceColor',p.col.darkGray,'EdgeColor','none');
scatter([1*ones(length(this_data_best),1);2*ones(length(this_data_worst),1)],[this_data_best;this_data_worst],'k')
ylim([-1,2])
xticks(1:2)
xticklabels({'Good seq','Bad seq'})
ylabel('Mean activity of neurons (zscore, rel to avg of all Y trials)')
title('AY trials, Y cells with X>A')

r=2; c=2; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_best = squeeze(nanmean(nanmean(nft_binned(selectivity_Y_AoverX,wdw,decoding_AY_best),2),1)) - nanmean(nanmean(avgTraces_all.Y(selectivity_Y_AoverX,wdw),2),1);
this_data_worst = squeeze(nanmean(nanmean(nft_binned(selectivity_Y_AoverX,wdw,decoding_AY_worst),2),1)) - nanmean(nanmean(avgTraces_all.Y(selectivity_Y_AoverX,wdw),2),1);
h=bar(diag([nanmean(this_data_best),nanmean(this_data_worst)]),'stacked','BaseValue',0,'FaceColor',p.col.darkGray,'EdgeColor','none');
scatter([1*ones(length(this_data_best),1);2*ones(length(this_data_worst),1)],[this_data_best;this_data_worst],'k')
ylim([-1,2])
xticks(1:2)
xticklabels({'Good seq','Bad seq'})
ylabel('Mean activity of neurons (zscore, rel to avg of all Y trials)')
title('AY trials, Y cells with A>X')

r=2; c=3; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_best = squeeze(nanmean(nanmean(nft_binned(selectivity_Y_XoverA,wdw,decoding_XY_best),2),1)) - nanmean(nanmean(avgTraces_all.Y(selectivity_Y_XoverA,wdw),2),1);
this_data_worst = squeeze(nanmean(nanmean(nft_binned(selectivity_Y_XoverA,wdw,decoding_XY_worst),2),1)) - nanmean(nanmean(avgTraces_all.Y(selectivity_Y_XoverA,wdw),2),1);
h=bar(diag([nanmean(this_data_best),nanmean(this_data_worst)]),'stacked','BaseValue',0,'FaceColor',p.col.darkGray,'EdgeColor','none');
scatter([1*ones(length(this_data_best),1);2*ones(length(this_data_worst),1)],[this_data_best;this_data_worst],'k')
ylim([-1,2])
xticks(1:2)
xticklabels({'Good seq','Bad seq'})
ylabel('Mean activity of neurons (zscore, rel to avg of all Y trials)')
title('XY trials, Y cells with X>A')

r=2; c=4; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_best = squeeze(nanmean(nanmean(nft_binned(selectivity_Y_AoverX,wdw,decoding_XY_best),2),1)) - nanmean(nanmean(avgTraces_all.Y(selectivity_Y_AoverX,wdw),2),1);
this_data_worst = squeeze(nanmean(nanmean(nft_binned(selectivity_Y_AoverX,wdw,decoding_XY_worst),2),1)) - nanmean(nanmean(avgTraces_all.Y(selectivity_Y_AoverX,wdw),2),1);
h=bar(diag([nanmean(this_data_best),nanmean(this_data_worst)]),'stacked','BaseValue',0,'FaceColor',p.col.darkGray,'EdgeColor','none');
scatter([1*ones(length(this_data_best),1);2*ones(length(this_data_worst),1)],[this_data_best;this_data_worst],'k')
ylim([-1,2])
xticks(1:2)
xticklabels({'Good seq','Bad seq'})
ylabel('Mean activity of neurons (zscore, rel to avg of all Y trials)')
title('XY trials, Y cells with A>X')

suptitle([num2str(prctile_dec),'% of best/worst decoding trials, ',num2str(prctile_sel),'% of most/least context-selective 2nd odour cells'])


%% Figure: seq2int_corr

nrows=2; ncols=4;
F = default_figure([20,0.5,20,9.9]);

r=1; c=1; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_x = decodingErrors(trials_all.stimuli.AB);
this_data_y = squeeze(nanmean(nanmean(nft_binned(selectivity_B_AoverX,wdw,trials_all.stimuli.AB),2),1)) - nanmean(nanmean(avgTraces_all.B(selectivity_B_AoverX,wdw),2),1);
scatter(this_data_x,this_data_y,'b')
[rho,pval] = fitLine(this_data_x,this_data_y,p.col.black);
ylim([-1,2]); xlim([0,3]);
xlabel('Decoding error (s)')
ylabel('Mean activity of neurons (rel to B trials avg)')
title(['AB trials, B cells (A>X), rho=',num2str(rho,1),', p=',num2str(pval,1)])

r=1; c=2; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_x = decodingErrors(trials_all.stimuli.AB);
this_data_y = squeeze(nanmean(nanmean(nft_binned(selectivity_B_XoverA,wdw,trials_all.stimuli.AB),2),1)) - nanmean(nanmean(avgTraces_all.B(selectivity_B_XoverA,wdw),2),1);
scatter(this_data_x,this_data_y,'b')
[rho,pval] = fitLine(this_data_x,this_data_y,p.col.black);
ylim([-1,2]); xlim([0,3]);
xlabel('Decoding error (s)')
ylabel('Mean activity of neurons (rel to B trials avg)')
title(['AB trials, B cells (X>A), rho=',num2str(rho,1),', p=',num2str(pval,1)])

r=1; c=3; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_x = decodingErrors(trials_all.stimuli.XB);
this_data_y = squeeze(nanmean(nanmean(nft_binned(selectivity_B_AoverX,wdw,trials_all.stimuli.XB),2),1)) - nanmean(nanmean(avgTraces_all.B(selectivity_B_AoverX,wdw),2),1);
scatter(this_data_x,this_data_y,'b')
[rho,pval] = fitLine(this_data_x,this_data_y,p.col.black);
ylim([-1,2]); xlim([0,3]);
xlabel('Decoding error (s)')
ylabel('Mean activity of neurons (rel to B trials avg)')
title(['XB trials, B cells (A>X), rho=',num2str(rho,1),', p=',num2str(pval,1)])

r=1; c=4; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_x = decodingErrors(trials_all.stimuli.XB);
this_data_y = squeeze(nanmean(nanmean(nft_binned(selectivity_B_XoverA,wdw,trials_all.stimuli.XB),2),1)) - nanmean(nanmean(avgTraces_all.B(selectivity_B_XoverA,wdw),2),1);
scatter(this_data_x,this_data_y,'b')
[rho,pval] = fitLine(this_data_x,this_data_y,p.col.black);
ylim([-1,2]); xlim([0,3]);
xlabel('Decoding error (s)')
ylabel('Mean activity of neurons (rel to B trials avg)')
title(['XB trials, B cells (X>A), rho=',num2str(rho,1),', p=',num2str(pval,1)])

r=2; c=1; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_x = decodingErrors(trials_all.stimuli.AY);
this_data_y = squeeze(nanmean(nanmean(nft_binned(selectivity_Y_XoverA,wdw,trials_all.stimuli.AB),2),1)) - nanmean(nanmean(avgTraces_all.Y(selectivity_Y_XoverA,wdw),2),1);
scatter(this_data_x,this_data_y,'b')
[rho,pval] = fitLine(this_data_x,this_data_y,p.col.black);
ylim([-1,2]); xlim([0,3]);
xlabel('Decoding error (s)')
ylabel('Mean activity of neurons (rel to B trials avg)')
title(['AY trials, B cells (X>A), rho=',num2str(rho,1),', p=',num2str(pval,1)])

r=2; c=2; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_x = decodingErrors(trials_all.stimuli.AY);
this_data_y = squeeze(nanmean(nanmean(nft_binned(selectivity_Y_AoverX,wdw,trials_all.stimuli.AB),2),1)) - nanmean(nanmean(avgTraces_all.Y(selectivity_Y_AoverX,wdw),2),1);
scatter(this_data_x,this_data_y,'b')
[rho,pval] = fitLine(this_data_x,this_data_y,p.col.black);
ylim([-1,2]); xlim([0,3]);
xlabel('Decoding error (s)')
ylabel('Mean activity of neurons (rel to B trials avg)')
title(['AY trials, B cells (A>X), rho=',num2str(rho,1),', p=',num2str(pval,1)])

r=2; c=3; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_x = decodingErrors(trials_all.stimuli.XY);
this_data_y = squeeze(nanmean(nanmean(nft_binned(selectivity_Y_XoverA,wdw,trials_all.stimuli.XB),2),1)) - nanmean(nanmean(avgTraces_all.Y(selectivity_Y_XoverA,wdw),2),1);
scatter(this_data_x,this_data_y,'b')
[rho,pval] = fitLine(this_data_x,this_data_y,p.col.black);
ylim([-1,2]); xlim([0,3]);
xlabel('Decoding error (s)')
ylabel('Mean activity of neurons (rel to B trials avg)')
title(['XY trials, B cells (X>A), rho=',num2str(rho,1),', p=',num2str(pval,1)])

r=2; c=4; subplot(nrows,ncols,(r-1)*ncols+c); hold on;
this_data_x = decodingErrors(trials_all.stimuli.XY);
this_data_y = squeeze(nanmean(nanmean(nft_binned(selectivity_Y_AoverX,wdw,trials_all.stimuli.XY),2),1)) - nanmean(nanmean(avgTraces_all.Y(selectivity_Y_AoverX,wdw),2),1);
scatter(this_data_x,this_data_y,'b')
[rho,pval] = fitLine(this_data_x,this_data_y,p.col.black);
ylim([-1,2]); xlim([0,3]);
xlabel('Decoding error (s)')
ylabel('Mean activity of neurons (rel to B trials avg)')
title(['XY trials, B cells (A>X), rho=',num2str(rho,1),', p=',num2str(pval,1)])

suptitle([num2str(prctile_sel),'% of most/least context-selective 2nd odour cells'])





