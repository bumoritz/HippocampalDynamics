%function inhibitionGridFigure

load('F:\SniffinHippo\AnalysisX\Stanage\Stanage_20210924\Stanage_20210924_nem_all.mat');

% Get cell indices

this_data = nem_all.testGroupResidual{10}.coefs(:,1+(1:5)); %
idcs_A_pos = find(this_data(:,1)>0);
idcs_A_pos = idcs_A_pos(randperm(length(idcs_A_pos)));
idcs_A_zero = find(this_data(:,1)==0);
idcs_A_zero = idcs_A_zero(randperm(length(idcs_A_zero)));
idcs_A_neg = find(this_data(:,1)<0);
idcs_A_neg = idcs_A_neg(randperm(length(idcs_A_neg)));
idcs_A_other = setdiff(setdiff(setdiff(find(iscell==1),idcs_A_pos),idcs_A_zero),idcs_A_neg);
idcs_A_sorted = [idcs_A_pos;idcs_A_zero;idcs_A_neg;idcs_A_other];

this_data = nem_all.testGroupResidual{13}.coefs(:,1+(6:10)); %
idcs_X_pos = find(this_data(:,1)>0);
idcs_X_pos = idcs_X_pos(randperm(length(idcs_X_pos)));
idcs_X_zero = find(this_data(:,1)==0);
idcs_X_zero = idcs_X_zero(randperm(length(idcs_X_zero)));
idcs_X_neg = find(this_data(:,1)<0);
idcs_X_neg = idcs_X_neg(randperm(length(idcs_X_neg)));
idcs_X_other = setdiff(setdiff(setdiff(find(iscell==1),idcs_X_pos),idcs_X_zero),idcs_X_neg);
idcs_X_sorted = [idcs_X_pos;idcs_X_zero;idcs_X_neg;idcs_X_other];

this_data = nem_all.testGroupResidual{16}.coefs(:,1+(11:15)); %
idcs_B_pos = find(this_data(:,1)>0);
idcs_B_pos = idcs_B_pos(randperm(length(idcs_B_pos)));
idcs_B_zero = find(this_data(:,1)==0);
idcs_B_zero = idcs_B_zero(randperm(length(idcs_B_zero)));
idcs_B_neg = find(this_data(:,1)<0);
idcs_B_neg = idcs_B_neg(randperm(length(idcs_B_neg)));
idcs_B_other = setdiff(setdiff(setdiff(find(iscell==1),idcs_B_pos),idcs_B_zero),idcs_B_neg);
idcs_B_sorted = [idcs_B_pos;idcs_B_zero;idcs_B_neg;idcs_B_other];

this_data = nem_all.testGroupResidual{19}.coefs(:,1+(16:20)); %
idcs_Y_pos = find(this_data(:,1)>0);
idcs_Y_pos = idcs_Y_pos(randperm(length(idcs_Y_pos)));
idcs_Y_zero = find(this_data(:,1)==0);
idcs_Y_zero = idcs_Y_zero(randperm(length(idcs_Y_zero)));
idcs_Y_neg = find(this_data(:,1)<0);
idcs_Y_neg = idcs_Y_neg(randperm(length(idcs_Y_neg)));
idcs_Y_other = setdiff(setdiff(setdiff(find(iscell==1),idcs_Y_pos),idcs_Y_zero),idcs_Y_neg);
idcs_Y_sorted = [idcs_Y_pos;idcs_Y_zero;idcs_Y_neg;idcs_Y_other];


%% Preparations

baselineWindow_1st = 1:10; %1:10;
responseWindow_1st = 16:20;
plottingWindow_1st = 1:25;

baselineWindow_2nd = 27:36; %27:36;
responseWindow_2nd = 42:46;
plottingWindow_2nd = 27:51;



act = dFF_beh; %F_hs_bs_beh; 


p = get_p; sync_beh = paq_beh; 
[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.inh,p,sync_beh,iscell);
trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
[traces_all,avgTraces_all,normAvgTraces_all] = createTracesStructs(nft_binned,trials_all);


%% Get normalised traces

% numTrials = 1:50; %1:500;

these_avgTraces = avgTraces_all.A(:,:);
this_baseline_mean = nanmean(these_avgTraces(:,baselineWindow_1st),2);
this_baseline_std = nanstd(these_avgTraces(:,baselineWindow_1st),[],2);
normAvgTraces_A = (these_avgTraces - this_baseline_mean) ./ this_baseline_std;
% normAvgTraces_A = these_avgTraces;

these_avgTraces = avgTraces_all.X(:,:);
this_baseline_mean = nanmean(these_avgTraces(:,baselineWindow_1st),2);
this_baseline_std = nanstd(these_avgTraces(:,baselineWindow_1st),[],2);
normAvgTraces_X = (these_avgTraces - this_baseline_mean) ./ this_baseline_std;
% normAvgTraces_X = these_avgTraces;

these_avgTraces = avgTraces_all.B(:,:);
this_baseline_mean = nanmean(these_avgTraces(:,baselineWindow_2nd),2);
this_baseline_std = nanstd(these_avgTraces(:,baselineWindow_2nd),[],2);
normAvgTraces_B = (these_avgTraces - this_baseline_mean) ./ this_baseline_std;
% normAvgTraces_B = these_avgTraces;

these_avgTraces = avgTraces_all.Y(:,:);
this_baseline_mean = nanmean(these_avgTraces(:,baselineWindow_2nd),2);
this_baseline_std = nanstd(these_avgTraces(:,baselineWindow_2nd),[],2);
normAvgTraces_Y = (these_avgTraces - this_baseline_mean) ./ this_baseline_std;
% normAvgTraces_Y = these_avgTraces;


%% Make figure - A cells

plt.clim = [-10,10]; plt.colormap = redblue;
%plt.clim = [0.1,0.3]; plt.colormap = parula;

nrows = 1; ncols = 4;
F = default_figure([20,0.5,20,9.9]);

r=1; c=1; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_A(idcs_A_sorted,plottingWindow_1st),NaN,[],p,info,plt);
ylabel('Cell (sorted by A-encoding)'); title('A trials');

r=1; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_X(idcs_A_sorted,plottingWindow_1st),NaN,[],p,info,plt);
ylabel('Cell (sorted by A-encoding)'); title('X trials');

r=1; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_B(idcs_A_sorted,plottingWindow_2nd),NaN,[],p,info,plt);
ylabel('Cell (sorted by A-encoding)'); title('B trials');

r=1; c=4; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_Y(idcs_A_sorted,plottingWindow_2nd),NaN,[],p,info,plt);
ylabel('Cell (sorted by A-encoding)'); title('Y trials');








%% Make figure - X cells

nrows = 1; ncols = 4;
F = default_figure([20,0.5,20,9.9]);

r=1; c=1; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_A(idcs_X_sorted,plottingWindow_1st),NaN,[],p,info,plt);
ylabel('Cell (sorted by X-encoding)'); title('A trials');

r=1; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_X(idcs_X_sorted,plottingWindow_1st),NaN,[],p,info,plt);
ylabel('Cell (sorted by X-encoding)'); title('X trials');

r=1; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_B(idcs_X_sorted,plottingWindow_2nd),NaN,[],p,info,plt);
ylabel('Cell (sorted by X-encoding)'); title('B trials');

r=1; c=4; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_Y(idcs_X_sorted,plottingWindow_2nd),NaN,[],p,info,plt);
ylabel('Cell (sorted by X-encoding)'); title('Y trials');


%% Make figure - B cells

nrows = 1; ncols = 4;
F = default_figure([20,0.5,20,9.9]);

r=1; c=1; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_A(idcs_B_sorted,plottingWindow_1st),NaN,[],p,info,plt);
ylabel('Cell (sorted by B-encoding)'); title('A trials');

r=1; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_X(idcs_B_sorted,plottingWindow_1st),NaN,[],p,info,plt);
ylabel('Cell (sorted by B-encoding)'); title('X trials');

r=1; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_B(idcs_B_sorted,plottingWindow_2nd),NaN,[],p,info,plt);
ylabel('Cell (sorted by B-encoding)'); title('B trials');

r=1; c=4; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_Y(idcs_B_sorted,plottingWindow_2nd),NaN,[],p,info,plt);
ylabel('Cell (sorted by B-encoding)'); title('Y trials');


%% Make figure - Y cells

nrows = 1; ncols = 4;
F = default_figure([20,0.5,20,9.9]);

r=1; c=1; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_A(idcs_Y_sorted,plottingWindow_1st),NaN,[],p,info,plt);
ylabel('Cell (sorted by Y-encoding)'); title('A trials');

r=1; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_X(idcs_Y_sorted,plottingWindow_1st),NaN,[],p,info,plt);
ylabel('Cell (sorted by Y-encoding)'); title('X trials');

r=1; c=3; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_B(idcs_Y_sorted,plottingWindow_2nd),NaN,[],p,info,plt);
ylabel('Cell (sorted by Y-encoding)'); title('B trials');

r=1; c=4; subplot(nrows,ncols,(r-1)*ncols+c);
heatMap_task(normAvgTraces_Y(idcs_Y_sorted,plottingWindow_2nd),NaN,[],p,info,plt);
ylabel('Cell (sorted by Y-encoding)'); title('Y trials');


%%

nrows = 2; ncols = 2;
F = default_figure;
plt={}; plt.c_abs=15; plt.c_label='Activity (z-score)'; 

r=1; c=1; subplot(nrows,ncols,(r-1)*ncols+c);
this_property = nanmean(normAvgTraces_A(:,responseWindow_1st),2);
plt.titleText='Response to odour A';
propertyMap(s2p_meta,iscell,this_property,1,plt);

r=1; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
this_property = nanmean(normAvgTraces_X(:,responseWindow_1st),2);
plt.titleText='Response to odour X';
propertyMap(s2p_meta,iscell,this_property,1,plt);

r=2; c=1; subplot(nrows,ncols,(r-1)*ncols+c);
this_property = nanmean(normAvgTraces_B(:,responseWindow_2nd),2);
plt.titleText='Response to odour B';
propertyMap(s2p_meta,iscell,this_property,1,plt);

r=2; c=2; subplot(nrows,ncols,(r-1)*ncols+c);
this_property = nanmean(normAvgTraces_Y(:,responseWindow_2nd),2);
plt.titleText='Response to odour Y';
propertyMap(s2p_meta,iscell,this_property,1,plt);



% they all do kind of the same thing for all odours, just stronger for 2nd odours
% link 2p stim inhibition to imaging inhibition story
% pool across all 4 odours, then 





