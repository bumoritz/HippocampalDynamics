%% ExampleImprintingSnake

info = get_info;

this_data_catch = [vertcat(avgTraces_recruited_seq.A_var0__Acatchonly_A{:,[2,4]});vertcat(avgTraces_recruited_seq.X_var0__Xcatchonly_X{:,[2,4]})];
this_data_stim = [vertcat(avgTraces_recruited_seq.A_var1__Acatchonly_A{:,[2,4]});vertcat(avgTraces_recruited_seq.X_var1__Xcatchonly_X{:,[2,4]})];
this_data_bsl = [vertcat(avgTraces_recruited_seq.nft_bl_avg_binned__Acatchonly_A{:,[2,4]});vertcat(avgTraces_recruited_seq.nft_bl_avg_binned__Xcatchonly_X{:,[2,4]})];

this_data_stimTime = [vertcat(stimTime_seq.Acatchonly_A{[1:22,24:end],[2,4]});vertcat(stimTime_seq.Xcatchonly_X{[1:22,24:end],[2,4]})];
this_data_peakTime = [vertcat(peakTime_seq.Acatchonly_A{[1:22,24:end],[2,4]});vertcat(peakTime_seq.Xcatchonly_X{[1:22,24:end],[2,4]})];
these_idcs = find(abs(this_data_stimTime-this_data_peakTime)<2);

% comment out if necessary
% this_data_catch = this_data_catch(these_idcs,:);
% this_data_stim = this_data_stim(these_idcs,:);
% this_data_bsl = this_data_bsl(these_idcs,:);


%% a) Recruited cells. Stim trials sorted by stim trials.

this_data = this_data_stim;

this_normalisationData = this_data(:,p.general.bins_analysisWindow);
this_normalisationData_bl = this_data(:,p.general.bins_baselineWindow);
this_lower = nanmean(this_data_bsl,2); %nanmin(this_normalisationData,[],2);
this_upper = nanmax(this_normalisationData,[],2);
this_data_norm = (this_data-this_lower) ./ (this_upper-this_lower);

[~,temp] = nanmax(this_data_norm(:,p.general.bins_analysisWindow),[],2);
[~,this_order] = sort(temp);
this_order = flip(this_order);

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; 
heatMap_task(this_data_norm(this_order,:),NaN,[],p,info,plt);
xlim([6,53])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off
title('Recruited cells. Stim trials sorted by stim trials.')


%% b) Recruited cells. Catch trials sorted by stim trials.

this_normalisationData = this_data_catch(:,p.general.bins_analysisWindow);
this_normalisationData_bl = this_data_catch(:,p.general.bins_baselineWindow);
this_lower = nanmean(this_data_bsl,2); %nanmin(this_normalisationData,[],2);
this_upper = nanmax(this_normalisationData,[],2);
this_data_norm = (this_data_catch-this_lower) ./ (this_upper-this_lower);

[~,temp] = nanmax(this_data_stim(:,p.general.bins_analysisWindow),[],2);
[~,this_order] = sort(temp);
this_order = flip(this_order);

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; 
heatMap_task(this_data_norm(this_order,:),NaN,[],p,info,plt);
xlim([6,53])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off
title('Recruited cells. Catch trials sorted by stim trials.')


%% c) Recruited cells. Catch trials sorted by catch trials.

this_data = this_data_catch;

this_normalisationData = this_data(:,p.general.bins_analysisWindow);
this_normalisationData_bl = this_data(:,p.general.bins_baselineWindow);
this_lower = nanmean(this_data_bsl,2); %nanmin(this_normalisationData,[],2);
this_upper = nanmax(this_normalisationData,[],2);
this_data_norm = (this_data-this_lower) ./ (this_upper-this_lower);

[~,temp] = nanmax(this_data_norm(:,p.general.bins_analysisWindow),[],2);
[~,this_order] = sort(temp);
this_order = flip(this_order);

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; 
heatMap_task(this_data_norm(this_order,:),NaN,[],p,info,plt);
xlim([6,53])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off
title('Recruited cells. Catch trials sorted by catch trials.')


%% d) Recruited cells. Stim trials sorted by catch trials.

this_normalisationData = this_data_stim(:,p.general.bins_analysisWindow);
this_normalisationData_bl = this_data_stim(:,p.general.bins_baselineWindow);
this_lower = nanmean(this_data_bsl,2); %nanmin(this_normalisationData,[],2);
this_upper = nanmax(this_normalisationData,[],2);
this_data_norm = (this_data_stim-this_lower) ./ (this_upper-this_lower);

[~,temp] = nanmax(this_data_catch(:,p.general.bins_analysisWindow),[],2);
[~,this_order] = sort(temp);
this_order = flip(this_order);

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; 
heatMap_task(this_data_norm(this_order,:),NaN,[],p,info,plt);
xlim([6,53])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off
title('Recruited cells. Stim trials sorted by catch trials.')


%% Fig2_ExampleImprintingSnake_catchTrials

this_normalisationData_stim = this_data_stim(:,p.general.bins_analysisWindow);
this_normalisationData_catch = this_data_catch(:,p.general.bins_analysisWindow);
this_normalisationData_pooled = [this_normalisationData_stim,this_normalisationData_catch];

this_data = this_data_catch;

this_normalisationData = this_data(:,p.general.bins_analysisWindow);
this_normalisationData_bl = this_data(:,p.general.bins_baselineWindow);
this_lower = nanmean(this_data_bsl,2); %nanmin(this_normalisationData,[],2);
this_upper = nanmax(this_normalisationData,[],2); %nanmax(this_normalisationData_pooled,[],2);
this_data_norm = (this_data-this_lower) ./ (this_upper-this_lower);

[~,temp] = nanmax(this_data_norm(:,p.general.bins_analysisWindow),[],2);
[these_peak_locations,this_order] = sort(temp);
this_order = flip(this_order);

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1;
temp = this_data_norm(this_order,:);
heatMap_task(temp(1:59,:),NaN,[],p,info,plt); hold on;
h=plot(flip(these_peak_locations(end-58:end))+p.general.bins_analysisWindow(1)-1,1:length(these_peak_locations(end-58:end)),'w.','MarkerSize',5);

xlim([6,53])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off
%title('Recruited cells. Catch trials sorted by catch trials.')

% save plot
savefig(F,[save_root_fig,'\Fig2_ExampleImprintingSnake_catchTrials.fig']);
saveas(F,[save_root_png,'\Fig2_ExampleImprintingSnake_catchTrials.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleImprintingSnake_catchTrials.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_ExampleImprintingSnake_stimTrials

this_normalisationData_stim = this_data_stim(:,p.general.bins_analysisWindow);
this_normalisationData_catch = this_data_catch(:,p.general.bins_analysisWindow);
this_normalisationData_pooled = [this_normalisationData_stim,this_normalisationData_catch];

this_normalisationData = this_data_stim(:,p.general.bins_analysisWindow);
this_normalisationData_bl = this_data_stim(:,p.general.bins_baselineWindow);
this_lower = nanmean(this_data_bsl,2); % nanmin(this_normalisationData,[],2);
this_upper = nanmax(this_normalisationData,[],2); %nanmax(this_normalisationData_pooled,[],2);
this_data_norm = (this_data_stim-this_lower) ./ (this_upper-this_lower);

[~,temp] = nanmax(this_data_catch(:,p.general.bins_analysisWindow),[],2);
[a,this_order] = sort(temp);
this_order = flip(this_order);

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; 
temp = this_data_norm(this_order,:);
heatMap_task(temp(1:59,:),NaN,[],p,info,plt); hold on;
h=plot(flip(these_peak_locations(end-58:end))+p.general.bins_analysisWindow(1)-1,1:length(these_peak_locations(end-58:end)),'w.','MarkerSize',5);

xlim([6,53])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off
%title('Recruited cells. Stim trials sorted by catch trials. Cropped.')

% save plot
savefig(F,[save_root_fig,'\Fig2_ExampleImprintingSnake_stimTrials.fig']);
saveas(F,[save_root_png,'\Fig2_ExampleImprintingSnake_stimTrials.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleImprintingSnake_stimTrials.pdf']); set(gcf,'Color',[1,1,1])



%% Fig2_ExampleImprintingSnake_CherryPickVersion

numCells_total = 233;
this_selection = [3,5,6,9,10,12,13,14,19,...
    27,30,33,35,39,40,44,45,46,50,53,66,72,92,...
    210,212,214,217,227,228,230];
numCells_selected = length(this_selection);

this_normalisationData_stim = this_data_stim(:,p.general.bins_analysisWindow);
this_normalisationData_catch = this_data_catch(:,p.general.bins_analysisWindow);
this_normalisationData_pooled = [this_normalisationData_stim,this_normalisationData_catch];

this_data = this_data_catch;

this_normalisationData = this_data(:,p.general.bins_analysisWindow);
this_normalisationData_bl = this_data(:,p.general.bins_baselineWindow);
this_lower = nanmean(this_data_bsl,2); %nanmin(this_normalisationData,[],2);
this_upper = nanmax(this_normalisationData,[],2); %nanmax(this_normalisationData_pooled,[],2);
this_data_norm = (this_data-this_lower) ./ (this_upper-this_lower);

[~,temp] = nanmax(this_data_norm(:,p.general.bins_analysisWindow),[],2);
[these_peak_locations,this_order] = sort(temp);
this_order = flip(this_order);

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1;
temp = this_data_norm(this_order,:);
heatMap_task(temp(this_selection,:),NaN,[],p,info,plt); hold on;
temp = flip(these_peak_locations);
h=plot(temp(this_selection)+p.general.bins_analysisWindow(1)-1,1:numCells_selected,'w.','MarkerSize',5);

xlim([6,53])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off
%title('Recruited cells. Catch trials sorted by catch trials.')

% save plot
% savefig(F,[save_root_fig,'\Fig2_ExampleImprintingSnake_CherryPickVersion_catchTrials.fig']);
% saveas(F,[save_root_png,'\Fig2_ExampleImprintingSnake_CherryPickVersion_catchTrials.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleImprintingSnake_CherryPickVersion_catchTrials.pdf']); set(gcf,'Color',[1,1,1])


this_normalisationData_stim = this_data_stim(:,p.general.bins_analysisWindow);
this_normalisationData_catch = this_data_catch(:,p.general.bins_analysisWindow);
this_normalisationData_pooled = [this_normalisationData_stim,this_normalisationData_catch];

this_normalisationData = this_data_stim(:,p.general.bins_analysisWindow);
this_normalisationData_bl = this_data_stim(:,p.general.bins_baselineWindow);
this_lower = nanmean(this_data_bsl,2); % nanmin(this_normalisationData,[],2);
this_upper = nanmax(this_normalisationData,[],2); %nanmax(this_normalisationData_pooled,[],2);
this_data_norm = (this_data_stim-this_lower) ./ (this_upper-this_lower);

[~,temp] = nanmax(this_data_catch(:,p.general.bins_analysisWindow),[],2);
[a,this_order] = sort(temp);
this_order = flip(this_order);

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; 
temp = this_data_norm(this_order,:);
heatMap_task(temp(this_selection,:),NaN,[],p,info,plt); hold on;
temp = flip(these_peak_locations);
h=plot(temp(this_selection)+p.general.bins_analysisWindow(1)-1,1:numCells_selected,'w.','MarkerSize',5);

xlim([6,53])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off
%title('Recruited cells. Stim trials sorted by catch trials. Cropped.')

% save plot
% savefig(F,[save_root_fig,'\Fig2_ExampleImprintingSnake_CherryPickVersion_stimTrials.fig']);
% saveas(F,[save_root_png,'\Fig2_ExampleImprintingSnake_CherryPickVersion_stimTrials.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleImprintingSnake_CherryPickVersion_stimTrials.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_ExampleImprintingSnake_CherryPickVersion - NEW

% 25 dots -> select 1 per bin
% highest idcs are at very beginning of sequence

% full
% this_selection = this_order_saved([6,11,21,32,38,...
%     [],[],56,60,65,...
%     [],85,97,102,108,...
%     114,120,136,154,171,...
%     176,202,218,230,232]);

% without first two
this_selection = this_order_saved([[],[],21,32,38,...
    [],[],56,60,65,...
    [],85,96,[],108,...
    114,120,136,154,171,...
    176,202,218,230,232]); % remove first two

% correct:
% this_selection = [];

numCells_total = 233;
numCells_selected = length(this_selection);

this_data_stim_sel = this_data_stim(this_selection,:);
this_data_catch_sel = this_data_catch(this_selection,:);
this_data_bsl_sel = this_data_bsl(this_selection,:);

%%% a) Recruited cells. Stim trials sorted by stim trials.

this_normalisationData = this_data_stim_sel(:,p.general.bins_analysisWindow);
this_lower = nanmean(this_data_bsl_sel,2); %nanmin(this_normalisationData,[],2);
this_upper = nanmax(this_normalisationData,[],2);
this_data_norm = (this_data_stim_sel-this_lower) ./ (this_upper-this_lower);

[~,temp] = nanmax(this_data_stim_sel(:,p.general.bins_analysisWindow),[],2);
[these_peak_locations,this_order] = sort(temp);
this_order = flip(this_order);

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; 
heatMap_task(this_data_norm(this_order,:),NaN,[],p,info,plt); hold on;
h=plot(flip(these_peak_locations)+p.general.bins_analysisWindow(1)-1,1:numCells_selected,'w.','MarkerSize',5);
xlim([6,53])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off
%title('Recruited cells. Stim trials sorted by stim trials.')

% save plot
savefig(F,[save_root_fig,'\Fig2_ExampleImprintingSnake_CherryPickVersion_stimTrials.fig']);
saveas(F,[save_root_png,'\Fig2_ExampleImprintingSnake_CherryPickVersion_stimTrials.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleImprintingSnake_CherryPickVersion_stimTrials.pdf']); set(gcf,'Color',[1,1,1])


%%% b) Recruited cells. Catch trials sorted by stim trials.

this_normalisationData = this_data_catch_sel(:,p.general.bins_analysisWindow);
this_lower = nanmean(this_data_bsl_sel,2); %nanmin(this_normalisationData,[],2);
this_upper = nanmax(this_normalisationData,[],2);
this_data_norm = (this_data_catch_sel-this_lower) ./ (this_upper-this_lower);

% [~,temp] = nanmax(this_data_stim_sel(:,p.general.bins_analysisWindow),[],2);
% [~,this_order] = sort(temp);
% this_order = flip(this_order);

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; 
heatMap_task(this_data_norm(this_order,:),NaN,[],p,info,plt); hold on;
h=plot(flip(these_peak_locations)+p.general.bins_analysisWindow(1)-1,1:numCells_selected,'w.','MarkerSize',5);
xlim([6,53])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off
%title('Recruited cells. Catch trials sorted by stim trials.')

% save plot
savefig(F,[save_root_fig,'\Fig2_ExampleImprintingSnake_CherryPickVersion_catchTrials.fig']);
saveas(F,[save_root_png,'\Fig2_ExampleImprintingSnake_CherryPickVersion_catchTrials.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleImprintingSnake_CherryPickVersion_catchTrials.pdf']); set(gcf,'Color',[1,1,1])



%% d3) Recruited cells. Stim trials sorted by catch trials. Cropped. Then sorted by stim trials.

% this_normalisationData_stim = this_data_stim(:,p.general.bins_analysisWindow);
% this_normalisationData_catch = this_data_catch(:,p.general.bins_analysisWindow);
% this_normalisationData_pooled = [this_normalisationData_stim,this_normalisationData_catch];
% 
% this_normalisationData = this_data_stim(:,p.general.bins_analysisWindow);
% this_normalisationData_bl = this_data_stim(:,p.general.bins_baselineWindow);
% this_lower = nanmean(this_data_bsl,2); % nanmin(this_normalisationData,[],2);
% this_upper = nanmax(this_normalisationData,[],2); %nanmax(this_normalisationData_pooled,[],2);
% this_data_norm = (this_data_stim-this_lower) ./ (this_upper-this_lower);
% 
% [~,temp] = nanmax(this_data_catch(:,p.general.bins_analysisWindow),[],2);
% [a,this_order] = sort(temp);
% this_order = flip(this_order);
% 
% F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;
% plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; 
% temp = this_data_norm(this_order,:);
% temp0 = temp(1:59,:);
% 
% [~,temp] = nanmax(temp0(:,p.general.bins_analysisWindow),[],2);
% [a,this_order] = sort(temp);
% this_order_X = flip(this_order);
% 
% heatMap_task(temp0(this_order_X,:),NaN,[],p,info,plt);
% xlim([6,53])
% set(gca,'xtick',[]);
% set(gca,'xticklabel',[]);
% set(gca,'ytick',[]);
% set(gca,'yticklabel',[]);
% set(gca,'box','off');
% set(gca,'xcolor','none');
% set(gca,'ycolor','none');
% box off
% title('Recruited cells, cropped. Stim trials sorted by stim trials.')


%% c3) Recruited cells. Catch trials sorted by catch trials. Cropped. Then sorted by stim trials.

% this_normalisationData_stim = this_data_stim(:,p.general.bins_analysisWindow);
% this_normalisationData_catch = this_data_catch(:,p.general.bins_analysisWindow);
% this_normalisationData_pooled = [this_normalisationData_stim,this_normalisationData_catch];
% 
% this_data = this_data_catch;
% 
% this_normalisationData = this_data(:,p.general.bins_analysisWindow);
% this_normalisationData_bl = this_data(:,p.general.bins_baselineWindow);
% this_lower = nanmean(this_data_bsl,2); %nanmin(this_normalisationData,[],2);
% this_upper = nanmax(this_normalisationData,[],2); %nanmax(this_normalisationData_pooled,[],2);
% this_data_norm = (this_data-this_lower) ./ (this_upper-this_lower);
% 
% [~,temp] = nanmax(this_data_norm(:,p.general.bins_analysisWindow),[],2);
% [~,this_order] = sort(temp);
% this_order = flip(this_order);
% 
% F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;
% plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1;
% temp = this_data_norm(this_order,:);
% temp0 = temp(1:59,:);
% 
% heatMap_task(temp0(this_order_X,:),NaN,[],p,info,plt);
% xlim([6,53])
% set(gca,'xtick',[]);
% set(gca,'xticklabel',[]);
% set(gca,'ytick',[]);
% set(gca,'yticklabel',[]);
% set(gca,'box','off');
% set(gca,'xcolor','none');
% set(gca,'ycolor','none');
% box off
% title('Recruited cells. Catch trials sorted by catch trials.')


%% Plot only those where a cell formed within 1 s of the stim location





