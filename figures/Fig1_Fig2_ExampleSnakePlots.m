%% Fig_ExampleSnakePlots

%% Example pooled snakes

% for example non-runner snake:
% for spks: import data using Analyses_Master with ops.do_tuningAnalysis = true and ops.tng.do_allTrials = true
% for dFF: import data using Analyses_Master with ops.do_tuningAnalysisDff = true and ops.tngn.do_allTrials = true
% use BullyBoy_20220214 as non-runner and Ao_20220402 as runner
% use Frida_20210409 (runner) as Fig1 example

inputType = 'tng'; %'tngn';
act = getActivityMeasure(p.(inputType),act_struct);
load([path.filepart_out,[inputType,'_all.mat']]);
this_struct = eval([inputType,'_all']);

save_root_fig = [path.root_summary,'figures\Fig1_fig\'];
save_root_png = [path.root_summary,'figures\Fig1_png\'];
save_root_pdf = [path.root_summary,'figures\Fig1_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig1_txt\'];

[prop,nf_binned,nft_binned] = preprocessActivityMeasure(act,p.(inputType),p,paq_beh,iscell);
bl_trial_frames_binned = [prop.trial_frames_binned(1,:)-[1:5]';prop.trial_frames_binned(1:10,:)];
nft_bl_binned = reshape(nf_binned(:,bl_trial_frames_binned),[],size(bl_trial_frames_binned,1),prop.numTrials);
nft_bl_avg_binned = nanmean(nft_bl_binned,3);
trials_all = createTrialsStruct(task,1:prop.numTrials_incl);
[traces_all,avgTraces_all,normAvgTraces_all] = createTracesStructs(nft_binned,trials_all);


%% Fig2_ExampleSnakePlots_NonRunner

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;

% pool data from A and X
these_refs_sort = 1; 
this_tracesStruct = traces_all; these_idcs_unsorted = {find(this_struct.passed.AW.Aonly_sel==1),find(this_struct.passed.AW.Xonly_sel==1)}; these_idcs_sorted = []; these_refs_sort = [1;4]; these_refs_norm = [];
these_avgTraces = {}; these_avgTraces{1} = avgTraces_all.A; these_avgTraces{2} = avgTraces_all.X;
this_data_A = these_avgTraces{1}(these_idcs_unsorted{1},:);
this_data_X = these_avgTraces{2}(these_idcs_unsorted{2},:);
this_data = [this_data_A;this_data_X];

% normalise data
this_normalisationData = this_data(:,p.general.bins_analysisWindow);
this_normalisationData_bl = this_data(:,p.general.bins_baselineWindow);
this_lower = nanmean(nft_bl_avg_binned([these_idcs_unsorted{1};these_idcs_unsorted{2}],:),2);
this_upper = nanmax(this_normalisationData,[],2);
this_data_norm = (this_data-this_lower) ./ (this_upper-this_lower);

% sort data
[~,temp] = nanmax(this_data_norm(:,p.general.bins_analysisWindow),[],2);
[~,this_order] = sort(temp);
this_order = flip(this_order);

% plot data
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; heatMap_task(this_data_norm(this_order,:),NaN,[],p,info,plt);
xlim([6,53])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off

% save plot
savefig(F,[save_root_fig,'\Fig2_ExampleSnakePlots_NonRunner.fig']);
saveas(F,[save_root_png,'\Fig2_ExampleSnakePlots_NonRunner.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleSnakePlots_NonRunner.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_ExampleSnakePlots_NonRunner.txt'],'wt');
fprintf(fid,['\nExampleSnakePlots_NonRunner\nn(seq)=',num2str(size(this_data_norm,1)),'\nn(cells)=',num2str(sum(prop.iscell==1)),...
    '\nn(seqA)=',num2str(length(these_idcs_unsorted{1})),'\nn(seqX)=',num2str(length(these_idcs_unsorted{2})),'\n']); fclose(fid);
fprintf(['\nExampleSnakePlots_NonRunner\nn(seq)=',num2str(size(this_data_norm,1)),'\nn(cells)=',num2str(sum(prop.iscell==1)),...
    '\nn(seqA)=',num2str(length(these_idcs_unsorted{1})),'\nn(seqX)=',num2str(length(these_idcs_unsorted{2})),'\n']);


%% Fig2_ExampleSnakePlots_Runner

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;

% pool data from A and X
these_refs_sort = 1; 
this_tracesStruct = traces_all; these_idcs_unsorted = {find(this_struct.passed.AW.Aonly_sel==1),find(this_struct.passed.AW.Xonly_sel==1)}; these_idcs_sorted = []; these_refs_sort = [1;4]; these_refs_norm = [];
these_avgTraces = {}; these_avgTraces{1} = avgTraces_all.A; these_avgTraces{2} = avgTraces_all.X;
this_data_A = these_avgTraces{1}(these_idcs_unsorted{1},:);
this_data_X = these_avgTraces{2}(these_idcs_unsorted{2},:);
this_data = [this_data_A;this_data_X];

% normalise data
this_normalisationData = this_data(:,p.general.bins_analysisWindow);
this_normalisationData_bl = this_data(:,p.general.bins_baselineWindow);
this_lower = nanmean(nft_bl_avg_binned([these_idcs_unsorted{1};these_idcs_unsorted{2}],:),2);
this_upper = nanmax(this_normalisationData,[],2);
this_data_norm = (this_data-this_lower) ./ (this_upper-this_lower);

% sort data
[~,temp] = nanmax(this_data_norm(:,p.general.bins_analysisWindow),[],2);
[~,this_order] = sort(temp);
this_order = flip(this_order);

% plot data
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; heatMap_task(this_data_norm(this_order,:),NaN,[],p,info,plt);
xlim([6,53])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off

% save plot
savefig(F,[save_root_fig,'\Fig2_ExampleSnakePlots_Runner.fig']);
saveas(F,[save_root_png,'\Fig2_ExampleSnakePlots_Runner.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleSnakePlots_Runner.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_ExampleSnakePlots_Runner.txt'],'wt');
fprintf(fid,['\nExampleSnakePlots_Runner\nn(seq)=',num2str(size(this_data_norm,1)),'\nn(cells)=',num2str(sum(prop.iscell==1)),...
    '\nn(seqA)=',num2str(length(these_idcs_unsorted{1})),'\nn(seqX)=',num2str(length(these_idcs_unsorted{2})),'\n']); fclose(fid);
fprintf(['\nExampleSnakePlots_Runner\nn(seq)=',num2str(size(this_data_norm,1)),'\nn(cells)=',num2str(sum(prop.iscell==1)),...
    '\nn(seqA)=',num2str(length(these_idcs_unsorted{1})),'\nn(seqX)=',num2str(length(these_idcs_unsorted{2})),'\n']);


%% Fig1_ExampleSnakePlots_1tile

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.333*34)]); hold on;

% pool data from A and X
these_refs_sort = 1; 
this_tracesStruct = traces_all; these_idcs_unsorted = {find(this_struct.passed.AW.Aonly_sel==1),find(this_struct.passed.AW.Xonly_sel==1)}; these_idcs_sorted = []; these_refs_sort = [1;4]; these_refs_norm = [];
these_avgTraces = {}; these_avgTraces{1} = avgTraces_all.A; these_avgTraces{2} = avgTraces_all.X;
this_data_A = these_avgTraces{1}(these_idcs_unsorted{1},:);
this_data_X = these_avgTraces{2}(these_idcs_unsorted{2},:);
this_data = [this_data_A;this_data_X];

% normalise data
this_normalisationData = this_data(:,p.general.bins_analysisWindow);
this_normalisationData_bl = this_data(:,p.general.bins_baselineWindow);
this_lower = nanmean(nft_bl_avg_binned([these_idcs_unsorted{1};these_idcs_unsorted{2}],:),2);
this_upper = nanmax(this_normalisationData,[],2);
this_data_norm = (this_data-this_lower) ./ (this_upper-this_lower);

% sort data
[~,temp] = nanmax(this_data_norm(:,p.general.bins_analysisWindow),[],2);
[~,this_order] = sort(temp);
this_order = flip(this_order);

% plot data
plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; heatMap_task(this_data_norm(this_order,:),NaN,[],p,info,plt);
xlim([6,53])
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'box','off');
set(gca,'xcolor','none');
set(gca,'ycolor','none');
box off

% save plot
savefig(F,[save_root_fig,'\Fig1_ExampleSnakePlots_1tile.fig']);
saveas(F,[save_root_png,'\Fig1_ExampleSnakePlots_1tile.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_ExampleSnakePlots_1tile.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig1_ExampleSnakePlots_1tile.txt'],'wt');
fprintf(fid,['\nExampleSnakePlots_1tile\nn(seq)=',num2str(size(this_data_norm,1)),'\nn(cells)=',num2str(sum(prop.iscell==1)),...
    '\nn(seqA)=',num2str(length(these_idcs_unsorted{1})),'\nn(seqX)=',num2str(length(these_idcs_unsorted{2})),'\n']); fclose(fid);
fprintf(['\nExampleSnakePlots_1tile\nn(seq)=',num2str(size(this_data_norm,1)),'\nn(cells)=',num2str(sum(prop.iscell==1)),...
    '\nn(seqA)=',num2str(length(these_idcs_unsorted{1})),'\nn(seqX)=',num2str(length(these_idcs_unsorted{2})),'\n']);


%% Fig1_ExampleSnakePlots_8tile

F = paper_figure([0,0.5,mm2inch(3*34),mm2inch(3*34)]); hold on;

this_lower = nanmean(nft_bl_avg_binned(:,2)); %%% THIS LOOKS LIKE AN ERROR
plt = struct(); plt.h = F; plt.split = 'odd_even'; plt.lower = this_lower; %{this_lower_A,this_lower_X};
F = snakePlots(traces_all,{find(this_struct.passed.AW.Aonly==1),find(this_struct.passed.AW.Xonly==1)},[],[1;3],[],p,info,plt);
% suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
%     char(p.(inputType).activityMeasure),'-smo',num2str(p.(inputType).smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
%     'passed=AXonly, all trials'])

for i=1:8
    subplot(2,4,i)
    xlim([6,53])
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    set(gca,'box','off');
    set(gca,'xcolor','none');
    set(gca,'ycolor','none');
    box off
end
% plot data
% plt = struct(); plt.clim = [0,1]; plt.lineWidth = 1; heatMap_task(this_data_norm(this_order,:),NaN,[],p,info,plt);
% xlim([6,53])
% set(gca,'xtick',[]);
% set(gca,'xticklabel',[]);
% set(gca,'ytick',[]);
% set(gca,'yticklabel',[]);
% set(gca,'box','off');
% set(gca,'xcolor','none');
% set(gca,'ycolor','none');
% box off

% save plot
savefig(F,[save_root_fig,'\Fig1_ExampleSnakePlots_8tile.fig']);
saveas(F,[save_root_png,'\Fig1_ExampleSnakePlots_8tile.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_ExampleSnakePlots_8tile.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig1_ExampleSnakePlots_8tile.txt'],'wt');
fprintf(fid,['\Fig1_ExampleSnakePlots_8tile\nn(seq)=',num2str(size(this_data_norm,1)),'\nn(cells)=',num2str(sum(prop.iscell==1)),...
    '\nn(seqA)=',num2str(length(these_idcs_unsorted{1})),'\nn(seqX)=',num2str(length(these_idcs_unsorted{2})),'\n']); fclose(fid);
fprintf(['\Fig1_ExampleSnakePlots_8tile\nn(seq)=',num2str(size(this_data_norm,1)),'\nn(cells)=',num2str(sum(prop.iscell==1)),...
    '\nn(seqA)=',num2str(length(these_idcs_unsorted{1})),'\nn(seqX)=',num2str(length(these_idcs_unsorted{2})),'\n']);


%% Colorbar

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

imagesc([0,0;1,1],[0,1])
colormap('parula');
temp=colorbar;
    
savefig(F,[save_root_fig,'\Fig1_ExampleSnakePlots_Colorbar.fig']);
saveas(F,[save_root_png,'\Fig1_ExampleSnakePlots_Colorbar.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_ExampleSnakePlots_Colorbar.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_ExampleInfoPlots_Runner

p = get_p; info = get_info;
path.root_summary = 'C:\SniffinHippo\Summary\';
save_root_fig = [path.root_summary,'figures\Fig2_fig\'];
save_root_png = [path.root_summary,'figures\Fig2_png\'];
save_root_pdf = [path.root_summary,'figures\Fig2_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig2_txt\'];

bcon = load('C:\SniffinHippo\Analysis\Ao\Ao_20220402\Ao_20220402_bcon.mat');
bcon = bcon.bcon;

nrows = 3; ncols = 1; m=0;
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.2*34)]); hold on;

m = m+1; subplot(nrows,ncols,m); hold on;
taskLines(p,info); xlabel(''); xticklabels({});
shadedErrorBar(1:p.general.numBins,nanmean(bcon.binwise_full.distance,1),nansem(bcon.binwise_full.distance,1),'lineProps',p.col.black)
xlim([6,53])
ylim([-200,400])
yticks([-200,0,400])
%ylabel('Distance (cm)')

m = m+1; subplot(nrows,ncols,m); hold on;
taskLines(p,info); xlabel(''); xticklabels({});
shadedErrorBar(1:p.general.numBins,nanmean(bcon.binwise_full.velocity,1),nansem(bcon.binwise_full.velocity,1),'lineProps',p.col.black)
xlim([6,53])
ylim([0,100])
yticks([0,100])
%ylabel('Velocity (cm/s)')

m = m+1; subplot(nrows,ncols,m); hold on;
taskLines(p,info);
shadedErrorBar(1:p.general.numBins,nanmean(bcon.binwise_full.acceleration,1),nansem(bcon.binwise_full.acceleration,1),'lineProps',p.col.black)
xlim([6,53])
ylim([-40,40])
yticks([-40,0,40])
%ylabel('Acceleration (cm/s2)')

% save plot
savefig(F,[save_root_fig,'\Fig2_ExampleInfoPlots_Runner.fig']);
saveas(F,[save_root_png,'\Fig2_ExampleInfoPlots_Runner.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleInfoPlots_Runner.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_ExampleInfoPlots_Nonrunner

p = get_p; info = get_info;
path.root_summary = 'C:\SniffinHippo\Summary\';
save_root_fig = [path.root_summary,'figures\Fig2_fig\'];
save_root_png = [path.root_summary,'figures\Fig2_png\'];
save_root_pdf = [path.root_summary,'figures\Fig2_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig2_txt\'];

bcon = load('C:\SniffinHippo\Analysis\BullyBoy\BullyBoy_20220214\BullyBoy_20220214_bcon.mat');
bcon = bcon.bcon;

nrows = 3; ncols = 1; m=0;
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.2*34)]); hold on;

m = m+1; subplot(nrows,ncols,m); hold on;
taskLines(p,info); xlabel(''); xticklabels({});
shadedErrorBar(1:p.general.numBins,nanmean(bcon.binwise_full.distance,1),nansem(bcon.binwise_full.distance,1),'lineProps',p.col.black)
xlim([6,53])
ylim([-200,400])
yticks([-200,0,400])
%ylabel('Distance (cm)')

m = m+1; subplot(nrows,ncols,m); hold on;
taskLines(p,info); xlabel(''); xticklabels({});
shadedErrorBar(1:p.general.numBins,nanmean(bcon.binwise_full.velocity,1),nansem(bcon.binwise_full.velocity,1),'lineProps',p.col.black)
xlim([6,53])
ylim([0,100])
yticks([0,100])
%ylabel('Velocity (cm/s)')

m = m+1; subplot(nrows,ncols,m); hold on;
taskLines(p,info);
shadedErrorBar(1:p.general.numBins,nanmean(bcon.binwise_full.acceleration,1),nansem(bcon.binwise_full.acceleration,1),'lineProps',p.col.black)
xlim([6,53])
ylim([-40,40])
yticks([-40,0,40])
%ylabel('Acceleration (cm/s2)')

% save plot
savefig(F,[save_root_fig,'\Fig2_ExampleInfoPlots_Nonrunner.fig']);
saveas(F,[save_root_png,'\Fig2_ExampleInfoPlots_Nonrunner.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleInfoPlots_Nonrunner.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_ExampleFitPlots_Runner

p = get_p; info = get_info;
path.root_summary = 'C:\SniffinHippo\Summary\';
save_root_fig = [path.root_summary,'figures\Fig2_fig\'];
save_root_png = [path.root_summary,'figures\Fig2_png\'];
save_root_pdf = [path.root_summary,'figures\Fig2_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig2_txt\'];

warp = load('F:\SniffinHippo\AnalysisX\Ao\Ao_20220402\Ao_20220402_warp_all.mat');
warp = warp.warp_all;

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(0.7*34)]); hold on;

this_x = 0:0.01:5;
this_y = warp.peak.exp1.a*exp(warp.peak.exp1.b*this_x);
plot(warp.peak.hist.edges(1:end-1),warp.peak.hist.prob,'k.')
plot(this_x,this_y,'k-')

xlim([0,5])
xticks([0:1:5])
xticklabels({'0','','','','','5'})
xlabel('Firing field peak (s)')
ylim([0,0.2])
yticks([0,0.2])
ylabel({'Proportion of','sequence cells'})

% save plot
savefig(F,[save_root_fig,'\Fig2_ExampleFitPlots_Runner.fig']);
saveas(F,[save_root_png,'\Fig2_ExampleFitPlots_Runner.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleFitPlots_Runner.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_ExampleFitPlots_Nonrunner

p = get_p; info = get_info;
path.root_summary = 'C:\SniffinHippo\Summary\';
save_root_fig = [path.root_summary,'figures\Fig2_fig\'];
save_root_png = [path.root_summary,'figures\Fig2_png\'];
save_root_pdf = [path.root_summary,'figures\Fig2_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig2_txt\'];

warp = load('F:\SniffinHippo\AnalysisX\BullyBoy\BullyBoy_20220214\BullyBoy_20220214_warp_all.mat');
warp = warp.warp_all;

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(0.7*34)]); hold on;

this_x = 0:0.01:5;
this_y = warp.peak.exp1.a*exp(warp.peak.exp1.b*this_x);
plot(warp.peak.hist.edges(1:end-1),warp.peak.hist.prob,'k.')
plot(this_x,this_y,'k-')

xlim([0,5])
xticks([0:1:5])
xticklabels({'0','','','','','5'})
xlabel('Firing field peak (s)')
ylim([0,0.2])
yticks([0,0.2])
ylabel({'Proportion of','sequence cells'})

% save plot
savefig(F,[save_root_fig,'\Fig2_ExampleFitPlots_Nonrunner.fig']);
saveas(F,[save_root_png,'\Fig2_ExampleFitPlots_Nonrunner.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleFitPlots_Nonrunner.pdf']); set(gcf,'Color',[1,1,1])

