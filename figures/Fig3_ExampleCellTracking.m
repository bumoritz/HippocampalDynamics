%% Fig3_ExampleCellTracking
% run Analyses_Master and tuningAnalysis for Elrond_20210601

path.root_summary = 'C:\SniffinHippo\Summary\';
save_root_fig = [path.root_summary,'figures\Fig3_fig\'];
save_root_png = [path.root_summary,'figures\Fig3_png\'];
save_root_pdf = [path.root_summary,'figures\Fig3_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig3_txt\'];


%% Preparations

% MODIFY TRIAL BLOCK IN plt.lower!!!

% Aonly, 100t, ktile
these_idcs = {};
for i=1:length(passed_100t)
    these_idcs{i} = find(passed_100t{i}.AW.Aonly==1);
end
plt = struct(); plt.split = 'none'; plt.normalisationRef = 'plotwise'; plt.sequences = "A"; plt.trialTypes = "A"; plt.lower = nanmean(nanmean(nft_bl_binned(:,:,701:800),3),2); %plt.lower = nanmean(nft_bl_avg_binned,2);
[F,these_normAvgTraces_A,these_idcs_sorted_A,plt] = snakePlots(traces_100t,these_idcs,[],[1:length(passed_100t)],[],p,info,plt);
suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
    this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
    'passed=Aonly, 100t, block-wise'])
savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Aonly_100t_ktile.fig']);
saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Aonly_100t_ktile.png']);

% Xonly, 100t, ktile
these_idcs = {};
for i=1:length(passed_100t)
    these_idcs{i} = find(passed_100t{i}.AW.Xonly==1);
end
plt = struct(); plt.split = 'none'; plt.normalisationRef = 'plotwise'; plt.sequences = "X"; plt.trialTypes = "X"; plt.lower = nanmean(nanmean(nft_bl_binned(:,:,701:800),3),2); %plt.lower = nanmean(nft_bl_avg_binned,2);
[F,these_normAvgTraces_X,these_idcs_sorted_X,plt] = snakePlots(traces_100t,these_idcs,[],[1:length(passed_100t)],[],p,info,plt);
suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
    this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
    'passed=Xonly, 100t, block-wise'])
savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Xonly_100t_ktile.fig']);
saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Xonly_100t_ktile.png']);


%% Fig3_ExampleSnakePlots_Elrond_d2_copy

%Aonlx and Xonly, ktile
F = paper_figure([0,0.5,mm2inch(10*34),mm2inch(3*34)]); hold on;
nrows = 1; ncols = length(passed_100t); numBlocks = length(passed_100t); these_refs_sort = [1:length(passed_100t)];
for r=1:nrows
    for c=1:ncols
        subplot(nrows,ncols,(r-1)*ncols+c)
        this_data_A = these_normAvgTraces_A{r,c}(these_idcs_sorted_A{c},:);
        this_data_X = these_normAvgTraces_X{r,c}(these_idcs_sorted_X{c},:);
        this_data = [this_data_A;this_data_X];

        [~,temp] = nanmax(this_data(:,p.general.bins_analysisWindow),[],2);
        [~,temp2] = sort(temp);
        
        if c==1
            this_sorting_1 = temp2;
        elseif c==8
            this_sorting_8 = temp2;
        end

        heatMap_task(this_data(temp2,:),NaN,[],p,info,plt); xlim([6,53]);
        set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]); set(gca,'box','off'); set(gca,'xcolor','none'); set(gca,'ycolor','none'); box off;
    end
end

% save plot
% savefig(F,[save_root_fig,'\Fig3_ExampleSnakePlots_Elrond_d2_copy-8.fig']);
% saveas(F,[save_root_png,'\Fig3_ExampleSnakePlots_Elrond_d2_copy-8.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_ExampleSnakePlots_Elrond_d2_copy-8.pdf']); set(gcf,'Color',[1,1,1])
% savefig(F,[save_root_fig,'\Fig3_ExampleSnakePlots_Elrond_d3_copy-1.fig']);
% saveas(F,[save_root_png,'\Fig3_ExampleSnakePlots_Elrond_d3_copy-1.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_ExampleSnakePlots_Elrond_d3_copy-1.pdf']); set(gcf,'Color',[1,1,1])


%% Preparations

% Aonly, 100t, ktile
these_idcs = {};
for i=1:length(passed_100t)
    these_idcs{i} = find(passed_100t{i}.AW.Aonly==1);
end
plt = struct(); plt.split = 'none'; plt.normalisationRef = 'rowwise_ref'; plt.sequences = "A"; plt.trialTypes = "A";    
[F,these_normAvgTraces_A,these_idcs_sorted_A,plt] = snakePlots(traces_100t,these_idcs,[],[1:length(passed_100t)],[1],p,info,plt);
suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
    this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
    'passed=Aonly, 100t, block-wise'])
savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Aonly_100t_ktile.fig']);
saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Aonly_100t_ktile.png']);

% Xonly, 100t, ktile
these_idcs = {};
for i=1:length(passed_100t)
    these_idcs{i} = find(passed_100t{i}.AW.Xonly==1);
end
plt = struct(); plt.split = 'none'; plt.normalisationRef = 'rowwise_ref'; plt.sequences = "X"; plt.trialTypes = "X";    
[F,these_normAvgTraces_X,these_idcs_sorted_X,plt] = snakePlots(traces_100t,these_idcs,[],[1:length(passed_100t)],[1],p,info,plt);
suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
    this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
    'passed=Xonly, 100t, block-wise'])
savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Xonly_100t_ktile.fig']);
saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Xonly_100t_ktile.png']);


%% Fig3_ExampleSnakePlots_Elrond_d2_sorted_1

%Aonlx and Xonly, ktile
F = paper_figure([0,0.5,mm2inch(10*34),mm2inch(3*34)]); hold on;
nrows = 1; ncols = length(passed_100t); numBlocks = length(passed_100t); these_refs_sort = [1:length(passed_100t)];
for r=1:nrows
    for c=1:ncols
        subplot(nrows,ncols,(r-1)*ncols+c)
        this_data_A = these_normAvgTraces_A{r,c}(these_idcs_sorted_A{ 1 },:);
        this_data_X = these_normAvgTraces_X{r,c}(these_idcs_sorted_X{ 1 },:);
        this_data = [this_data_A;this_data_X];

        heatMap_task(this_data(this_sorting_1,:),NaN,[],p,info,plt); xlim([6,53]);
        set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]); set(gca,'box','off'); set(gca,'xcolor','none'); set(gca,'ycolor','none'); box off;
    end
end

% save plot
% savefig(F,[save_root_fig,'\Fig3_ExampleSnakePlots_Elrond_d2_sorted_1.fig']);
% saveas(F,[save_root_png,'\Fig3_ExampleSnakePlots_Elrond_d2_sorted_1.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_ExampleSnakePlots_Elrond_d2_sorted_1.pdf']); set(gcf,'Color',[1,1,1])


%% Preparations

% Aonly, 100t, ktile
these_idcs = {};
for i=1:length(passed_100t)
    these_idcs{i} = find(passed_100t{i}.AW.Aonly==1);
end
plt = struct(); plt.split = 'none'; plt.normalisationRef = 'rowwise_ref'; plt.sequences = "A"; plt.trialTypes = "A";    
[F,these_normAvgTraces_A,these_idcs_sorted_A,plt] = snakePlots(traces_100t,these_idcs,[],[1:length(passed_100t)],[8],p,info,plt);
suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
    this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
    'passed=Aonly, 100t, block-wise'])
savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Aonly_100t_ktile.fig']);
saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Aonly_100t_ktile.png']);

% Xonly, 100t, ktile
these_idcs = {};
for i=1:length(passed_100t)
    these_idcs{i} = find(passed_100t{i}.AW.Xonly==1);
end
plt = struct(); plt.split = 'none'; plt.normalisationRef = 'rowwise_ref'; plt.sequences = "X"; plt.trialTypes = "X";    
[F,these_normAvgTraces_X,these_idcs_sorted_X,plt] = snakePlots(traces_100t,these_idcs,[],[1:length(passed_100t)],[8],p,info,plt);
suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
    this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
    'passed=Xonly, 100t, block-wise'])
savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Xonly_100t_ktile.fig']);
saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Xonly_100t_ktile.png']);


%% Fig3_ExampleSnakePlots_Elrond_d2_sorted_8

%Aonlx and Xonly, ktile
F = paper_figure([0,0.5,mm2inch(10*34),mm2inch(3*34)]); hold on;
nrows = 1; ncols = length(passed_100t); numBlocks = length(passed_100t); these_refs_sort = [1:length(passed_100t)];
for r=1:nrows
    for c=1:ncols
        subplot(nrows,ncols,(r-1)*ncols+c)
        this_data_A = these_normAvgTraces_A{r,c}(these_idcs_sorted_A{ 8 },:);
        this_data_X = these_normAvgTraces_X{r,c}(these_idcs_sorted_X{ 8 },:);
        this_data = [this_data_A;this_data_X];

        heatMap_task(this_data(this_sorting_8,:),NaN,[],p,info,plt); xlim([6,53]);
        set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]); set(gca,'box','off'); set(gca,'xcolor','none'); set(gca,'ycolor','none'); box off;
    end
end

% save plot
% savefig(F,[save_root_fig,'\Fig3_ExampleSnakePlots_Elrond_d2_sorted_8.fig']);
% saveas(F,[save_root_png,'\Fig3_ExampleSnakePlots_Elrond_d2_sorted_8.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_ExampleSnakePlots_Elrond_d2_sorted_8.pdf']); set(gcf,'Color',[1,1,1])


%% Preparations

% Aonly, 100t, ktile
these_idcs = {};
for i=1:length(passed_100t)
    these_idcs{i} = find(passed_100t{i}.AW.Aonly==1);
end
plt = struct(); plt.split = 'none'; plt.normalisationRef = 'rowwise_all'; plt.sequences = "A"; plt.trialTypes = "A"; plt.lower = nanmean(nft_bl_avg_binned,2);
[F,these_normAvgTraces_A,these_idcs_sorted_A,plt] = snakePlots(traces_100t,these_idcs,[],[1:length(passed_100t)],[],p,info,plt);
suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
    this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
    'passed=Aonly, 100t, block-wise'])
savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Aonly_100t_ktile.fig']);
saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Aonly_100t_ktile.png']);

% Xonly, 100t, ktile
these_idcs = {};
for i=1:length(passed_100t)
    these_idcs{i} = find(passed_100t{i}.AW.Xonly==1);
end
plt = struct(); plt.split = 'none'; plt.normalisationRef = 'rowwise_all'; plt.sequences = "X"; plt.trialTypes = "X"; plt.lower = nanmean(nft_bl_avg_binned,2);
[F,these_normAvgTraces_X,these_idcs_sorted_X,plt] = snakePlots(traces_100t,these_idcs,[],[1:length(passed_100t)],[],p,info,plt);
suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', ',...
    this_activity_measure,'-smo',num2str(p.tng.smoothingSd_preBinning),'-bin',num2str(p.general.binSize),', ',...
    'passed=Xonly, 100t, block-wise'])
savefig(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Xonly_100t_ktile.fig']);
saveas(F,[path.filepart_outX,'plots/SnakePlots/',info.animal,'_',info.date,'_',this_activity_measure,'_','snakePlot_Xonly_100t_ktile.png']);

% pool data from A and X
this_struct = tng_all;
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
%this_order = flip(this_order);
    

%% Fig3_ExampleSnakePlots_Elrond_d2_sorted_full

%Aonlx and Xonly, ktile
F = paper_figure([0,0.5,mm2inch(10*34),mm2inch(3*34)]); hold on;
nrows = 1; ncols = length(passed_100t); numBlocks = length(passed_100t); these_refs_sort = [1:length(passed_100t)];
for r=1:nrows
    for c=1:ncols
        subplot(nrows,ncols,(r-1)*ncols+c);
        this_data_A = these_normAvgTraces_A{r,c}(these_idcs_unsorted{1},:);
        this_data_X = these_normAvgTraces_X{r,c}(these_idcs_unsorted{2},:);
        this_data = [this_data_A;this_data_X];

        heatMap_task(this_data(this_order,:),NaN,[],p,info,plt); xlim([6,53]); hold on;
        set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]); set(gca,'box','off'); set(gca,'xcolor','none'); set(gca,'ycolor','none'); box off;
    end
end

% save plot
savefig(F,[save_root_fig,'\Fig3_ExampleSnakePlots_Elrond_d2_sorted_full.fig']);
saveas(F,[save_root_png,'\Fig3_ExampleSnakePlots_Elrond_d2_sorted_full.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_ExampleSnakePlots_Elrond_d2_sorted_full.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_ExampleSnakePlots_Elrond_d2_sorted_full_peaks

%Aonlx and Xonly, ktile
F = paper_figure([0,0.5,mm2inch(10*34),mm2inch(3*34)]); hold on;
nrows = 1; ncols = length(passed_100t); numBlocks = length(passed_100t); these_refs_sort = [1:length(passed_100t)];
for r=1:nrows
    for c=1:ncols
        subplot(nrows,ncols,(r-1)*ncols+c);
        this_data_A = these_normAvgTraces_A{r,c}(these_idcs_unsorted{1},:);
        this_data_X = these_normAvgTraces_X{r,c}(these_idcs_unsorted{2},:);
        this_data = [this_data_A;this_data_X];

        heatMap_task(this_data(this_order,:),NaN,[],p,info,plt); xlim([6,53]); hold on;
        set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]); set(gca,'box','off'); set(gca,'xcolor','none'); set(gca,'ycolor','none'); box off;

        % highlight peak times
        this_data_sorted = this_data(this_order,p.general.bins_analysisWindow);
        [~,these_peakTimes] = nanmax(this_data_sorted,[],2);
        h=plot(these_peakTimes+p.general.bins_analysisWindow(1)-1,1:length(these_peakTimes),'w.','MarkerSize',5);
    end
end

% save plot
savefig(F,[save_root_fig,'\Fig3_ExampleSnakePlots_Elrond_d2_sorted_full_peaks.fig']);
saveas(F,[save_root_png,'\Fig3_ExampleSnakePlots_Elrond_d2_sorted_full_peaks.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_ExampleSnakePlots_Elrond_d2_sorted_full_peaks.pdf']); set(gcf,'Color',[1,1,1])






