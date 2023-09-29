%% Fig2_SequenceAnalysis_OUT

p = get_p; info = get_info;


%% --- Main figure ---

%% Fig2_f_SequenceAnalysis_SeqCellProfile

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

this_data_nonrunners = propByPeakTime_perSeqCells.pref(d_info.running(:,1)==0,:);
this_data_runners = propByPeakTime_perSeqCells.pref(d_info.running(:,1)==1,:);
shadedErrorBar(1:size(this_data_nonrunners,2),nanmean(this_data_nonrunners',2)*100,nansem(this_data_nonrunners',2)*100,'lineProps',p.col.nonrunner);
shadedErrorBar(1:size(this_data_runners,2),nanmean(this_data_runners',2)*100,nansem(this_data_runners',2)*100,'lineProps',p.col.runner);
xline(5*5/3,'k:'); xline(2*5*5/3,'k:');
xlim([1,size(this_data_nonrunners,2)+1]-0.5)
xticks([0:5:length(binEdges)]+0.5)
xticklabels({'0','1','2','3','4','5'})
ytickformat('percentage')
ylim([0,15])
yticks([0:5:15])
xlabel('Firing field peak (s)')
ylabel('Proportion of sequence cells')

savefig(F,[save_root_fig,'\Fig2_f_SequenceAnalysis_SeqCellProfile.fig']);
saveas(F,[save_root_png,'\Fig2_f_SequenceAnalysis_SeqCellProfile.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_SequenceAnalysis_SeqCellProfile.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_SequenceAnalysis_SeqCellProfile.txt'],'wt');


%% Fig2_f_SequenceAnalysis_SeqCellProfile_abs

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

this_data_nonrunners = propByPeakTime_perIscell.pref(d_info.running(:,1)==0,:);
this_data_runners = propByPeakTime_perIscell.pref(d_info.running(:,1)==1,:);
shadedErrorBar(1:size(this_data_nonrunners,2),nanmean(this_data_nonrunners',2)*100,nansem(this_data_nonrunners',2)*100,'lineProps',p.col.nonrunner);
shadedErrorBar(1:size(this_data_runners,2),nanmean(this_data_runners',2)*100,nansem(this_data_runners',2)*100,'lineProps',p.col.runner);
xline(5*5/3,'k:'); xline(2*5*5/3,'k:');
xlim([1,size(this_data_nonrunners,2)+1]-0.5)
xticks([0:5:length(binEdges)]+0.5)
xticklabels({'0','1','2','3','4','5'})
ytickformat('percentage')
ylim([0,15])
yticks([0:5:15])
xlabel('Firing field peak (s)')
ylabel('Proportion of cells')

savefig(F,[save_root_fig,'\Fig2_f_SequenceAnalysis_SeqCellProfile_abs.fig']);
saveas(F,[save_root_png,'\Fig2_f_SequenceAnalysis_SeqCellProfile_abs.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_SequenceAnalysis_SeqCellProfile_abs.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_SequenceAnalysis_SeqCellProfile_abs.txt'],'wt');


%% Fig2_f_SequenceAnalysis_SequenceCells

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,8);
this_data_nonrunners = propByCellType_perSeqCells.pref(d_info.running(:,1)==0,1)*100;
this_data_runners = propByCellType_perSeqCells.pref(d_info.running(:,1)==1,1)*100;
this_data(1:length(this_data_nonrunners),1) = this_data_nonrunners;
this_data(1:length(this_data_runners),2) = this_data_runners;
this_data_nonrunners = propByCellType_perSeqCells.pref(d_info.running(:,1)==0,2)*100;
this_data_runners = propByCellType_perSeqCells.pref(d_info.running(:,1)==1,2)*100;
this_data(1:length(this_data_nonrunners),4) = this_data_nonrunners;
this_data(1:length(this_data_runners),5) = this_data_runners;
this_data_nonrunners = propByCellType_perSeqCells.pref(d_info.running(:,1)==0,3)*100;
this_data_runners = propByCellType_perSeqCells.pref(d_info.running(:,1)==1,3)*100;
this_data(1:length(this_data_nonrunners),7) = this_data_nonrunners;
this_data(1:length(this_data_runners),8) = this_data_runners;

v = violinplot(this_data,labels);
for i=1:8
    v(i).BoxColor = 'k'; v(i).ViolinPlot.EdgeColor = 'none'; v(i).ScatterPlot.Marker = '.'; v(i).ScatterPlot.SizeData = 100; v(i).BoxPlot.LineWidth = 1; v(i).WhiskerPlot.Color = 'none'; v(i).MedianPlot.SizeData = 15;
    if mod(i,3)==0
        v(i).ViolinColor = p.col.white; v(i).ScatterPlot.MarkerEdgeColor =  p.col.white;
    elseif mod(i,3)==1
        v(i).ViolinColor = p.col.nonrunner; v(i).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner;
    elseif mod(i,3)==2
        v(i).ViolinColor = p.col.runner; v(i).ScatterPlot.MarkerEdgeColor =  p.col.runner;
    end
end

xlim([0,9])
xticks([1.5,4.5,7.5])
xticklabels({'Early','Intermediate','Late'})
xlabel('Firing field peak')

ytickformat('percentage')
ylim([0,100])
yticks([0:50:100])
ylabel({'Proportion of sequence cells'});

savefig(F,[save_root_fig,'\Fig2_f_SequenceAnalysis_SequenceCells.fig']);
saveas(F,[save_root_png,'\Fig2_f_SequenceAnalysis_SequenceCells.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_SequenceAnalysis_SequenceCells.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_SequenceAnalysis_SequenceCells.txt'],'wt');


%% Fig2_f_SequenceAnalysis_SequenceCells_abs

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,8);
this_data_nonrunners = propByCellType_perIscell.pref(d_info.running(:,1)==0,1)*100;
this_data_runners = propByCellType_perIscell.pref(d_info.running(:,1)==1,1)*100;
this_data(1:length(this_data_nonrunners),1) = this_data_nonrunners;
this_data(1:length(this_data_runners),2) = this_data_runners;
this_data_nonrunners = propByCellType_perIscell.pref(d_info.running(:,1)==0,2)*100;
this_data_runners = propByCellType_perIscell.pref(d_info.running(:,1)==1,2)*100;
this_data(1:length(this_data_nonrunners),4) = this_data_nonrunners;
this_data(1:length(this_data_runners),5) = this_data_runners;
this_data_nonrunners = propByCellType_perIscell.pref(d_info.running(:,1)==0,3)*100;
this_data_runners = propByCellType_perIscell.pref(d_info.running(:,1)==1,3)*100;
this_data(1:length(this_data_nonrunners),7) = this_data_nonrunners;
this_data(1:length(this_data_runners),8) = this_data_runners;

v = violinplot(this_data,labels);
for i=1:8
    v(i).BoxColor = 'k'; v(i).ViolinPlot.EdgeColor = 'none'; v(i).ScatterPlot.Marker = '.'; v(i).ScatterPlot.SizeData = 100; v(i).BoxPlot.LineWidth = 1; v(i).WhiskerPlot.Color = 'none'; v(i).MedianPlot.SizeData = 15;
    if mod(i,3)==0
        v(i).ViolinColor = p.col.white; v(i).ScatterPlot.MarkerEdgeColor =  p.col.white;
    elseif mod(i,3)==1
        v(i).ViolinColor = p.col.nonrunner; v(i).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner;
    elseif mod(i,3)==2
        v(i).ViolinColor = p.col.runner; v(i).ScatterPlot.MarkerEdgeColor =  p.col.runner;
    end
end

xlim([0,9])
xticks([1.5,4.5,7.5])
xticklabels({'Early','Intermediate','Late'})
xlabel('Firing field peak')

ytickformat('percentage')
ylim([0,100])
yticks([0:50:100])
ylabel({'Proportion of cells'});

savefig(F,[save_root_fig,'\Fig2_f_SequenceAnalysis_SequenceCells_abs.fig']);
saveas(F,[save_root_png,'\Fig2_f_SequenceAnalysis_SequenceCells_abs.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_SequenceAnalysis_SequenceCells_abs.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_SequenceAnalysis_SequenceCells_abs.txt'],'wt');


%% Fig2_f_SequenceAnalysis_SequenceCells_early

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,2);
this_data_nonrunners = propByCellType_perSeqCells.pref(d_info.running(:,1)==0,1)*100;
this_data_runners = propByCellType_perSeqCells.pref(d_info.running(:,1)==1,1)*100;
this_data(1:length(this_data_nonrunners),1) = this_data_nonrunners;
this_data(1:length(this_data_runners),2) = this_data_runners;

v = violinplot(this_data,labels); 
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner; v(2).ScatterPlot.MarkerEdgeColor = p.col.runner; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15;
ytickformat('percentage')
ylim([0,100])
yticks([0:50:100])
ylabel({'Proportion of sequence cells'});

[pval,~,stats] = ranksum(this_data_nonrunners,this_data_runners)

savefig(F,[save_root_fig,'\Fig2_f_SequenceAnalysis_SequenceCells_early.fig']);
saveas(F,[save_root_png,'\Fig2_f_SequenceAnalysis_SequenceCells_early.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_SequenceAnalysis_SequenceCells_early.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_SequenceAnalysis_SequenceCells_early.txt'],'wt');


%% Fig2_f_SequenceAnalysis_SequenceCells_intermediate

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,2);
this_data_nonrunners = propByCellType_perSeqCells.pref(d_info.running(:,1)==0,2)*100;
this_data_runners = propByCellType_perSeqCells.pref(d_info.running(:,1)==1,2)*100;
this_data(1:length(this_data_nonrunners),1) = this_data_nonrunners;
this_data(1:length(this_data_runners),2) = this_data_runners;

v = violinplot(this_data,labels); 
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner; v(2).ScatterPlot.MarkerEdgeColor = p.col.runner; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15;
ytickformat('percentage')
ylim([0,100])
yticks([0:50:100])
ylabel({'Proportion of sequence cells'});

[pval,~,stats] = ranksum(this_data_nonrunners,this_data_runners)

savefig(F,[save_root_fig,'\Fig2_f_SequenceAnalysis_SequenceCells_intermediate.fig']);
saveas(F,[save_root_png,'\Fig2_f_SequenceAnalysis_SequenceCells_intermediate.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_SequenceAnalysis_SequenceCells_intermediate.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_SequenceAnalysis_SequenceCells_intermediate.txt'],'wt');


%% Fig2_f_SequenceAnalysis_SequenceCells_late

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,2);
this_data_nonrunners = propByCellType_perSeqCells.pref(d_info.running(:,1)==0,3)*100;
this_data_runners = propByCellType_perSeqCells.pref(d_info.running(:,1)==1,3)*100;
this_data(1:length(this_data_nonrunners),1) = this_data_nonrunners;
this_data(1:length(this_data_runners),2) = this_data_runners;

v = violinplot(this_data,labels); 
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner; v(2).ScatterPlot.MarkerEdgeColor = p.col.runner; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15;
ytickformat('percentage')
ylim([0,100])
yticks([0:50:100])
ylabel({'Proportion of sequence cells'});

[pval,~,stats] = ranksum(this_data_nonrunners,this_data_runners)

savefig(F,[save_root_fig,'\Fig2_f_SequenceAnalysis_SequenceCells_late.fig']);
saveas(F,[save_root_png,'\Fig2_f_SequenceAnalysis_SequenceCells_late.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_SequenceAnalysis_SequenceCells_late.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_SequenceAnalysis_SequenceCells_late.txt'],'wt');


%% Fig2_f_SequenceAnalysis_ExpB

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,2);
this_data_nonrunners = exp_b(d_info.running(:,1)==0);
this_data_runners = exp_b(d_info.running(:,1)==1);
this_data(1:length(this_data_nonrunners),1) = this_data_nonrunners;
this_data(1:length(this_data_runners),2) = this_data_runners;

v = violinplot(this_data,labels); 
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner; v(2).ScatterPlot.MarkerEdgeColor = p.col.runner; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15;
ylim([-1.5,0])
yticks([-1.5:0.5:0])
ylabel({'Sequence shape','(decay rate, s-1)'});

[pval,~,stats] = ranksum(this_data_nonrunners,this_data_runners)

savefig(F,[save_root_fig,'\Fig2_f_SequenceAnalysis_ExpB.fig']);
saveas(F,[save_root_png,'\Fig2_f_SequenceAnalysis_ExpB.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_SequenceAnalysis_ExpB.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_SequenceAnalysis_ExpB.txt'],'wt');


%% Fig2_f_SequenceAnalysis_Selectivity

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,2);
this_data_nonrunners = selectivity.pref(d_info.running(:,1)==0);
this_data_runners = selectivity.pref(d_info.running(:,1)==1);
this_data(1:length(this_data_nonrunners),1) = this_data_nonrunners;
this_data(1:length(this_data_runners),2) = this_data_runners;

v = violinplot(this_data,labels); 
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner; v(2).ScatterPlot.MarkerEdgeColor = p.col.runner; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15;
ylim([0,0.3])
yticks([0:0.1:0.3])
ylabel({'Selectivity'});

[pval,~,stats] = ranksum(this_data_nonrunners,this_data_runners)

savefig(F,[save_root_fig,'\Fig2_f_SequenceAnalysis_Selectivity.fig']);
saveas(F,[save_root_png,'\Fig2_f_SequenceAnalysis_Selectivity.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_SequenceAnalysis_Selectivity.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_SequenceAnalysis_Selectivity.txt'],'wt');


%% Fig2_f_DecodingAnalysis_typeDecodingCorrect

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

this_data_nonrunners = typeDecodingCorrect(d_info.running(:,1)==0,:);
this_data_runners = typeDecodingCorrect(d_info.running(:,1)==1,:);
shadedErrorBar(1:size(this_data_nonrunners,2),nanmean(this_data_nonrunners',2)*100,nansem(this_data_nonrunners',2)*100,'lineProps',p.col.nonrunner);
shadedErrorBar(1:size(this_data_runners,2),nanmean(this_data_runners',2)*100,nansem(this_data_runners',2)*100,'lineProps',p.col.runner);

xline(5*5/3,'k:'); xline(2*5*5/3,'k:');
xlim([0,25]); xticks([0:5:25]); xticklabels({'0','1','2','3','4','5'});
xlabel('Time (s)')
ytickformat('percentage')
ylim([50,100]); yticks([50,75,100]);
ylabel('1st Odor decoding accuracy')

savefig(F,[save_root_fig,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect.fig']);
saveas(F,[save_root_png,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect.txt'],'wt');


%% Fig2_f_DecodingAnalysis_typeDecodingCorrect_EIL

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,8);
this_data_nonrunners = nanmean(typeDecodingCorrect(d_info.running(:,1)==0,1:8),2)*100; %0:5/25:5 split by 0,1.66,3.33,5
this_data_runners = nanmean(typeDecodingCorrect(d_info.running(:,1)==1,1:8),2)*100;
this_data(1:length(this_data_nonrunners),1) = this_data_nonrunners;
this_data(1:length(this_data_runners),2) = this_data_runners;
this_data_nonrunners = nanmean(typeDecodingCorrect(d_info.running(:,1)==0,9:16),2)*100; %0:5/25:5 split by 0,1.66,3.33,5
this_data_runners = nanmean(typeDecodingCorrect(d_info.running(:,1)==1,9:16),2)*100;
this_data(1:length(this_data_nonrunners),4) = this_data_nonrunners;
this_data(1:length(this_data_runners),5) = this_data_runners;
this_data_nonrunners = nanmean(typeDecodingCorrect(d_info.running(:,1)==0,17:24),2)*100; %0:5/25:5 split by 0,1.66,3.33,5
this_data_runners = nanmean(typeDecodingCorrect(d_info.running(:,1)==1,17:24),2)*100;
this_data(1:length(this_data_nonrunners),7) = this_data_nonrunners;
this_data(1:length(this_data_runners),8) = this_data_runners;

v = violinplot(this_data,labels);
for i=1:8
    v(i).BoxColor = 'k'; v(i).ViolinPlot.EdgeColor = 'none'; v(i).ScatterPlot.Marker = '.'; v(i).ScatterPlot.SizeData = 100; v(i).BoxPlot.LineWidth = 1; v(i).WhiskerPlot.Color = 'none'; v(i).MedianPlot.SizeData = 15;
    if mod(i,3)==0
        v(i).ViolinColor = p.col.white; v(i).ScatterPlot.MarkerEdgeColor =  p.col.white;
    elseif mod(i,3)==1
        v(i).ViolinColor = p.col.nonrunner; v(i).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner;
    elseif mod(i,3)==2
        v(i).ViolinColor = p.col.runner; v(i).ScatterPlot.MarkerEdgeColor =  p.col.runner;
    end
end

xlim([0,9])
xticks([1.5,4.5,7.5])
xticklabels({'Early','Intermediate','Late'})
xlabel('Time')

ytickformat('percentage')
ylim([50,100]); yticks([50,75,100]);
ylabel('1st Odor decoding accuracy')

savefig(F,[save_root_fig,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_EIL.fig']);
saveas(F,[save_root_png,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_EIL.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_EIL.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_EIL.txt'],'wt');



%% Fig2_f_DecodingAnalysis_typeDecodingCorrect_early

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,2);
this_data_nonrunners = nanmean(typeDecodingCorrect(d_info.running(:,1)==0,1:8),2)*100; %0:5/25:5 split by 0,1.66,3.33,5
this_data_runners = nanmean(typeDecodingCorrect(d_info.running(:,1)==1,1:8),2)*100;
this_data(1:length(this_data_nonrunners),1) = this_data_nonrunners;
this_data(1:length(this_data_runners),2) = this_data_runners;

v = violinplot(this_data,labels); 
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner; v(2).ScatterPlot.MarkerEdgeColor = p.col.runner; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15;
ytickformat('percentage')
ylim([50,100]); yticks([50,75,100]);
ylabel('1st odor decoding accuracy')

[pval,~,stats] = ranksum(this_data_nonrunners,this_data_runners)

savefig(F,[save_root_fig,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_early.fig']);
saveas(F,[save_root_png,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_early.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_early.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_early.txt'],'wt');


%% Fig2_f_DecodingAnalysis_typeDecodingCorrect_middle

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,2);
this_data_nonrunners = nanmean(typeDecodingCorrect(d_info.running(:,1)==0,9:16),2)*100; %0:5/25:5 split by 0,1.66,3.33,5
this_data_runners = nanmean(typeDecodingCorrect(d_info.running(:,1)==1,9:16),2)*100;
this_data(1:length(this_data_nonrunners),1) = this_data_nonrunners;
this_data(1:length(this_data_runners),2) = this_data_runners;

v = violinplot(this_data,labels); 
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner; v(2).ScatterPlot.MarkerEdgeColor = p.col.runner; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15;
ytickformat('percentage')
ylim([50,100]); yticks([50,75,100]);
ylabel('1st odor decoding accuracy')

[pval,~,stats] = ranksum(this_data_nonrunners,this_data_runners)

savefig(F,[save_root_fig,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_middle.fig']);
saveas(F,[save_root_png,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_middle.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_middle.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_middle.txt'],'wt');


%% Fig2_f_DecodingAnalysis_typeDecodingCorrect_late

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,2);
this_data_nonrunners = nanmean(typeDecodingCorrect(d_info.running(:,1)==0,17:24),2)*100; %0:5/25:5 split by 0,1.66,3.33,5
this_data_runners = nanmean(typeDecodingCorrect(d_info.running(:,1)==1,17:24),2)*100;
this_data(1:length(this_data_nonrunners),1) = this_data_nonrunners;
this_data(1:length(this_data_runners),2) = this_data_runners;

v = violinplot(this_data,labels); 
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner; v(2).ScatterPlot.MarkerEdgeColor = p.col.runner; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15;
ytickformat('percentage')
ylim([50,100]); yticks([50,75,100]);
ylabel('1st odor decoding accuracy')

[pval,~,stats] = ranksum(this_data_nonrunners,this_data_runners)

savefig(F,[save_root_fig,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_late.fig']);
saveas(F,[save_root_png,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_late.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_late.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_DecodingAnalysis_typeDecodingCorrect_late.txt'],'wt');


%% --- Supplementary figure ---

%% Fig2_f_SequenceAnalysis_NumCells

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,2);
this_data_nonrunners = numCells(d_info.running(:,1)==0);
this_data_runners = numCells(d_info.running(:,1)==1);
this_data(1:length(this_data_nonrunners),1) = this_data_nonrunners;
this_data(1:length(this_data_runners),2) = this_data_runners;

v = violinplot(this_data,labels); 
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner; v(2).ScatterPlot.MarkerEdgeColor = p.col.runner; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15;
ylim([0,1500])
yticks([0,500,1000,1500])
ylabel(['Number of cells']);

[pval,~,stats] = ranksum(this_data_nonrunners,this_data_runners)

savefig(F,[save_root_fig,'\Fig2_f_SequenceAnalysis_NumCells.fig']);
saveas(F,[save_root_png,'\Fig2_f_SequenceAnalysis_NumCells.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_SequenceAnalysis_NumCells.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_SequenceAnalysis_NumCells.txt'],'wt');


%% Fig2_f_SequenceAnalysis_PropSeqCells

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,2);
this_data_nonrunners = propSeqCells.pref(d_info.running(:,1)==0)*100;
this_data_runners = propSeqCells.pref(d_info.running(:,1)==1)*100;
this_data(1:length(this_data_nonrunners),1) = this_data_nonrunners;
this_data(1:length(this_data_runners),2) = this_data_runners;

v = violinplot(this_data,labels); 
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner; v(2).ScatterPlot.MarkerEdgeColor = p.col.runner; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15;
ytickformat('percentage')
ylim([0,50])
yticks([0:25:50])
ylabel({'Sequence participation','probability'});

[pval,~,stats] = ranksum(this_data_nonrunners,this_data_runners)

savefig(F,[save_root_fig,'\Fig2_f_SequenceAnalysis_PropSeqCells.fig']);
saveas(F,[save_root_png,'\Fig2_f_SequenceAnalysis_PropSeqCells.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_SequenceAnalysis_PropSeqCells.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_SequenceAnalysis_PropSeqCells.txt'],'wt');


%% Fig2_f_SequenceAnalysis_ExpA

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data = nan(d_info.numAnimals,2);
this_data_nonrunners = exp_a(d_info.running(:,1)==0);
this_data_runners = exp_a(d_info.running(:,1)==1);
this_data(1:length(this_data_nonrunners),1) = this_data_nonrunners;
this_data(1:length(this_data_runners),2) = this_data_runners;

v = violinplot(this_data,labels); 
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(1).ViolinPlot.EdgeColor = 'none'; v(2).ViolinPlot.EdgeColor = 'none';
v(1).ScatterPlot.Marker = '.'; v(2).ScatterPlot.Marker = '.'; v(1).ScatterPlot.MarkerEdgeColor =  p.col.nonrunner; v(2).ScatterPlot.MarkerEdgeColor = p.col.runner; v(1).ScatterPlot.SizeData = 100; v(2).ScatterPlot.SizeData = 100;
v(1).BoxPlot.LineWidth = 1; v(2).BoxPlot.LineWidth = 1; v(1).WhiskerPlot.Color = 'none'; v(2).WhiskerPlot.Color = 'none';
v(1).MedianPlot.SizeData = 15; v(2).MedianPlot.SizeData = 15;
ylim([0,0.3])
yticks([0:0.1:0.3])
ylabel({'Sequence shape','(initial value)'});

[pval,~,stats] = ranksum(this_data_nonrunners,this_data_runners)

savefig(F,[save_root_fig,'\Fig2_f_SequenceAnalysis_ExpA.fig']);
saveas(F,[save_root_png,'\Fig2_f_SequenceAnalysis_ExpA.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_f_SequenceAnalysis_ExpA.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_f_SequenceAnalysis_ExpA.txt'],'wt');


%% --- In preparation ---

%% Fig2_Bcon_Distance

% trigger misalignment for animal 9
distance(9,:) = NaN;

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;

taskLines(p,info);
yline(0,'k-');
shadedErrorBar(1:p.general.numBins,nanmean(distance(find(d_info.running(:,1)==0),:),1),nansem(distance(find(d_info.running(:,1)==0),:),1),'lineProps',p.col.nonrunner)
shadedErrorBar(1:p.general.numBins,nanmean(distance(find(d_info.running(:,1)==1),:),1),nansem(distance(find(d_info.running(:,1)==1),:),1),'lineProps',p.col.runner)
ylim([-400,800])
yticks([-400:400:800])
ylabel('Distance (cm)')

savefig(F,[save_root_fig,'\Fig2_Bcon_Distance.fig']);
saveas(F,[save_root_png,'\Fig2_Bcon_Distance.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_Bcon_Distance.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_Bcon_Distance_AX

% trigger misalignment for animal 9
distance_X(9,:) = NaN;
distance_A(9,:) = NaN;

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;

taskLines(p,info);
yline(0,'k-');
shadedErrorBar(1:p.general.numBins,nanmean(distance_X(find(d_info.running(:,1)==1),:),1),nansem(distance_X(find(d_info.running(:,1)==1),:),1),'lineProps',p.col.X)
shadedErrorBar(1:p.general.numBins,nanmean(distance_A(find(d_info.running(:,1)==1),:),1),nansem(distance_A(find(d_info.running(:,1)==1),:),1),'lineProps',p.col.A)
ylim([-400,800])
yticks([-400:400:800])
ylabel('Distance (cm)')

savefig(F,[save_root_fig,'\Fig2_Bcon_Distance_AX.fig']);
saveas(F,[save_root_png,'\Fig2_Bcon_Distance_AX.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_Bcon_Distance_AX.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_Bcon_Velocity

% trigger misalignment for animal 9
velocity(9,:) = NaN;

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;

taskLines(p,info);
shadedErrorBar(1:p.general.numBins,nanmean(velocity(find(d_info.running(:,1)==0),:),1),nansem(velocity(find(d_info.running(:,1)==0),:),1),'lineProps',p.col.nonrunner)
shadedErrorBar(1:p.general.numBins,nanmean(velocity(find(d_info.running(:,1)==1),:),1),nansem(velocity(find(d_info.running(:,1)==1),:),1),'lineProps',p.col.runner)
ylim([0,120])
yticks([0:40:120])
ylabel('Velocity (cm/s)')

savefig(F,[save_root_fig,'\Fig2_Bcon_Velocity.fig']);
saveas(F,[save_root_png,'\Fig2_Bcon_Velocity.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_Bcon_Velocity.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_Bcon_Velocity_AX

% trigger misalignment for animal 9
velocity_X(9,:) = NaN;
velocity_A(9,:) = NaN;

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;

taskLines(p,info);
shadedErrorBar(1:p.general.numBins,nanmean(velocity_X(find(d_info.running(:,1)==1),:),1),nansem(velocity_X(find(d_info.running(:,1)==1),:),1),'lineProps',p.col.X)
shadedErrorBar(1:p.general.numBins,nanmean(velocity_A(find(d_info.running(:,1)==1),:),1),nansem(velocity_A(find(d_info.running(:,1)==1),:),1),'lineProps',p.col.A)
ylim([0,120])
yticks([0:40:120])
ylabel('Velocity (cm/s)')

savefig(F,[save_root_fig,'\Fig2_Bcon_Velocity_AX.fig']);
saveas(F,[save_root_png,'\Fig2_Bcon_Velocity_AX.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_Bcon_Velocity_AX.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_Bcon_Acceleration

% trigger misalignment for animal 9
acceleration(9,:) = NaN;

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;

taskLines(p,info);
yline(0,'k-');
shadedErrorBar(1:p.general.numBins,nanmean(acceleration(find(d_info.running(:,1)==0),:),1),nansem(acceleration(find(d_info.running(:,1)==0),:),1),'lineProps',p.col.nonrunner)
shadedErrorBar(1:p.general.numBins,nanmean(acceleration(find(d_info.running(:,1)==1),:),1),nansem(acceleration(find(d_info.running(:,1)==1),:),1),'lineProps',p.col.runner)
ylim([-80,80])
yticks([-80:40:80])
ylabel('Acceleration (cm/s2)')

savefig(F,[save_root_fig,'\Fig2_Bcon_Acceleration.fig']);
saveas(F,[save_root_png,'\Fig2_Bcon_Acceleration.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_Bcon_Acceleration.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_Bcon_Acceleration_AX

% trigger misalignment for animal 9
acceleration_X(9,:) = NaN;
acceleration_A(9,:) = NaN;

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;

taskLines(p,info);
yline(0,'k-');
shadedErrorBar(1:p.general.numBins,nanmean(acceleration_X(find(d_info.running(:,1)==1),:),1),nansem(acceleration_X(find(d_info.running(:,1)==1),:),1),'lineProps',p.col.X)
shadedErrorBar(1:p.general.numBins,nanmean(acceleration_A(find(d_info.running(:,1)==1),:),1),nansem(acceleration_A(find(d_info.running(:,1)==1),:),1),'lineProps',p.col.A)
ylim([-80,80])
yticks([-80:40:80])
ylabel('Acceleration (cm/s2)')

savefig(F,[save_root_fig,'\Fig2_Bcon_Acceleration_AX.fig']);
saveas(F,[save_root_png,'\Fig2_Bcon_Acceleration_AX.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_Bcon_Acceleration_AX.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_SequenceAnalysis_PlaceCellProfile

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

%this_data_nonrunners = propByPeakDistance_perSeqCells(d_info.running(:,1)==0,:);
this_data_runners = propByPeakDistance_perSeqCells(d_info.running(:,1)==1,:);
%shadedErrorBar(1:size(this_data_nonrunners,2),nanmean(this_data_nonrunners',2)*100,nansem(this_data_nonrunners',2)*100,'lineProps',p.col.nonrunner);
shadedErrorBar(1:size(this_data_runners,2),nanmean(this_data_runners',2)*100,nansem(this_data_runners',2)*100,'lineProps',p.col.runner);
xlim([1,size(this_data_runners,2)+1]-0.5)
xticks([0: 2 :length(distanceBinEdges)]+0.5)
xticklabels({'0','100','200','300','400','500','600'})
ytickformat('percentage')
ylim([0,20])
yticks([0:10:20])
xlabel('Firing field peak (cm)')
ylabel('Proportion of sequence cells')

savefig(F,[save_root_fig,'\Fig2_SequenceAnalysis_PlaceCellProfile.fig']);
saveas(F,[save_root_png,'\Fig2_SequenceAnalysis_PlaceCellProfile.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_SequenceAnalysis_PlaceCellProfile.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_SequenceAnalysis_PlaceCellProfile.txt'],'wt');


%% Fig2_DecodingAnalysis_timeDecodingError_s

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

this_data_nonrunners = timeDecodingError_s(d_info.running(:,1)==0,:);
this_data_runners = timeDecodingError_s(d_info.running(:,1)==1,:);
shadedErrorBar(1:size(this_data_nonrunners,2),nanmean(this_data_nonrunners',2),nansem(this_data_nonrunners',2),'lineProps',p.col.nonrunner);
shadedErrorBar(1:size(this_data_runners,2),nanmean(this_data_runners',2),nansem(this_data_runners',2),'lineProps',p.col.runner);

xlim([0,25]); xticks([0:5:25]); xticklabels({'0','1','2','3','4','5'})
xlabel('Time (s)')
ylim([0,3]); yticks([0:1:3]);
ylabel('Time decoding error (s)')

savefig(F,[save_root_fig,'\Fig2_DecodingAnalysis_timeDecodingError_s.fig']);
saveas(F,[save_root_png,'\Fig2_DecodingAnalysis_timeDecodingError_s.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_DecodingAnalysis_timeDecodingError_s.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_DecodingAnalysis_timeDecodingError_s.txt'],'wt');


%% Fig2_DecodingAnalysis_timeDecodingError_withinCat_s

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

this_data_nonrunners = timeDecodingError_withinCat_s(d_info.running(:,1)==0,:);
this_data_runners = timeDecodingError_withinCat_s(d_info.running(:,1)==1,:);
shadedErrorBar(1:size(this_data_nonrunners,2),nanmean(this_data_nonrunners',2),nansem(this_data_nonrunners',2),'lineProps',p.col.nonrunner);
shadedErrorBar(1:size(this_data_runners,2),nanmean(this_data_runners',2),nansem(this_data_runners',2),'lineProps',p.col.runner);

xlim([0,25]); xticks([0:5:25]); xticklabels({'0','1','2','3','4','5'})
xlabel('Time (s)')
ylim([0,3]); yticks([0:1:3]);
ylabel('Time decoding error [within] (s)')

savefig(F,[save_root_fig,'\Fig2_DecodingAnalysis_timeDecodingError_withinCat_s.fig']);
saveas(F,[save_root_png,'\Fig2_DecodingAnalysis_timeDecodingError_withinCat_s.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_DecodingAnalysis_timeDecodingError_withinCat_s.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_DecodingAnalysis_timeDecodingError_withinCat_s.txt'],'wt');


%% Fig2_DecodingAnalysis_timeDecodingError_s_allCells

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

this_data_nonrunners = timeDecodingError_s_allCells(d_info.running(:,1)==0,:);
this_data_runners = timeDecodingError_s_allCells(d_info.running(:,1)==1,:);
shadedErrorBar(1:size(this_data_nonrunners,2),nanmean(this_data_nonrunners',2),nansem(this_data_nonrunners',2),'lineProps',p.col.nonrunner);
shadedErrorBar(1:size(this_data_runners,2),nanmean(this_data_runners',2),nansem(this_data_runners',2),'lineProps',p.col.runner);

xlim([0,25]); xticks([0:5:25]); xticklabels({'0','1','2','3','4','5'})
xlabel('Time (s)')
ylim([0,3]); yticks([0:1:3]);
ylabel('Time decoding error (s)')

savefig(F,[save_root_fig,'\Fig2_DecodingAnalysis_timeDecodingError_s_allCells.fig']);
saveas(F,[save_root_png,'\Fig2_DecodingAnalysis_timeDecodingError_s_allCells.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_DecodingAnalysis_timeDecodingError_s_allCells.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_DecodingAnalysis_timeDecodingError_s_allCells.txt'],'wt');


%% Fig2_DecodingAnalysis_timeDecodingError_withinCat_s_allCells

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

this_data_nonrunners = timeDecodingError_withinCat_s_allCells(d_info.running(:,1)==0,:);
this_data_runners = timeDecodingError_withinCat_s_allCells(d_info.running(:,1)==1,:);
shadedErrorBar(1:size(this_data_nonrunners,2),nanmean(this_data_nonrunners',2),nansem(this_data_nonrunners',2),'lineProps',p.col.nonrunner);
shadedErrorBar(1:size(this_data_runners,2),nanmean(this_data_runners',2),nansem(this_data_runners',2),'lineProps',p.col.runner);

xlim([0,25]); xticks([0:5:25]); xticklabels({'0','1','2','3','4','5'})
xlabel('Time (s)')
ylim([0,3]); yticks([0:1:3]);
ylabel('Time decoding error [within] (s)')

savefig(F,[save_root_fig,'\Fig2_DecodingAnalysis_timeDecodingError_withinCat_s_allCells.fig']);
saveas(F,[save_root_png,'\Fig2_DecodingAnalysis_timeDecodingError_withinCat_s_allCells.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_DecodingAnalysis_timeDecodingError_withinCat_s_allCells.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_DecodingAnalysis_timeDecodingError_withinCat_s_allCells.txt'],'wt');


%% Fig2_DecodingAnalysis_typeDecodingCorrect_allCells

F = paper_figure([0,0.5,mm2inch(2*34),mm2inch(34)]); hold on;

this_data_nonrunners = typeDecodingCorrect_allCells(d_info.running(:,1)==0,:);
this_data_runners = typeDecodingCorrect_allCells(d_info.running(:,1)==1,:);
shadedErrorBar(1:size(this_data_nonrunners,2),nanmean(this_data_nonrunners',2)*100,nansem(this_data_nonrunners',2)*100,'lineProps',p.col.nonrunner);
shadedErrorBar(1:size(this_data_runners,2),nanmean(this_data_runners',2)*100,nansem(this_data_runners',2)*100,'lineProps',p.col.runner);

xlim([0,25]); xticks([0:5:25]); xticklabels({'0','1','2','3','4','5'});
xlabel('Time (s)')
ytickformat('percentage')
ylim([50,100]); yticks([50,75,100]);
ylabel('1st odor decoding accuracy')

savefig(F,[save_root_fig,'\Fig2_DecodingAnalysis_typeDecodingCorrect_allCells.fig']);
saveas(F,[save_root_png,'\Fig2_DecodingAnalysis_typeDecodingCorrect_allCells.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_DecodingAnalysis_typeDecodingCorrect_allCells.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_DecodingAnalysis_typeDecodingCorrect_allCells.txt'],'wt');


%% --- Performance correlations ---

%% Fig2_PerformanceCorr_timeDecodingError_withinCat_s

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_x_nonrunners = nanmean(timeDecodingError_withinCat_s(d_info.running(:,1)==0,:),2);
this_data_x_runners = nanmean(timeDecodingError_withinCat_s(d_info.running(:,1)==1,:),2);

this_data_y_nonrunners = correct(d_info.running(:,1)==0);
this_data_y_runners = correct(d_info.running(:,1)==1);
this_data_x = [this_data_x_nonrunners;this_data_x_runners];
this_data_y = [this_data_y_nonrunners;this_data_y_runners];

yline(50,'k:');
scatter(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.nonrunner)
scatter(this_data_x_runners(:),this_data_y_runners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k')

% xlim([50,100]); xticks([50,75,100]);
xlabel('Time decoding error (s)')

ytickformat('percentage')
ylim([35,100])
yticks([50,100])
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig2_PerformanceCorr_timeDecodingError_withinCat_s.fig']);
saveas(F,[save_root_png,'\Fig2_PerformanceCorr_timeDecodingError_withinCat_s.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_PerformanceCorr_timeDecodingError_withinCat_s.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_PerformanceCorr_timeDecodingError_withinCat_s_late

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_x_nonrunners = nanmean(timeDecodingError_withinCat_s(d_info.running(:,1)==0,17:24),2);
this_data_x_runners = nanmean(timeDecodingError_withinCat_s(d_info.running(:,1)==1,17:24),2);

this_data_y_nonrunners = correct(d_info.running(:,1)==0);
this_data_y_runners = correct(d_info.running(:,1)==1);
this_data_x = [this_data_x_nonrunners;this_data_x_runners];
this_data_y = [this_data_y_nonrunners;this_data_y_runners];

yline(50,'k:');
scatter(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.nonrunner)
scatter(this_data_x_runners(:),this_data_y_runners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k')

% xlim([50,100]); xticks([50,75,100]);
xlabel('Time decoding error (late; s)')

ytickformat('percentage')
ylim([35,100])
yticks([50,100])
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig2_PerformanceCorr_timeDecodingError_withinCat_s_late.fig']);
saveas(F,[save_root_png,'\Fig2_PerformanceCorr_timeDecodingError_withinCat_s_late.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_PerformanceCorr_timeDecodingError_withinCat_s_late.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_PerformanceCorr_typeDecodingCorrect

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_x_nonrunners = nanmean(typeDecodingCorrect(d_info.running(:,1)==0,:),2)*100;
this_data_x_runners = nanmean(typeDecodingCorrect(d_info.running(:,1)==1,:),2)*100;

this_data_y_nonrunners = correct(d_info.running(:,1)==0);
this_data_y_runners = correct(d_info.running(:,1)==1);
this_data_x = [this_data_x_nonrunners;this_data_x_runners];
this_data_y = [this_data_y_nonrunners;this_data_y_runners];

yline(50,'k:');
scatter(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.nonrunner)
scatter(this_data_x_runners(:),this_data_y_runners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k')

xlim([50,100]); xticks([50,75,100]);
xlabel('1st odor decoding accuracy')

ytickformat('percentage')
ylim([35,100])
yticks([50,100])
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig2_PerformanceCorr_typeDecodingCorrect.fig']);
saveas(F,[save_root_png,'\Fig2_PerformanceCorr_typeDecodingCorrect.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_PerformanceCorr_typeDecodingCorrect.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_PerformanceCorr_typeDecodingCorrect_late

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_x_nonrunners = nanmean(typeDecodingCorrect(d_info.running(:,1)==0,17:24),2)*100;
this_data_x_runners = nanmean(typeDecodingCorrect(d_info.running(:,1)==1,17:24),2)*100;

this_data_y_nonrunners = correct(d_info.running(:,1)==0);
this_data_y_runners = correct(d_info.running(:,1)==1);
this_data_x = [this_data_x_nonrunners;this_data_x_runners];
this_data_y = [this_data_y_nonrunners;this_data_y_runners];

yline(50,'k:');
scatter(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.nonrunner)
scatter(this_data_x_runners(:),this_data_y_runners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k')

xlim([50,100]); xticks([50,75,100]);
xlabel('1st odor decoding accuracy (late)')

ytickformat('percentage')
ylim([35,100])
yticks([50,100])
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig2_PerformanceCorr_typeDecodingCorrect_late.fig']);
saveas(F,[save_root_png,'\Fig2_PerformanceCorr_typeDecodingCorrect_late.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_PerformanceCorr_typeDecodingCorrect_late.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_PerformanceCorr_SequenceCells_late

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_x_nonrunners = propByCellType_perSeqCells.pref(d_info.running(:,1)==0,3)*100;
this_data_x_runners = propByCellType_perSeqCells.pref(d_info.running(:,1)==1,3)*100;

this_data_y_nonrunners = correct(d_info.running(:,1)==0);
this_data_y_runners = correct(d_info.running(:,1)==1);
this_data_x = [this_data_x_nonrunners;this_data_x_runners];
this_data_y = [this_data_y_nonrunners;this_data_y_runners];

yline(50,'k:');
scatter(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.nonrunner)
scatter(this_data_x_runners(:),this_data_y_runners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k')

xlim([0,100])
xticks([0:50:100])
xlabel({'Proportion of sequence cells (late)'});

ytickformat('percentage')
ylim([35,100])
yticks([50,100])
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig2_PerformanceCorr_SequenceCells_late.fig']);
saveas(F,[save_root_png,'\Fig2_PerformanceCorr_SequenceCells_late.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_PerformanceCorr_SequenceCells_late.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_PerformanceCorr_Selectivity

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_x_nonrunners = selectivity.pref(d_info.running(:,1)==0);
this_data_x_runners = selectivity.pref(d_info.running(:,1)==1);

this_data_y_nonrunners = correct(d_info.running(:,1)==0);
this_data_y_runners = correct(d_info.running(:,1)==1);
this_data_x = [this_data_x_nonrunners;this_data_x_runners];
this_data_y = [this_data_y_nonrunners;this_data_y_runners];

yline(50,'k:');
scatter(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.nonrunner)
scatter(this_data_x_runners(:),this_data_y_runners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k')

xlim([0,0.3])
xticks([0:0.1:0.3])
xlabel({'Selectivity'})
ytickformat('percentage')
ylim([35,100])
yticks([50,100])
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig2_PerformanceCorr_Selectivity.fig']);
saveas(F,[save_root_png,'\Fig2_PerformanceCorr_Selectivity.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_PerformanceCorr_Selectivity.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_PerformanceCorr_ExpB

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_x_nonrunners = exp_b(d_info.running(:,1)==0);
this_data_x_runners = exp_b(d_info.running(:,1)==1);

this_data_y_nonrunners = correct(d_info.running(:,1)==0);
this_data_y_runners = correct(d_info.running(:,1)==1);
this_data_x = [this_data_x_nonrunners;this_data_x_runners];
this_data_y = [this_data_y_nonrunners;this_data_y_runners];

yline(50,'k:');
scatter(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.nonrunner)
scatter(this_data_x_runners(:),this_data_y_runners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k')

xlim([-1.5,0])
xticks([-1.5:0.5:0])
xlabel({'Sequence shape','(decay rate, s-1)'});
ytickformat('percentage')
ylim([35,100])
yticks([50,100])
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig2_PerformanceCorr_ExpB.fig']);
saveas(F,[save_root_png,'\Fig2_PerformanceCorr_ExpB.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_PerformanceCorr_ExpB.pdf']); set(gcf,'Color',[1,1,1])

