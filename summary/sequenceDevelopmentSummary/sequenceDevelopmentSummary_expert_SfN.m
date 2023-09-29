%% SfN

% d_info = selectDataset(d_info,'-g02345678-d1-e01-r01-p01-l01-i00',sheet,path,ops);


%% Preparations

binEdges = [0,1.5,3,5.3];
numBins = length(binEdges)-1;
binLabels = {'Early peak (0 s - 1.5 s)','Middle peak (1.5 s - 3 s)','Late peak (3 s - 5.3 s)'};


%% Get data

correct = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~isempty(d{i,1})
        correct(i) = d{i,1}.perf.general.correct;
    end
end

warp_exp1_a = nan(d_info.numAnimals,1);
warp_exp1_b = nan(d_info.numAnimals,1);
numSeqCells = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if isfield(d{i,1},'warp_all')
        warp_exp1_a(i) = d{i,1}.warp_all.exp1.a;
        warp_exp1_b(i) = d{i,1}.warp_all.exp1.b;
        try
            numSeqCells(i) = d{i,1}.warp_all.input.numSequenceCells;
        catch
        end
    end
end


%% Get sequence cell by time

tng_struct = 'tngnn_all';

binEdges = [0,1.5,3,5.3];
numBins = length(binEdges)-1;
binLabels = {'Early peak (0 s - 1.5 s)','Middle peak (1.5 s - 3 s)','Late peak (3 s - 5.3 s)'};

seqCellsByTime_num = nan(d_info.numAnimals,numBins);
seqCellsByTime_prop = nan(d_info.numAnimals,numBins);
for i=1:d_info.numAnimals
    if isfield(d{i,1},tng_struct)
        
        % get idcs
        this_numCells = length(find(d{i,1}.(tng_struct).prop.iscell==1));
        these_idcs_A = find(d{i,1}.(tng_struct).passed.AW.Aonly==1);
        these_idcs_X = find(d{i,1}.(tng_struct).passed.AW.Xonly==1);

        % get peak locations
        these_peakTimes_A = d{i,1}.(tng_struct).firingField.A_AW.peakLocation_s(these_idcs_A);
        these_peakTimes_X = d{i,1}.(tng_struct).firingField.X_AW.peakLocation_s(these_idcs_X);
        these_peakTimes = [these_peakTimes_A; these_peakTimes_X];

        % bin peak locations
        temp = discretize(these_peakTimes,binEdges);
        for n=1:numBins
            seqCellsByTime_num(i,n) = nansum(temp==n);
            seqCellsByTime_prop(i,n) = nansum(temp==n) / this_numCells;
        end
    end
end


%% Sequence characterisation runners vs non-runners

nrows = 1; ncols = 5; m=0;
F = default_figure([0,3.5,18,2.5]);
these_labels = {'non-runners','runners'};

% get data
this_data_a = nan(d_info.numAnimals,2);
temp = warp_exp1_a(~d_info.excl(:,1) & d_info.running(:,1)==0);
this_data_a(1:length(temp),1) = temp;
temp = warp_exp1_a(~d_info.excl(:,1) & d_info.running(:,1)==1);
this_data_a(1:length(temp),2) = temp;
this_data_b = nan(d_info.numAnimals,2);
temp = warp_exp1_b(~d_info.excl(:,1) & d_info.running(:,1)==0);
this_data_b(1:length(temp),1) = temp;
temp = warp_exp1_b(~d_info.excl(:,1) & d_info.running(:,1)==1);
this_data_b(1:length(temp),2) = temp;

% a) coefficient a
m=m+1; subplot(nrows,ncols,m); hold on;
v = violinplot(this_data_a,these_labels);
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k';
xlim([0.6,2.4])
ylim([0,0.3])
yticks([0,0.3])
ylabel(['Coefficient a',newline,'(initial proportion)']);
title(['p=',num2str(ranksum(this_data_a(:,1),this_data_a(:,2)),2),' (ranksum)'],'FontSize',10,'FontWeight','normal');

% b) coefficient b
m=m+1; subplot(nrows,ncols,m); hold on;
v = violinplot(this_data_b,these_labels);
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k';
xlim([0.6,2.4])
ylim([-1.5,0])
yticks([-1.5,0])
ylabel(['Coefficient b',newline,'(decay rate, s-1)']);
title(['p=',num2str(ranksum(this_data_b(:,1),this_data_b(:,2)),2),' (ranksum)'],'FontSize',10,'FontWeight','normal');

% get data
this_data_early = nan(d_info.numAnimals,2);
temp = seqCellsByTime_prop(~d_info.excl(:,1) & d_info.running(:,1)==0, 1);
this_data_early(1:length(temp),1) = temp;
temp = seqCellsByTime_prop(~d_info.excl(:,1) & d_info.running(:,1)==1, 1);
this_data_early(1:length(temp),2) = temp;
this_data_middle = nan(d_info.numAnimals,2);
temp = seqCellsByTime_prop(~d_info.excl(:,1) & d_info.running(:,1)==0, 2);
this_data_middle(1:length(temp),1) = temp;
temp = seqCellsByTime_prop(~d_info.excl(:,1) & d_info.running(:,1)==1, 2);
this_data_middle(1:length(temp),2) = temp;
this_data_late = nan(d_info.numAnimals,2);
temp = seqCellsByTime_prop(~d_info.excl(:,1) & d_info.running(:,1)==0, 3);
this_data_late(1:length(temp),1) = temp;
temp = seqCellsByTime_prop(~d_info.excl(:,1) & d_info.running(:,1)==1, 3);
this_data_late(1:length(temp),2) = temp;

% a) early
m=m+1; subplot(nrows,ncols,m); hold on;
v = violinplot(this_data_early,these_labels);
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k';
xlim([0.6,2.4])
ylim([0,0.5])
yticks([0,0.5])
ylabel('Proportion of cells');
NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(NE(1),NE(2),['p=',num2str(ranksum(this_data_early(:,1),this_data_early(:,2)),2),' (ranksum)'],'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
title(binLabels{1});

% b) middle
m=m+1; subplot(nrows,ncols,m); hold on;
v = violinplot(this_data_middle,these_labels);
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k';
xlim([0.6,2.4])
ylim([0,0.2])
yticks([0,0.2])
ylabel('Proportion of cells');
NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(NE(1),NE(2),...
    ['p=',num2str(ranksum(this_data_middle(:,1),this_data_middle(:,2)),2),' (ranksum)'],'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
title(binLabels{2});

% c) late
m=m+1; subplot(nrows,ncols,m); hold on;
v = violinplot(this_data_late,these_labels);
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k';
xlim([0.6,2.4])
ylim([0,0.2])
yticks([0,0.2])
ylabel('Proportion of cells');
NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(NE(1),NE(2),...
    ['p=',num2str(ranksum(this_data_late(:,1),this_data_late(:,2)),2),' (ranksum)'],'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
title(binLabels{3});

savefig(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Expert Imaging\new\expert.fig']);
saveas(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Expert Imaging\new\expert.png']);





%% Number of sequence cells

nrows = 1; ncols = 1; m=0;
F = default_figure([0,3.5,5,3.5]);
these_labels = {'non-runners','runners'};

% get data
this_data_a = nan(d_info.numAnimals,2);
temp = numSeqCells(~d_info.excl(:,1) & d_info.running(:,1)==0);
this_data_a(1:length(temp),1) = temp;
temp = numSeqCells(~d_info.excl(:,1) & d_info.running(:,1)==1);
this_data_a(1:length(temp),2) = temp;

% a) coefficient a
m=m+1; subplot(nrows,ncols,m); hold on;
v = violinplot(this_data_a,these_labels);
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k';
xlim([0.6,2.4])
ylim([0,500])
yticks([0:100:500])
ylabel(['Number of sequence cells']);
title(['p=',num2str(ranksum(this_data_a(:,1),this_data_a(:,2)),2),' (ranksum)'],'FontSize',10,'FontWeight','normal');

savefig(F,'C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Supp\numSeqCells_expert.fig');
saveas(F,'C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Supp\numSeqCells_expert.png');



