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
        numSeqCells(i) = length(d{i,1}.warp_all.input.idcs_Aonly) + length(d{i,1}.warp_all.input.idcs_Xonly);
    end
end


%% Get sequence cell by time

binEdges = [0,1.5,3,5.3];
numBins = length(binEdges)-1;
binLabels = {'Early peak (0 s - 1.5 s)','Middle peak (1.5 s - 3 s)','Late peak (3 s - 5.3 s)'};

seqCellsByTime_num = nan(d_info.numAnimals,numBins);
seqCellsByTime_prop = nan(d_info.numAnimals,numBins);
for i=1:d_info.numAnimals
    if isfield(d{i,1},'tng_all')
        
        % get idcs
        this_numCells = length(find(d{i,1}.tng_all.prop.iscell==1));
        these_idcs_A = find(d{i,1}.tng_all.passed.AW.Aonly==1);
        these_idcs_X = find(d{i,1}.tng_all.passed.AW.Xonly==1);

        % get peak locations
        these_peakTimes_A = d{i,1}.tng_all.firingField.A_AW.peakLocation_s(these_idcs_A);
        these_peakTimes_X = d{i,1}.tng_all.firingField.X_AW.peakLocation_s(these_idcs_X);
        these_peakTimes = [these_peakTimes_A; these_peakTimes_X];

        % bin peak locations
        temp = discretize(these_peakTimes,binEdges);
        for n=1:numBins
            seqCellsByTime_num(i,n) = nansum(temp==n);
            seqCellsByTime_prop(i,n) = nansum(temp==n) / this_numCells;
        end
    end
end


%%

% select sessions
these_animals = 'd1-r0'; these_sessions = ~d_info.excl(:,1)  & d_info.running(:,1)==0; %  & d_info.running(:,1)==1

F = default_figure();

% select data
this_data_x = warp_exp1_a(these_sessions);
this_data_y = correct(these_sessions);

subplot(1,2,1); hold on;
yline(50,'k:');
plot(this_data_x(:),this_data_y(:)*100,'ko')
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k');
%xlim([0,0.25])
ylim([30,100])
xlabel('Warping (a)')
ylabel('Performance (%correct)')
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['\rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
title('a')

% select data
this_data_x = warp_exp1_b(these_sessions);
this_data_y = correct(these_sessions);

subplot(1,2,2); hold on;
yline(50,'k:');
plot(this_data_x(:),this_data_y(:)*100,'ko')
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k');
%xlim([-2,1])
ylim([30,100])
xlabel('Warping (b)')
ylabel('Performance (%correct)')
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['\rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
title('b')

suptitle(['Warping-performance correlation (',these_animals,')'])
savefig(F,[path.root_summary,'plots\warping\warping_performance_',these_animals,'.fig']);
saveas(F,[path.root_summary,'plots\warping\warping_performance_',these_animals,'.png']);








%% --- SFN FIGURES ---

% runners vs non-runners
% a) warp
% b) early, middle, late
% c) performance
% d) corr(warp,performance)
% e) corr(early, middle, late ,performance)


%% Sequence characterisation - performance - runners vs non-runners

nrows = 1; ncols = 1; m=0;
F = default_figure([-20,0.5,5,5]);
these_labels = {'non-runners','runners'};

% get data
this_data_correct = nan(d_info.numAnimals,2);
temp = correct(~d_info.excl(:,1) & d_info.running(:,1)==0);
this_data_correct(1:length(temp),1) = temp;
temp = correct(~d_info.excl(:,1) & d_info.running(:,1)==1);
this_data_correct(1:length(temp),2) = temp;

% performance
m=m+1; subplot(nrows,ncols,m); hold on;
v = violinplot(this_data_correct*100,these_labels);
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k';
ylim([35,100])
ylabel('Performance (%correct)');
title(['p=',num2str(ranksum(this_data_correct(:,1),this_data_correct(:,2)),2),' (ranksum)'],'FontSize',10,'FontWeight','normal');

savefig(F,[path.root_summary,'plots\warping\performance_d1_r01.fig']);
saveas(F,[path.root_summary,'plots\warping\performance_d1_r01.png']);


%% Sequence characterisation - exponential fit - runners vs non-runners

nrows = 1; ncols = 2; m=0;
F = default_figure([-20,0.5,10,5]);
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
ylim([0,0.25])
ylabel('Exponential coefficient a (initial proportion)');
title(['p=',num2str(ranksum(this_data_a(:,1),this_data_a(:,2)),2),' (ranksum)'],'FontSize',10,'FontWeight','normal');

% b) coefficient b
m=m+1; subplot(nrows,ncols,m); hold on;
v = violinplot(this_data_b,these_labels);
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k';
ylim([-1.2,0])
ylabel('Exponential coefficient b (decay rate)');
title(['p=',num2str(ranksum(this_data_b(:,1),this_data_b(:,2)),2),' (ranksum)'],'FontSize',10,'FontWeight','normal');

savefig(F,[path.root_summary,'plots\warping\warping_exponential_d1_r01.fig']);
saveas(F,[path.root_summary,'plots\warping\warping_exponential_d1_r01.png']);


%% Sequence characterisation - sequence cells by time - runners vs non-runners

nrows = 1; ncols = 3; m=0;
F = default_figure([-20,0.5,10,5]);
these_labels = {'non-runners','runners'};

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
ylim([0,0.5])
ylabel('Proportion of cells');
NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(NE(1),NE(2),['p=',num2str(ranksum(this_data_early(:,1),this_data_early(:,2)),2),' (ranksum)'],'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
title(binLabels{1});

% b) middle
m=m+1; subplot(nrows,ncols,m); hold on;
v = violinplot(this_data_middle,these_labels);
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k';
ylim([0,0.5])
ylabel('Proportion of cells');
NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(NE(1),NE(2),...
    ['p=',num2str(ranksum(this_data_middle(:,1),this_data_middle(:,2)),2),' (ranksum)'],'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
title(binLabels{2});

% c) late
m=m+1; subplot(nrows,ncols,m); hold on;
v = violinplot(this_data_late,these_labels);
v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k';
ylim([0,0.5])
ylabel('Proportion of cells');
NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(NE(1),NE(2),...
    ['p=',num2str(ranksum(this_data_late(:,1),this_data_late(:,2)),2),' (ranksum)'],'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
title(binLabels{3});

savefig(F,[path.root_summary,'plots\warping\warping_sequence_d1_r01.fig']);
saveas(F,[path.root_summary,'plots\warping\warping_sequence_d1_r01.png']);






