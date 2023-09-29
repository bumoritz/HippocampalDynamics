%% Get data

correct = nan(d_info.numAnimals,d_info.numDays);
correct_switch = nan(d_info.numAnimals,2);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if ~isempty(d{i,j})
            correct(i,j) = d{i,j}.perf.general.correct;
        end
    end
    correct_switch(i,1) = nanmean([correct(i,2),correct(i,4)]);
    correct_switch(i,2) = nanmean([correct(i,3),correct(i,5)]);
end

blocks = [5,8,5,8,5];
correct_best100t = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if ~isempty(d{i,j})
            temp = [];
            for k=1:blocks(j)
                these_blocks_20t = (k-1)*5+1:k*5;
                temp = [temp, nanmean(d{i,j}.perf.blocks_general.correct(these_blocks_20t))];
            end
            correct_best100t(i,j) = nanmax(temp);
        end
    end
end

blocks = [5,8,5,8,5]*5;
correct_best20t = nan(d_info.numAnimals,d_info.numDays);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if ~isempty(d{i,j})
            temp = [];
            for k=1:blocks(j)
                temp = [temp, nanmean(d{i,j}.perf.blocks_general.correct(k))];
            end
            correct_best20t(i,j) = nanmax(temp);
        end
    end
end

% speed = nan(d_info.numAnimals,d_info.numDays);
% speed_switch = nan(d_info.numAnimals,2);
% for i=1:d_info.numAnimals
%     for j=1:d_info.numDays
%         if ~isempty(d{i,j})
%             speed(i,j) = d{i,j}.perf.general.correct;
%         end
%     end
%     speed_switch(i,1) = nanmean([speed(i,2),speed(i,4)]);
%     speed_switch(i,2) = nanmean([speed(i,3),speed(i,5)]);
% end

warp_exp1_a = nan(d_info.numAnimals,d_info.numDays);
warp_exp1_b = nan(d_info.numAnimals,d_info.numDays);
warp_exp1_a_switch = nan(d_info.numAnimals,2);
warp_exp1_b_switch = nan(d_info.numAnimals,2);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if isfield(d{i,j},'warp_all')
            warp_exp1_a(i,j) = d{i,j}.warp_all.exp1.a;
            warp_exp1_b(i,j) = d{i,j}.warp_all.exp1.b;
        end
    end
    warp_exp1_a_switch(i,1) = nanmean([warp_exp1_a(i,2),warp_exp1_a(i,4)]);
    warp_exp1_a_switch(i,2) = nanmean([warp_exp1_a(i,3),warp_exp1_a(i,5)]);
    warp_exp1_b_switch(i,1) = nanmean([warp_exp1_b(i,2),warp_exp1_b(i,4)]);
    warp_exp1_b_switch(i,2) = nanmean([warp_exp1_b(i,3),warp_exp1_b(i,5)]);
end

nonrunners = ~d_info.excl & d_info.running==0;
runners = ~d_info.excl & d_info.running==1;
nonrunners_switch_d1 = false(d_info.numAnimals,1);
runners_switch_d1 = false(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~d_info.excl(i,2) && ~d_info.excl(i,4)
        runners_switch_d1(i) = runners(i,2) && runners(i,4);
        nonrunners_switch_d1(i) = nonrunners(i,2) && nonrunners(i,4);
    elseif ~d_info.excl(i,2)
        runners_switch_d1(i) = runners(i,2);
        nonrunners_switch_d1(i) = nonrunners(i,2);
    elseif ~d_info.excl(i,4)
        runners_switch_d1(i) = runners(i,4);
        nonrunners_switch_d1(i) = nonrunners(i,4);
    end
end
nonrunners_switch_d2 = false(d_info.numAnimals,1);
runners_switch_d2 = false(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~d_info.excl(i,3) && ~d_info.excl(i,5)
        runners_switch_d2(i) = runners(i,3) && runners(i,5);
        nonrunners_switch_d2(i) = nonrunners(i,3) && nonrunners(i,5);
    elseif ~d_info.excl(i,3)
        runners_switch_d2(i) = runners(i,3);
        nonrunners_switch_d2(i) = nonrunners(i,3);
    elseif ~d_info.excl(i,5)
        runners_switch_d2(i) = runners(i,5);
        nonrunners_switch_d2(i) = nonrunners(i,5);
    end
end
nonrunners_switch_d12 = false(d_info.numAnimals,1);
runners_switch_d12 = false(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if (~d_info.excl(i,2) || ~d_info.excl(i,4)) && (~d_info.excl(i,3) || ~d_info.excl(i,5))
        runners_switch_d12(i) = runners_switch_d1(i) && runners_switch_d2(i);
        nonrunners_switch_d12(i) = nonrunners_switch_d1(i) && nonrunners_switch_d2(i);       
    elseif (~d_info.excl(i,2) || ~d_info.excl(i,4))
        runners_switch_d12(i) = runners_switch_d1(i);
        nonrunners_switch_d12(i) = nonrunners_switch_d1(i);
    elseif (~d_info.excl(i,3) || ~d_info.excl(i,5))
        runners_switch_d12(i) = runners_switch_d2(i);
        nonrunners_switch_d12(i) = nonrunners_switch_d2(i);
    end
end


%% Get sequence cell by time

binEdges = [0,1.5,3,5.3];
numBins = length(binEdges)-1;
binLabels = {'Early peak (0 s - 1.5 s)','Middle peak (1.5 s - 3 s)','Late peak (3 s - 5.3 s)'};

seqCellsByTime_num = nan(d_info.numAnimals,d_info.numDays,numBins);
seqCellsByTime_prop = nan(d_info.numAnimals,d_info.numDays,numBins);
seqCellsByTime_num_switch = nan(d_info.numAnimals,2,numBins);
seqCellsByTime_prop_switch = nan(d_info.numAnimals,2,numBins);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if isfield(d{i,j},'tng_all')

            % get idcs
            this_numCells = length(find(d{i,j}.tng_all.prop.iscell==1));
            these_idcs_A = find(d{i,j}.tng_all.passed.AW.Aonly==1);
            these_idcs_X = find(d{i,j}.tng_all.passed.AW.Xonly==1);

            % get peak locations
            these_peakTimes_A = d{i,j}.tng_all.firingField.A_AW.peakLocation_s(these_idcs_A);
            these_peakTimes_X = d{i,j}.tng_all.firingField.X_AW.peakLocation_s(these_idcs_X);
            these_peakTimes = [these_peakTimes_A; these_peakTimes_X];

            % bin peak locations
            temp = discretize(these_peakTimes,binEdges);
            for n=1:numBins
                seqCellsByTime_num(i,j,n) = nansum(temp==n);
                seqCellsByTime_prop(i,j,n) = nansum(temp==n) / this_numCells;
            end
        end
    end
    for n=1:numBins
         seqCellsByTime_num_switch(i,1,n) = nanmean([seqCellsByTime_num(i,2,n),seqCellsByTime_num(i,4,n)]);
         seqCellsByTime_num_switch(i,2,n) = nanmean([seqCellsByTime_num(i,3,n),seqCellsByTime_num(i,5,n)]);
         seqCellsByTime_prop_switch(i,1,n) = nanmean([seqCellsByTime_prop(i,2,n),seqCellsByTime_prop(i,4,n)]);
         seqCellsByTime_prop_switch(i,2,n) = nanmean([seqCellsByTime_prop(i,3,n),seqCellsByTime_prop(i,5,n)]);
    end
end


%% Sequence characterisation - performance - runners vs non-runners - bar

nrows = 1; ncols = 5; m=0;
F = default_figure([-20,0.5,20,2.5]);
these_labels = categorical({'non-runners','runners'});
these_labels = reordercats(these_labels,{'non-runners','runners'});

for j=1:d_info.numDays
    % get data
    this_data_correct = nan(d_info.numAnimals,2);
    temp = correct(~d_info.excl(:,j) & d_info.running(:,j)==0, j);
    this_data_correct(1:length(temp),1) = temp;
    temp = correct(~d_info.excl(:,j) & d_info.running(:,j)==1, j);
    this_data_correct(1:length(temp),2) = temp;

    % performance
    m=m+1; subplot(nrows,ncols,m); hold on;
    v = bar(these_labels,nanmean(this_data_correct,1)*100);
    yline(50,'k:')
    for k=1:length(this_data_correct)
        plot(these_labels,this_data_correct*100,'ok')
    end
    v.FaceColor = 'flat'; v.CData(1,:) = p.col.nonrunner; v.CData(2,:) = p.col.runner;
    ylim([35,100])
    ytickformat('percentage')
    ylabel('Performance');
    text(2.5,98,['p=',num2str(ranksum(this_data_correct(:,1),this_data_correct(:,2)),2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title(['Day ',num2str(j)])
end

savefig(F,[path.root_summary,'plots\warping\performance_d2345_r01.fig']);
saveas(F,[path.root_summary,'plots\warping\performance_d2345_r01.png']);


%% DABEST


%% Sequence characterisation - performance - runners vs non-runners - bar

nrows = 1; ncols = 5; m=0;
F = default_figure([-20,0.5,20,5]);
these_labels = categorical({'non-runners','runners'});
these_labels = reordercats(these_labels,{'non-runners','runners'});

for j=1:d_info.numDays
    % get data
    this_data_correct = nan(d_info.numAnimals,2);
    temp = correct_best100t(~d_info.excl(:,j) & d_info.running(:,j)==0, j);
    this_data_correct(1:length(temp),1) = temp;
    temp = correct_best100t(~d_info.excl(:,j) & d_info.running(:,j)==1, j);
    this_data_correct(1:length(temp),2) = temp;

    % performance
    m=m+1; subplot(nrows,ncols,m); hold on;
    v = bar(these_labels,nanmean(this_data_correct,1)*100);
    yline(50,'k:')
    for k=1:length(this_data_correct)
        plot(these_labels,this_data_correct*100,'ok')
    end
    v.FaceColor = 'flat'; v.CData(1,:) = p.col.nonrunner; v.CData(2,:) = p.col.runner;
    ylim([35,100])
    ytickformat('percentage')
    ylabel('Performance (best block of 100 trials)');
    text(2.5,98,['p=',num2str(ranksum(this_data_correct(:,1),this_data_correct(:,2)),2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title(['Day ',num2str(j)])
end

savefig(F,[path.root_summary,'plots\warping\performance_best100t_d2345_r01.fig']);
saveas(F,[path.root_summary,'plots\warping\performance_best100t_d2345_r01.png']);


%% Sequence characterisation - performance - runners vs non-runners - bar

nrows = 1; ncols = 5; m=0;
F = default_figure([-20,0.5,20,5]);
these_labels = categorical({'non-runners','runners'});
these_labels = reordercats(these_labels,{'non-runners','runners'});

for j=1:d_info.numDays
    % get data
    this_data_correct = nan(d_info.numAnimals,2);
    temp = correct_best20t(~d_info.excl(:,j) & d_info.running(:,j)==0, j);
    this_data_correct(1:length(temp),1) = temp;
    temp = correct_best20t(~d_info.excl(:,j) & d_info.running(:,j)==1, j);
    this_data_correct(1:length(temp),2) = temp;

    % performance
    m=m+1; subplot(nrows,ncols,m); hold on;
    v = bar(these_labels,nanmean(this_data_correct,1)*100);
    yline(50,'k:')
    for k=1:length(this_data_correct)
        plot(these_labels,this_data_correct*100,'ok')
    end
    v.FaceColor = 'flat'; v.CData(1,:) = p.col.nonrunner; v.CData(2,:) = p.col.runner;
    ylim([35,100])
    ytickformat('percentage')
    ylabel('Performance (best block of 20 trials)');
    text(2.5,98,['p=',num2str(ranksum(this_data_correct(:,1),this_data_correct(:,2)),2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title(['Day ',num2str(j)])
end

savefig(F,[path.root_summary,'plots\warping\performance_best20t_d2345_r01.fig']);
saveas(F,[path.root_summary,'plots\warping\performance_best20t_d2345_r01.png']);


%% Sequence characterisation - sequence cells by time - runners vs non-runners - bar

nrows = 1; ncols = 3; m=0;
F = default_figure([-20,0.5,10,5]);
these_labels = {'non-runners','runners'};

% get data %%%%% right now it's set for comparing on the first switch days
this_data_early = nan(d_info.numAnimals,2);
temp = seqCellsByTime_prop_switch(nonrunners_switch_d1, 1, 1);
this_data_early(1:length(temp),1) = temp;
temp = seqCellsByTime_prop_switch(runners_switch_d1, 1, 1);
this_data_early(1:length(temp),2) = temp;

this_data_middle = nan(d_info.numAnimals,2);
temp = seqCellsByTime_prop_switch(nonrunners_switch_d1, 1, 2);
this_data_middle(1:length(temp),1) = temp;
temp = seqCellsByTime_prop_switch(runners_switch_d1, 1, 2);
this_data_middle(1:length(temp),2) = temp;

this_data_late = nan(d_info.numAnimals,2);
temp = seqCellsByTime_prop_switch(nonrunners_switch_d1, 1, 3);
this_data_late(1:length(temp),1) = temp;
temp = seqCellsByTime_prop_switch(runners_switch_d1, 1, 3);
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

savefig(F,[path.root_summary,'plots\warping\warping_sequence_d2345_r01.fig']);
saveas(F,[path.root_summary,'plots\warping\warping_sequence_d2345_r01.png']);



%% Sequence characterisation - sequence cells by time - runners vs non-runners - bar

nrows = 1; ncols = 3; m=0;
F = default_figure([-20,0.5,10,5]);
these_labels = {'non-runners','runners'};

% get data
% 1 - nonrunner, d1
% 2 - runner, d1
% 3 - nonrunner, d2
% 4 - runner, d2

this_data_early = nan(d_info.numAnimals,4);
temp = seqCellsByTime_prop_switch(nonrunners_switch_d1, 1, 1);
this_data_early(1:length(temp),1) = temp;
temp = seqCellsByTime_prop_switch(runners_switch_d1, 1, 1);
this_data_early(1:length(temp),2) = temp;
temp = seqCellsByTime_prop_switch(nonrunners_switch_d2, 2, 1);
this_data_early(1:length(temp),3) = temp;
temp = seqCellsByTime_prop_switch(runners_switch_d2, 2, 1);
this_data_early(1:length(temp),4) = temp;

this_data_middle = nan(d_info.numAnimals,4);
temp = seqCellsByTime_prop_switch(nonrunners_switch_d1, 1, 2);
this_data_middle(1:length(temp),1) = temp;
temp = seqCellsByTime_prop_switch(runners_switch_d1, 1, 2);
this_data_middle(1:length(temp),2) = temp;
temp = seqCellsByTime_prop_switch(nonrunners_switch_d2, 2, 2);
this_data_middle(1:length(temp),3) = temp;
temp = seqCellsByTime_prop_switch(runners_switch_d2, 2, 2);
this_data_middle(1:length(temp),4) = temp;

this_data_late = nan(d_info.numAnimals,4);
temp = seqCellsByTime_prop_switch(nonrunners_switch_d1, 1, 3);
this_data_late(1:length(temp),1) = temp;
temp = seqCellsByTime_prop_switch(runners_switch_d1, 1, 3);
this_data_late(1:length(temp),2) = temp;
temp = seqCellsByTime_prop_switch(nonrunners_switch_d2, 2, 3);
this_data_late(1:length(temp),3) = temp;
temp = seqCellsByTime_prop_switch(runners_switch_d2, 2, 3);
this_data_late(1:length(temp),4) = temp;


%%

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

savefig(F,[path.root_summary,'plots\warping\warping_sequence_d2345_r01.fig']);
saveas(F,[path.root_summary,'plots\warping\warping_sequence_d2345_r01.png']);






%% --- DISCARDED ---



%% Sequence characterisation - performance - runners vs non-runners

% nrows = 1; ncols = 2; m=0;
% F = default_figure([-20,0.5,10,5]);
% these_labels = {'non-runners','runners'};
% 
% for j=1:2
%     % get data
%     this_data_correct = nan(d_info.numAnimals,2);
%     temp = correct_switch((~d_info.excl(:,j+1) & d_info.running(:,j+1)==0) & (~d_info.excl(:,j+3) & d_info.running(:,j+3)==0), j);
%     this_data_correct(1:length(temp),1) = temp;
%     temp = correct_switch((~d_info.excl(:,j+1) & d_info.running(:,j+1)==1) | (~d_info.excl(:,j+3) & d_info.running(:,j+3)==1), j); % |
%     this_data_correct(1:length(temp),2) = temp;
% 
%     % performance
%     m=m+1; subplot(nrows,ncols,m); hold on;
%     v = violinplot(this_data_correct*100,these_labels);
%     v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k';
%     ylim([35,100])
%     ytickformat('percentage')
%     ylabel('Performance (%correct)');
%     SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
%     text(SE(1),SE(2),['p=',num2str(ranksum(this_data_correct(:,1),this_data_correct(:,2)),2),' (ranksum)'],...
%     'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
%     title(['Day ',num2str(j)])
% end
% 
% savefig(F,[path.root_summary,'plots\warping\performance_switches_r01.fig']);
% saveas(F,[path.root_summary,'plots\warping\performance_switches_r01.png']);



%% --- SFN FIGURES ---

% runners vs non-runners
% a) warp
% b) early, middle, late
% c) performance
% d) corr(warp,performance)
% e) corr(early, middle, late ,performance)


%% Sequence characterisation - performance - runners vs non-runners - violin

% nrows = 1; ncols = 5; m=0;
% F = default_figure([-20,0.5,20,5]);
% these_labels = {'non-runners','runners'};
% 
% for j=1:d_info.numDays
%     % get data
%     this_data_correct = nan(d_info.numAnimals,2);
%     temp = correct(~d_info.excl(:,j) & d_info.running(:,j)==0, j);
%     this_data_correct(1:length(temp),1) = temp;
%     temp = correct(~d_info.excl(:,j) & d_info.running(:,j)==1, j);
%     this_data_correct(1:length(temp),2) = temp;
% 
%     % performance
%     m=m+1; subplot(nrows,ncols,m); hold on;
%     v = violinplot(this_data_correct*100,these_labels);
%     v(1).ViolinColor = p.col.nonrunner; v(2).ViolinColor = p.col.runner; v(1).BoxColor = 'k'; v(2).BoxColor = 'k';
%     ylim([40,100])
%     ytickformat('percentage')
%     ylabel('Performance (%correct)');
%     SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
%     text(SE(1),SE(2),['p=',num2str(ranksum(this_data_correct(:,1),this_data_correct(:,2)),2),' (ranksum)'],...
%     'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
%     title(['Day ',num2str(j)])
% end
% 
% savefig(F,[path.root_summary,'plots\warping\performance_d2345_r01_violin.fig']);
% saveas(F,[path.root_summary,'plots\warping\performance_d2345_r01_violin.png']);


