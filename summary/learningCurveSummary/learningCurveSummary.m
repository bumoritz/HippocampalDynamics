function [lcs] = learningCurveSummary(d_info,d,ops,p,path)

if ~exist([path.root_summary,'plots\learningCurves'],'dir')
    mkdir([path.root_summary,'plots\learningCurves']);
end

temp = d_info.excl==false;
d = fillEmptyFields(d,'perf',[nanmin(find(temp(:,1))),1],[nanmin(find(temp(:,2))),2]);


%% Learning curves - all animals

% all animals, individuals
dat = {}; plt = struct(); plt.ylim = [35,100];
for j=1:d_info.numDays
    dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{1:d_info.numAnimals,j}}','Uniform',0))*100;
end
plt.traces = 'individuals';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_all_ind.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_all_ind.png']);

% all animals, mean+-sem
dat = {}; plt = struct(); plt.ylim = [35,100];
for j=1:d_info.numDays
    dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{1:d_info.numAnimals,j}}','Uniform',0))*100;
end
plt.traces = 'mean+-sem';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_all_avg.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_all_avg.png']);

% all animals, mean+-sem, go vs nogo
dat = {}; plt = struct();
for j=1:d_info.numDays
    dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').H,{d{1:d_info.numAnimals,j}}','Uniform',0))*100;
    dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').CR,{d{1:d_info.numAnimals,j}}','Uniform',0))*100;
    dat{3,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{1:d_info.numAnimals,j}}','Uniform',0))*100;
    plt.col{1,j} = p.col.reward; plt.col{2,j} = p.col.gray; plt.col{3,j} = p.col.black;
end
plt.traces = 'mean+-sem';
F = learningCurves(dat,plt,p);         
savefig(F,[path.root_summary,'plots\learningCurves\learning_all_avg_gonogo.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_all_avg_gonogo.png']);


%% Learning curves - no-stim animals

% no-stim animals, individuals
dat = {}; plt = struct(); plt.ylim = [35,100];
for j=1:d_info.numDays
    dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==1)],j}}','Uniform',0))*100;
    dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==2)],j}}','Uniform',0))*100;
    if j==1
        dat{3,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==0);find(d_info.group==3);find(d_info.group==4);find(d_info.group==5);find(d_info.group==6);find(d_info.group==7);find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{3,j} = nan(size(dat{1,j}));
    end
    plt.col{1,j} = p.col.odour; plt.col{2,j} = p.col.imaging; plt.col{3,j} = nanmean([p.col.imaging;p.col.white],1);
end
plt.traces = 'individuals';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_nostim_ind.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_nostim_ind.png']);

% no-stim animals, mean+-sem
dat = {}; plt = struct(); plt.ylim = [35,100];
for j=1:d_info.numDays
    dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==1)],j}}','Uniform',0))*100;
    dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==2)],j}}','Uniform',0))*100;
    if j==1
        dat{3,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==0);find(d_info.group==3);find(d_info.group==4);find(d_info.group==5);find(d_info.group==6);find(d_info.group==7);find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{3,j} = nan(size(dat{1,j}));
    end
    dat{4,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{:,j}}','Uniform',0))*100;
    plt.col{1,j} = p.col.odour; plt.col{2,j} = p.col.imaging; plt.col{3,j} = nanmean([p.col.imaging;p.col.white],1); plt.col{4,j} = p.col.darkGray; 
end
plt.traces = 'mean+-sem';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_nostim_avg.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_nostim_avg.png']);

% no-stim animals, mean+-sem, go vs nogo
dat = {}; plt = struct();
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').H,{d{:,j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').CR,{d{:,j}}','Uniform',0))*100;
        dat{3,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{:,j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').H,{d{[find(d_info.group==0);find(d_info.group==1);find(d_info.group==2)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').CR,{d{[find(d_info.group==0);find(d_info.group==1);find(d_info.group==2)],j}}','Uniform',0))*100;
        dat{3,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==0);find(d_info.group==1);find(d_info.group==2)],j}}','Uniform',0))*100;
    end
    plt.col{1,j} = p.col.reward; plt.col{2,j} = p.col.gray; plt.col{3,j} = p.col.darkGray;
end
plt.traces = 'mean+-sem';
F = learningCurves(dat,plt,p);         
savefig(F,[path.root_summary,'plots\learningCurves\learning_nostim_avg_gonogo.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_nostim_avg_gonogo.png']);


%% Learning curves - graded speed

% all animals, individuals, graded speed
if ops.lcs.runningMetrics
    dat = {}; plt = struct(); dat_col = {}; plt.ylim = [35,100];
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{1:d_info.numAnimals,j}}','Uniform',0))*100;
        dat_col{1,j} = nan(d_info.numAnimals,1);
        for i=1:d_info.numAnimals
            if isfield(d{i,j},'paq_beh')
                dat_col{1,j}(i) = nanmean(d{i,j}.paq_beh.speed);
            end
        end
    end
    plt.traces = 'individuals';
    F = learningCurves(dat,plt,p,dat_col);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_all_ind_speed.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_all_ind_speed.png']);
end

% no-stim animals, individuals, graded speed
if ops.lcs.runningMetrics
    dat = {}; plt = struct(); dat_col = {}; plt.ylim = [35,100];
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==1);find(d_info.group==2)],j}}','Uniform',0))*100;
        temp = [find(d_info.group==1);find(d_info.group==2)];
        dat_col{1,j} = nan(length(temp),1);
        for ii=1:length(temp)
            i = temp(ii);
            if isfield(d{i,j},'paq_beh')
                dat_col{1,j}(ii) = nanmean(d{i,j}.paq_beh.speed);
            end
        end
    end
    plt.traces = 'individuals';
    F = learningCurves(dat,plt,p,dat_col);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_nostim_ind_speed.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_nostim_ind_speed.png']);
end

% imaging animals, individuals, graded speed
if ops.lcs.runningMetrics
    dat = {}; plt = struct(); dat_col = {}; plt.ylim = [35,100];
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==2)],j}}','Uniform',0))*100;
        temp = [find(d_info.group==2)];
        dat_col{1,j} = nan(length(temp),1);
        for ii=1:length(temp)
            i = temp(ii);
            if isfield(d{i,j},'paq_beh')
                dat_col{1,j}(ii) = nanmean(d{i,j}.paq_beh.speed);
            end
        end
    end
    plt.traces = 'individuals';
    F = learningCurves(dat,plt,p,dat_col);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_imaging_ind_speed.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_imaging_ind_speed.png']);
end

% stim animals, individuals, graded speed
if ops.lcs.runningMetrics
    dat = {}; plt = struct(); dat_col = {}; plt.ylim = [35,100];
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==7);find(d_info.group==8)],j}}','Uniform',0))*100;
        temp = [find(d_info.group==7);find(d_info.group==8)];
        dat_col{1,j} = nan(length(temp),1);
        for ii=1:length(temp)
            i = temp(ii);
            if isfield(d{i,j},'paq_beh')
                dat_col{1,j}(ii) = nanmean(d{i,j}.paq_beh.speed);
            end
        end
    end
    plt.traces = 'individuals';
    F = learningCurves(dat,plt,p,dat_col);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_ind_speed.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_ind_speed.png']);
end


%% Learning curves - runner vs non-runner

% all animals, individuals, runner vs non-runner
if ops.lcs.runningMetrics
    dat = {}; plt = struct(); plt.ylim = [35,100];
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{find(d_info.running(:,j)==0),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{find(d_info.running(:,j)==1),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonrunner; plt.col{2,j} = p.col.runner;
    end
    plt.traces = 'individuals';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_all_ind_running.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_all_ind_running.png']);
end

% all animals, mean+-sem, runner vs non-runner
if ops.lcs.runningMetrics
    dat = {}; plt = struct(); plt.ylim = [35,100];
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{find(d_info.running(:,j)==0),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{find(d_info.running(:,j)==1),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonrunner; plt.col{2,j} = p.col.runner;
    end
    plt.traces = 'mean+-sem';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_all_avg_running.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_all_avg_running.png']);
end

% no-stim animals, individuals, runner vs non-runner
if ops.lcs.runningMetrics
    dat = {}; plt = struct(); plt.ylim = [35,100];
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.running(:,j)==0),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.running(:,j)==1),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonrunner; plt.col{2,j} = p.col.runner;
    end
    plt.traces = 'individuals';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_nostim_ind_running.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_nostim_ind_running.png']);
end

% no-stim animals, mean+-sem, runner vs non-runner
if ops.lcs.runningMetrics
    dat = {}; plt = struct(); plt.ylim = [35,100];
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.running(:,j)==0),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.running(:,j)==1),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonrunner; plt.col{2,j} = p.col.runner;
    end
    plt.traces = 'mean+-sem';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_nostim_avg_running.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_nostim_avg_running.png']);
end

% imaging animals, individuals, runner vs non-runner
if ops.lcs.runningMetrics
    dat = {}; plt = struct(); plt.ylim = [35,100];
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.running(:,j)==0),[find(d_info.group==2)]),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.running(:,j)==1),[find(d_info.group==2)]),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonrunner; plt.col{2,j} = p.col.runner;
    end
    plt.traces = 'individuals';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_imaging_ind_running.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_imaging_ind_running.png']);
end

% imaging animals, mean+-sem, runner vs non-runner
if ops.lcs.runningMetrics
    dat = {}; plt = struct(); plt.ylim = [35,100];
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.running(:,j)==0),[find(d_info.group==2)]),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.running(:,j)==1),[find(d_info.group==2)]),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonrunner; plt.col{2,j} = p.col.runner;
    end
    plt.traces = 'mean+-sem';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_imaging_avg_running.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_imaging_avg_running.png']);
end


%% Learning curves - learning vs non-learning

% no-stim animals, individuals, learning vs non-learning
if ops.lcs.runningMetrics
    dat = {}; plt = struct(); plt.ylim = [35,100];
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.learning(:,j)==0),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.learning(:,j)==1),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonlearning; plt.col{2,j} = p.col.learning;
    end
    plt.traces = 'individuals';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_nostim_ind_learning.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_nostim_ind_learning.png']);
end

% no-stim animals, mean+-sem, learning vs non-learning
if ops.lcs.runningMetrics
    dat = {}; plt = struct(); plt.ylim = [35,100];
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.learning(:,j)==0),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.learning(:,j)==1),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonlearning; plt.col{2,j} = p.col.learning;
    end
    plt.traces = 'mean+-sem';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_nostim_avg_learning.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_nostim_avg_learning.png']);
end


%% Learning curves - stim animals

% stim animals, individuals
dat = {}; plt = struct(); plt.ylim = [35,100];
for j=1:d_info.numDays
    dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
    dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    if j==1
        plt.col{1,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{2,j} = nanmean([p.col.ctrl;p.col.gray],1);
    elseif j==2 || j==3
        plt.col{1,j} = p.col.seq; plt.col{2,j} = p.col.ctrl;
    elseif j==4 || j==5
        plt.col{1,j} = p.col.ctrl; plt.col{2,j} = p.col.seq;
    end
end
plt.traces = 'individuals';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_ind.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_ind.png']);

% stim animals, mean+-sem
dat = {}; plt = struct(); plt.ylim = [35,100];
for j=1:d_info.numDays
    dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
    dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    if j==1
        plt.col{1,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{2,j} = nanmean([p.col.ctrl;p.col.gray],1);
    elseif j==2 || j==3
        plt.col{1,j} = p.col.seq; plt.col{2,j} = p.col.ctrl;
    elseif j==4 || j==5
        plt.col{1,j} = p.col.ctrl; plt.col{2,j} = p.col.seq;
    end
end
plt.traces = 'mean+-sem';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_avg.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_avg.png']);


%% Learning curves - stim animals (e1p1)

% stim animals (e1p1), individuals
dat = {}; plt = struct(); plt.ylim = [35,100];
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
    end
    if j==1
        plt.col{1,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{2,j} = nanmean([p.col.ctrl;p.col.gray],1);
    elseif j==2 || j==3
        plt.col{1,j} = p.col.seq; plt.col{2,j} = p.col.ctrl;
    elseif j==4 || j==5
        plt.col{1,j} = p.col.ctrl; plt.col{2,j} = p.col.seq;
    end
end
plt.traces = 'individuals';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_ind.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_ind.png']);

% stim animals (e1p1), mean+-sem
dat = {}; plt = struct(); plt.ylim = [35,100];
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
    end
    if j==1
        plt.col{1,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{2,j} = nanmean([p.col.ctrl;p.col.gray],1);
    elseif j==2 || j==3
        plt.col{1,j} = p.col.seq; plt.col{2,j} = p.col.ctrl;
    elseif j==4 || j==5
        plt.col{1,j} = p.col.ctrl; plt.col{2,j} = p.col.seq;
    end
end
plt.traces = 'mean+-sem';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_avg.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_avg.png']);

% stim animals (e1p1), individuals, switches pooled
dat_switches = {};
dat_switches{1,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq
dat_switches{2,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
plt_switches = plt; plt_switches.col = {};
plt_switches.col{1,1} = p.col.seq; plt_switches.col{1,2} = p.col.seq; % seq
plt_switches.col{2,1} = p.col.ctrl; plt_switches.col{2,2} = p.col.ctrl; % ctrl
plt_switches.traces = 'individuals';
F = learningCurves(dat_switches,plt_switches,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_ind_switches.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_ind_switches.png']);

% stim animals (e1p1), mean+-sem
plt_switches.traces = 'mean+-sem';
F = learningCurves(dat_switches,plt_switches,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_avg_switches.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_avg_switches.png']);


%% Learning curves - stim animals (e1p1) - stim trials

% stim animals (e1p1) - stim trials, individuals
dat = {}; plt = struct(); plt.ylim = [35,100];
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
    end
    if j==1
        plt.col{1,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{2,j} = nanmean([p.col.ctrl;p.col.gray],1);
    elseif j==2 || j==3
        plt.col{1,j} = p.col.seq; plt.col{2,j} = p.col.ctrl;
    elseif j==4 || j==5
        plt.col{1,j} = p.col.ctrl; plt.col{2,j} = p.col.seq;
    end
end
plt.traces = 'individuals';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stim_ind.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stim_ind.png']);

% stim animals (e1p1) - stim trials, mean+-sem
dat = {}; plt = struct(); plt.ylim = [35,100];
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
    end
    if j==1
        plt.col{1,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{2,j} = nanmean([p.col.ctrl;p.col.gray],1);
    elseif j==2 || j==3
        plt.col{1,j} = p.col.seq; plt.col{2,j} = p.col.ctrl;
    elseif j==4 || j==5
        plt.col{1,j} = p.col.ctrl; plt.col{2,j} = p.col.seq;
    end
end
plt.traces = 'mean+-sem';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stim_avg.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stim_avg.png']);

% stim animals (e1p1) - stim trials, individuals, switches pooled
dat_switches = {};
dat_switches{1,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq
dat_switches{2,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
plt_switches = plt; plt_switches.col = {};
plt_switches.col{1,1} = p.col.seq; plt_switches.col{1,2} = p.col.seq; % seq
plt_switches.col{2,1} = p.col.ctrl; plt_switches.col{2,2} = p.col.ctrl; % ctrl
plt_switches.traces = 'individuals';
F = learningCurves(dat_switches,plt_switches,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stim_ind_switches.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stim_ind_switches.png']);

% stim animals (e1p1) - stim trials, mean+-sem
plt_switches.traces = 'mean+-sem';
F = learningCurves(dat_switches,plt_switches,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stim_avg_switches.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stim_avg_switches.png']);


%% Learning curves - stim animals (e1p1) - catch trials

% stim animals (e1p1) - catch trials, individuals
dat = {}; plt = struct(); plt.ylim = [35,100];
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
    end
    if j==1
        plt.col{1,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{2,j} = nanmean([p.col.ctrl;p.col.gray],1);
    elseif j==2 || j==3
        plt.col{1,j} = p.col.seq; plt.col{2,j} = p.col.ctrl;
    elseif j==4 || j==5
        plt.col{1,j} = p.col.ctrl; plt.col{2,j} = p.col.seq;
    end
end
plt.traces = 'individuals';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_catch_ind.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_catch_ind.png']);

% stim animals (e1p1) - catch trials, mean+-sem
dat = {}; plt = struct(); plt.ylim = [35,100];
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
    end
    if j==1
        plt.col{1,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{2,j} = nanmean([p.col.ctrl;p.col.gray],1);
    elseif j==2 || j==3
        plt.col{1,j} = p.col.seq; plt.col{2,j} = p.col.ctrl;
    elseif j==4 || j==5
        plt.col{1,j} = p.col.ctrl; plt.col{2,j} = p.col.seq;
    end
end
plt.traces = 'mean+-sem';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_catch_avg.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_catch_avg.png']);

% stim animals (e1p1) - catch trials, individuals, switches pooled
dat_switches = {};
dat_switches{1,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq
dat_switches{2,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
plt_switches = plt; plt_switches.col = {};
plt_switches.col{1,1} = p.col.seq; plt_switches.col{1,2} = p.col.seq; % seq
plt_switches.col{2,1} = p.col.ctrl; plt_switches.col{2,2} = p.col.ctrl; % ctrl
plt_switches.traces = 'individuals';
F = learningCurves(dat_switches,plt_switches,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_catch_ind_switches.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_catch_ind_switches.png']);

% stim animals (e1p1) - catch trials, mean+-sem
plt_switches.traces = 'mean+-sem';
F = learningCurves(dat_switches,plt_switches,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_catch_avg_switches.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_catch_avg_switches.png']);


%% Learning curves - stim animals (e1p1) - catch vs stim trials

% stim animals (e1p1) - catch trials, individuals
dat = {}; plt = struct(); plt.ylim = [35,100];
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
        dat{3,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{4,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
        dat{3,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
        dat{4,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
    end
    if j==1
        plt.col{1,j} = p.col.gray; plt.col{2,j} = p.col.gray; plt.col{3,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{4,j} = nanmean([p.col.ctrl;p.col.gray],1);
    elseif j==2 || j==3
        plt.col{1,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{2,j} = nanmean([p.col.ctrl;p.col.gray],1); plt.col{3,j} = p.col.seq; plt.col{4,j} = p.col.ctrl;
    elseif j==4 || j==5
        plt.col{1,j} = nanmean([p.col.ctrl;p.col.gray],1); plt.col{2,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{3,j} = p.col.ctrl; plt.col{4,j} = p.col.seq;
    end
end
plt.traces = 'individuals';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stimcatch_ind.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stimccatch_ind.png']);

% stim animals (e1p1) - catch trials, mean+-sem
dat = {}; plt = struct(); plt.ylim = [35,100];
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
        dat{3,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{4,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
        dat{3,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
        dat{4,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
    end
    if j==1
        plt.col{1,j} = p.col.gray; plt.col{2,j} = p.col.gray; plt.col{3,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{4,j} = nanmean([p.col.ctrl;p.col.gray],1);
    elseif j==2 || j==3
        plt.col{1,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{2,j} = nanmean([p.col.ctrl;p.col.gray],1); plt.col{3,j} = p.col.seq; plt.col{4,j} = p.col.ctrl;
    elseif j==4 || j==5
        plt.col{1,j} = nanmean([p.col.ctrl;p.col.gray],1); plt.col{2,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{3,j} = p.col.ctrl; plt.col{4,j} = p.col.seq;
    end
end
plt.traces = 'mean+-sem';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stimcatch_avg.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stimcatch_avg.png']);

% stim animals (e1p1) - catch trials, individuals, switches pooled
dat_switches = {};
dat_switches{1,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq
dat_switches{2,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
dat_switches{3,1} = [dat{3,2};dat{4,4}]; dat_switches{3,2} = [dat{3,3};dat{4,5}]; % seq
dat_switches{4,1} = [dat{4,2};dat{3,4}]; dat_switches{4,2} = [dat{4,3};dat{3,5}]; % ctrl
plt_switches = plt; plt_switches.col = {};
plt_switches.col{1,1} = p.col.seq; plt_switches.col{1,2} = p.col.seq; % seq
plt_switches.col{2,1} = p.col.ctrl; plt_switches.col{2,2} = p.col.ctrl; % ctrl
plt_switches.col{3,1} = nanmean([p.col.seq;p.col.gray],1); plt_switches.col{3,2} = nanmean([p.col.seq;p.col.gray],1); % seq
plt_switches.col{4,1} = nanmean([p.col.ctrl;p.col.gray],1); plt_switches.col{4,2} = nanmean([p.col.ctrl;p.col.gray],1); % ctrl
plt_switches.traces = 'individuals';
F = learningCurves(dat_switches,plt_switches,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stimcatch_ind_switches.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stimcatch_ind_switches.png']);

% stim animals (e1p1) - catch trials, mean+-sem
plt_switches.traces = 'mean+-sem';
F = learningCurves(dat_switches,plt_switches,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stimcatch_avg_switches.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_stimcatch_avg_switches.png']);


%% Learning curves - stim animals (e1p1) - delta

% stim animals (e1p1) - delta, individuals
dat = {}; plt = struct(); plt.ylim = [-20,20];
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==7)],j}}','Uniform',0))*100 - cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==8)],j}}','Uniform',0))*100 - cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100 - cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100 - cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
    end
    if j==1
        plt.col{1,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{2,j} = nanmean([p.col.ctrl;p.col.gray],1);
    elseif j==2 || j==3
        plt.col{1,j} = p.col.seq; plt.col{2,j} = p.col.ctrl;
    elseif j==4 || j==5
        plt.col{1,j} = p.col.ctrl; plt.col{2,j} = p.col.seq;
    end
end
plt.traces = 'individuals';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_delta_ind.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_delta_ind.png']);

% stim animals (e1p1) - delta, mean+-sem
dat = {}; plt = struct(); plt.ylim = [-20,20];
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==7)],j}}','Uniform',0))*100 - cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==8)],j}}','Uniform',0))*100 - cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100 - cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100 - cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
    end
    if j==1
        plt.col{1,j} = nanmean([p.col.seq;p.col.gray],1); plt.col{2,j} = nanmean([p.col.ctrl;p.col.gray],1);
    elseif j==2 || j==3
        plt.col{1,j} = p.col.seq; plt.col{2,j} = p.col.ctrl;
    elseif j==4 || j==5
        plt.col{1,j} = p.col.ctrl; plt.col{2,j} = p.col.seq;
    end
end
plt.traces = 'mean+-sem';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_delta_avg.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_delta_avg.png']);

% stim animals (e1p1) - delta, individuals, switches pooled
dat_switches = {};
dat_switches{1,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq
dat_switches{2,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
plt_switches = plt; plt_switches.col = {};
plt_switches.col{1,1} = p.col.seq; plt_switches.col{1,2} = p.col.seq; % seq
plt_switches.col{2,1} = p.col.ctrl; plt_switches.col{2,2} = p.col.ctrl; % ctrl
plt_switches.traces = 'individuals';
F = learningCurves(dat_switches,plt_switches,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_delta_ind_switches.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_delta_ind_switches.png']);

% stim animals (e1p1) - delta, mean+-sem
plt_switches.traces = 'mean+-sem';
F = learningCurves(dat_switches,plt_switches,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_delta_avg_switches.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_stim_e1p1_delta_avg_switches.png']);


