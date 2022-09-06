function [lcs] = learningCurveSummary(d_info,d,ops,p,path)

if ~exist([path.root_summary,'plots\learningCurves'],'dir')
    mkdir([path.root_summary,'plots\learningCurves']);
end

d = fillEmptyFields(d,'perf');


%% Build binary running vector (temporary)

if ops.lcs.runningMetrics
    running = nan(d_info.numAnimals,1);
    for i=1:d_info.numAnimals
        if isfield(d{i,1},'paq_beh') && (~isnan(nanmean(d{i,1}.paq_beh.speed)))
            running(i) = mean(d{i,1}.paq_beh.speed) > 20;
        end
    end
end
running(16) = true; % Arwen
running(27) = false; % Stanage
running(36:42) = 0 % Python, correct after data is imported


%% Learning curves - all animals

% all animals, individuals
dat = {}; plt = struct();
for j=1:d_info.numDays
    dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{1:d_info.numAnimals,j}}','Uniform',0))*100;
end
plt.traces = 'individuals';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves\learning_all_ind.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_all_ind.png']);

% all animals, mean+-sem
dat = {}; plt = struct();
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
dat = {}; plt = struct();
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
dat = {}; plt = struct();
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


%% Learning curves - stim animals

% stim animals, individuals
dat = {}; plt = struct();
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
savefig(F,[path.root_summary,'plots\learningCurves\learning_avg_ind.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_avg_ind.png']);

% stim animals, mean+-sem
dat = {}; plt = struct();
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
savefig(F,[path.root_summary,'plots\learningCurves\learning_avg_ind.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\learning_avg_ind.png']);


%% Learning curves - graded speed

% all animals, individuals, graded speed
if ops.lcs.runningMetrics
    dat = {}; plt = struct(); dat_col = {};
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
    dat = {}; plt = struct(); dat_col = {};
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
    dat = {}; plt = struct(); dat_col = {};
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
    dat = {}; plt = struct(); dat_col = {};
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
    dat = {}; plt = struct();
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{find(~running),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{find(running),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonrunner; plt.col{2,j} = p.col.runner;
    end
    plt.traces = 'individuals';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_all_ind_running.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_all_ind_running.png']);
end

% all animals, mean+-sem, runner vs non-runner
if ops.lcs.runningMetrics
    dat = {}; plt = struct();
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{find(~running),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{find(running),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonrunner; plt.col{2,j} = p.col.runner;
    end
    plt.traces = 'mean+-sem';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_all_avg_running.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_all_avg_running.png']);
end

% no-stim animals, individuals, runner vs non-runner
if ops.lcs.runningMetrics
    dat = {}; plt = struct();
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(~running),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(running),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonrunner; plt.col{2,j} = p.col.runner;
    end
    plt.traces = 'individuals';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_nostim_ind_running.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_nostim_ind_running.png']);
end

% no-stim animals, mean+-sem, runner vs non-runner
if ops.lcs.runningMetrics
    dat = {}; plt = struct();
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(~running),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(running),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonrunner; plt.col{2,j} = p.col.runner;
    end
    plt.traces = 'mean+-sem';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_nostim_avg_running.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_nostim_avg_running.png']);
end

% imaging animals, individuals, runner vs non-runner
if ops.lcs.runningMetrics
    dat = {}; plt = struct();
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(~running),[find(d_info.group==2)]),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(running),[find(d_info.group==2)]),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonrunner; plt.col{2,j} = p.col.runner;
    end
    plt.traces = 'individuals';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_imaging_ind_running.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_imaging_ind_running.png']);
end

% imaging animals, mean+-sem, runner vs non-runner
if ops.lcs.runningMetrics
    dat = {}; plt = struct();
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(~running),[find(d_info.group==2)]),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(running),[find(d_info.group==2)]),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonrunner; plt.col{2,j} = p.col.runner;
    end
    plt.traces = 'mean+-sem';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves\learning_imaging_avg_running.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves\learning_imaging_avg_running.png']);
end


%% TEMP

% no-stim animals, mean+-sem, go vs nogo
dat = {}; plt = struct(); 
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{:,j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==0);find(d_info.group==1);find(d_info.group==2)],j}}','Uniform',0))*100;
    end
    plt.col{1,j} = p.col.darkGray;
end
plt.traces = 'mean+-sem';
F = learningCurves(dat,plt,p);         
savefig(F,[path.root_summary,'plots\learningCurves\TEMP1.fig']);
saveas(F,[path.root_summary,'plots\learningCurves\TEMP1.png']);


%end