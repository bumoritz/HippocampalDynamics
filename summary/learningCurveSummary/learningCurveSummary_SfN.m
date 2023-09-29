%% SfN

% learningCurveSummary
% d_info = selectDataset(d_info,'-g12-d12345-e01-r01-p01-l01-i00',sheet,path,ops); % for runner vs non-runner
% d_info = selectDataset(d_info,'-g0123456789101112-d12345-e01-r01-p01-l01-i00',sheet,path,ops); % for general
% d_info = d_info = selectDataset(d_info,'-g0123456789101112-d12345-e1-r01-p01-l01-i00',sheet,path,ops); % for general


%% Preparations

if ~exist([path.root_summary,'plots\learningCurves_SfN'],'dir')
    mkdir([path.root_summary,'plots\learningCurves_SfN']);
end

temp = d_info.excl==false;
d = fillEmptyFields(d,'perf',[nanmin(find(temp(:,1))),1],[nanmin(find(temp(:,2))),2]);


%% Learning curves - no-stim animals

% no-stim animals, individuals
dat = {}; plt = struct(); plt.ylim = [35,100]; plt.flag = 'SfN'; plt.xlabel = 'Trial block (in blocks of 100 trials)';
for j=1:d_info.numDays
    dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==1)],j}}','Uniform',0))*100;
    dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==2)],j}}','Uniform',0))*100;
    if j==1
        dat{3,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==0);find(d_info.group==3);find(d_info.group==4);find(d_info.group==5);find(d_info.group==6);find(d_info.group==7);find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{3,j} = nan(size(dat{1,j}));
    end
    %plt.col{1,j} = p.col.odour; plt.col{2,j} = p.col.imaging; plt.col{3,j} = nanmean([p.col.imaging;p.col.white],1);
    plt.col{1,j} = p.col.black; plt.col{2,j} = p.col.black; plt.col{3,j} = p.col.black;
end
plt.traces = 'individuals';
F = learningCurves(dat,plt,p);
savefig(F,[path.root_summary,'plots\learningCurves_SfN\learning_nostim_ind.fig']);
saveas(F,[path.root_summary,'plots\learningCurves_SfN\learning_nostim_ind.png']);

% no-stim animals, mean+-sem, go vs nogo
dat = {}; plt = struct(); plt.flag = 'SfN'; plt.xlabel = 'Trial block (in blocks of 100 trials)';
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
savefig(F,[path.root_summary,'plots\learningCurves_SfN\learning_nostim_avg_gonogo.fig']);
saveas(F,[path.root_summary,'plots\learningCurves_SfN\learning_nostim_avg_gonogo.png']);


%% Learning curves - runner vs non-runner

% no-stim animals, individuals, runner vs non-runner
if ops.lcs.runningMetrics
    dat = {}; plt = struct(); plt.ylim = [35,100]; plt.flag = 'SfN'; plt.xlabel = 'Trial block (in blocks of 100 trials)';
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.running(:,j)==0),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.running(:,j)==1),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonrunner; plt.col{2,j} = p.col.runner;
    end
    plt.traces = 'individuals';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves_SfN\learning_nostim_ind_running.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves_SfN\learning_nostim_ind_running.png']);
end

% no-stim animals, mean+-sem, runner vs non-runner
if ops.lcs.runningMetrics
    dat = {}; plt = struct(); plt.ylim = [35,100]; plt.flag = 'SfN'; plt.xlabel = 'Trial block (in blocks of 100 trials)';
    for j=1:d_info.numDays
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.running(:,j)==0),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{intersect(find(d_info.running(:,j)==1),[find(d_info.group==1);find(d_info.group==2)]),j}}','Uniform',0))*100;
        plt.col{1,j} = p.col.nonrunner; plt.col{2,j} = p.col.runner;
    end
    plt.traces = 'mean+-sem';
    F = learningCurves(dat,plt,p);
    savefig(F,[path.root_summary,'plots\learningCurves_SfN\learning_nostim_avg_running.fig']);
    saveas(F,[path.root_summary,'plots\learningCurves_SfN\learning_nostim_avg_running.png']);
end


%% Learning curves - stim animals (e1p1)

% stim animals (e1p1), individuals
dat = {}; plt = struct(); plt.ylim = [35,100]; plt.flag = 'SfN'; plt.xlabel = 'Trial block (in blocks of 100 trials)';
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
savefig(F,[path.root_summary,'plots\learningCurves_SfN\learning_stim_e1p1_ind.fig']);
saveas(F,[path.root_summary,'plots\learningCurves_SfN\learning_stim_e1p1_ind.png']);

% stim animals (e1p1), mean+-sem
dat = {}; plt = struct(); plt.ylim = [35,100]; plt.flag = 'SfN'; plt.xlabel = 'Trial block (in blocks of 100 trials)';
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
savefig(F,[path.root_summary,'plots\learningCurves_SfN\learning_stim_e1p1_avg.fig']);
saveas(F,[path.root_summary,'plots\learningCurves_SfN\learning_stim_e1p1_avg.png']);

% stim animals (e1p1), individuals, switches pooled
dat_switches = {}; 
dat_switches{1,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq
dat_switches{2,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
plt_switches = plt; plt_switches.col = {}; plt.flag = 'SfN'; plt.xlabel = 'Trial block (in blocks of 100 trials)';
plt_switches.col{1,1} = p.col.seq; plt_switches.col{1,2} = p.col.seq; % seq
plt_switches.col{2,1} = p.col.ctrl; plt_switches.col{2,2} = p.col.ctrl; % ctrl
plt_switches.traces = 'individuals';
F = learningCurves(dat_switches,plt_switches,p);
savefig(F,[path.root_summary,'plots\learningCurves_SfN\learning_stim_e1p1_ind_switches.fig']);
saveas(F,[path.root_summary,'plots\learningCurves_SfN\learning_stim_e1p1_ind_switches.png']);

% stim animals (e1p1), mean+-sem
plt_switches.traces = 'mean+-sem';
F = learningCurves(dat_switches,plt_switches,p);
savefig(F,[path.root_summary,'plots\learningCurves_SfN\learning_stim_e1p1_avg_switches.fig']);
saveas(F,[path.root_summary,'plots\learningCurves_SfN\learning_stim_e1p1_avg_switches.png']);

