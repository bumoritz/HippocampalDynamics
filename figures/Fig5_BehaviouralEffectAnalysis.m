%% Fig5_BehaviouralEffectAnalysis

% import data using Summary_Master with ops.do_learningCurveSummary = true;
% d_info = selectDataset(d_info,'-g78-d2345-e01-r01-p01-l01-i00',sheet,path,ops);

save_root_fig = [path.root_summary,'figures\Fig5_fig\'];
save_root_png = [path.root_summary,'figures\Fig5_png\'];
save_root_pdf = [path.root_summary,'figures\Fig5_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig5_txt\'];


%% Get data

temp = d_info.excl==false;
d = fillEmptyFields(d,'perf',[nanmin(find(temp(:,2))),5],[nanmin(find(temp(:,2))),2]);

    
%% --- Main figure ---

%% Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r01

% get data
dat = {};
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),find(d_info.presponsive(:,j)==1)))],j}}','Uniform',0))*100;
    end
end
dat_switches = {}; 
dat_switches{1,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
dat_switches{2,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq

% make plot
plt_switches = struct(); plt_switches.ylim = [35,100]; plt_switches.flag = 'paper'; plt_switches.traces = 'individuals and mean'; plt_switches.xlabel = 'Trial block'; plt_switches.col = {};
plt_switches.col{1,1} = p.col.ctrl; plt_switches.col{1,2} = p.col.ctrl;
plt_switches.col{2,1} = p.col.seq; plt_switches.col{2,2} = p.col.seq;
F = learningCurves(dat_switches,plt_switches,p);

% save plot
savefig(F,[save_root_fig,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r01.fig']);
saveas(F,[save_root_png,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r01.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r01.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r01.txt'],'wt');
fprintf(fid,['\nBehaviouralEffectAnalysis_LearningCurve_p1e1r01\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']); fclose(fid);
fprintf(['\nBehaviouralEffectAnalysis_LearningCurve_p1e1r01\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']);


%% --- Supplementary figure ---

%% --- Reserve ---

%% Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r01

% get data
dat = {};
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==7),find(d_info.presponsive(:,j)==1))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==8),find(d_info.presponsive(:,j)==1))],j}}','Uniform',0))*100;    
    end
end
dat_switches = {}; 
dat_switches{1,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
dat_switches{2,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq

% make plot
plt_switches = struct(); plt_switches.ylim = [35,100]; plt_switches.flag = 'paper'; plt_switches.traces = 'individuals and mean'; plt_switches.xlabel = 'Trial block'; plt_switches.col = {};
plt_switches.col{1,1} = p.col.ctrl; plt_switches.col{1,2} = p.col.ctrl;
plt_switches.col{2,1} = p.col.seq; plt_switches.col{2,2} = p.col.seq;
F = learningCurves(dat_switches,plt_switches,p);

% save plot
savefig(F,[save_root_fig,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r01.fig']);
saveas(F,[save_root_png,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r01.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r01.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r01.txt'],'wt');
fprintf(fid,['\nBehaviouralEffectAnalysis_LearningCurve_p1e01r01\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']); fclose(fid);
fprintf(['\nBehaviouralEffectAnalysis_LearningCurve_p1e01r01\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']);


%% Fig5_BehaviouralEffectAnalysis_LearningCurve_p01e01r01

% get data
dat = {};
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;    
    end
end
dat_switches = {}; 
dat_switches{1,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
dat_switches{2,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq

% make plot
plt_switches = struct(); plt_switches.ylim = [35,100]; plt_switches.flag = 'paper'; plt_switches.traces = 'individuals and mean'; plt_switches.xlabel = 'Trial block'; plt_switches.col = {};
plt_switches.col{1,1} = p.col.ctrl; plt_switches.col{1,2} = p.col.ctrl;
plt_switches.col{2,1} = p.col.seq; plt_switches.col{2,2} = p.col.seq;
F = learningCurves(dat_switches,plt_switches,p);

% save plot
savefig(F,[save_root_fig,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p01e01r01.fig']);
saveas(F,[save_root_png,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p01e01r01.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p01e01r01.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p01e01r01.txt'],'wt');
fprintf(fid,['\nBehaviouralEffectAnalysis_LearningCurve_p01e01r01\nn(d1_ctrl)=',num2str(size(rmmissing(dat_switches{1,1},'MinNumMissing',5),1)),'\nn(d2_ctrl)=',num2str(size(rmmissing(dat_switches{1,2},'MinNumMissing',5),1)),'\nn(d1_seq)=',num2str(size(rmmissing(dat_switches{2,1},'MinNumMissing',5),1)),'\nn(d2_seq)=',num2str(size(rmmissing(dat_switches{2,2},'MinNumMissing',5),1)),'\n']); fclose(fid);
fprintf(['\nBehaviouralEffectAnalysis_LearningCurve_p01e01r01\nn(d1_ctrl)=',num2str(size(rmmissing(dat_switches{1,1},'MinNumMissing',5),1)),'\nn(d2_ctrl)=',num2str(size(rmmissing(dat_switches{1,2},'MinNumMissing',5),1)),'\nn(d1_seq)=',num2str(size(rmmissing(dat_switches{2,1},'MinNumMissing',5),1)),'\nn(d2_seq)=',num2str(size(rmmissing(dat_switches{2,2},'MinNumMissing',5),1)),'\n']);


%% Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r1 - faulty original

% get data
dat = {};
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1))))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1))))],j}}','Uniform',0))*100;
    end
end
dat_switches = {}; 
dat_switches{1,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
dat_switches{2,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq

% make plot
plt_switches = struct(); plt_switches.ylim = [35,100]; plt_switches.flag = 'paper'; plt_switches.traces = 'individuals and mean'; plt_switches.xlabel = 'Trial block'; plt_switches.col = {};
plt_switches.col{1,1} = p.col.ctrl; plt_switches.col{1,2} = p.col.ctrl;
plt_switches.col{2,1} = p.col.seq; plt_switches.col{2,2} = p.col.seq;
F = learningCurves(dat_switches,plt_switches,p);

% save plot
savefig(F,[save_root_fig,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r1.fig']);
saveas(F,[save_root_png,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r1.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r1.txt'],'wt');
fprintf(fid,['\nBehaviouralEffectAnalysis_LearningCurve_p1e1r1\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']); fclose(fid);
fprintf(['\nBehaviouralEffectAnalysis_LearningCurve_p1e1r1\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']);


%% Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r1 - corrected

% get data
dat = {};
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1))))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1))))],j}}','Uniform',0))*100;
    end
end
dat_switches = {}; 
dat_switches{1,1} = [dat{2,2};dat{1,4}]; dat_switches{1,2} = [dat{2,3};dat{1,5}]; % ctrl
dat_switches{2,1} = [dat{1,2};dat{2,4}]; dat_switches{2,2} = [dat{1,3};dat{2,5}]; % seq

% make plot
plt_switches = struct(); plt_switches.ylim = [35,100]; plt_switches.flag = 'paper'; plt_switches.traces = 'individuals and mean'; plt_switches.xlabel = 'Trial block'; plt_switches.col = {};
plt_switches.col{1,1} = p.col.ctrl; plt_switches.col{1,2} = p.col.ctrl;
plt_switches.col{2,1} = p.col.seq; plt_switches.col{2,2} = p.col.seq;
F = learningCurves(dat_switches,plt_switches,p);

% save plot
savefig(F,[save_root_fig,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r1.fig']);
saveas(F,[save_root_png,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r1.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r1.txt'],'wt');
fprintf(fid,['\nBehaviouralEffectAnalysis_LearningCurve_p1e1r1\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']); fclose(fid);
fprintf(['\nBehaviouralEffectAnalysis_LearningCurve_p1e1r1\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']);



%% Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r0 - faulty original

% get data
dat = {};
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==0))))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==0))))],j}}','Uniform',0))*100;
    end
end
dat_switches = {}; 
dat_switches{1,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
dat_switches{2,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq

% make plot
plt_switches = struct(); plt_switches.ylim = [35,100]; plt_switches.flag = 'paper'; plt_switches.traces = 'individuals and mean'; plt_switches.xlabel = 'Trial block'; plt_switches.col = {};
plt_switches.col{1,1} = p.col.ctrl; plt_switches.col{1,2} = p.col.ctrl;
plt_switches.col{2,1} = p.col.seq; plt_switches.col{2,2} = p.col.seq;
F = learningCurves(dat_switches,plt_switches,p);

% save plot
savefig(F,[save_root_fig,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r0.fig']);
saveas(F,[save_root_png,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r0.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r0.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r0.txt'],'wt');
fprintf(fid,['\nBehaviouralEffectAnalysis_LearningCurve_p1e1r0\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']); fclose(fid);
fprintf(['\nBehaviouralEffectAnalysis_LearningCurve_p1e1r0\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']);


%% Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r0 - corrected

% get data
dat = {};
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==7),intersect(find(d_info.engagement(:,j)==1),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==0))))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==8),intersect(find(d_info.engagement(:,j)==1),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==0))))],j}}','Uniform',0))*100;
    end
end
dat_switches = {}; 
dat_switches{1,1} = [dat{2,2};dat{1,4}]; dat_switches{1,2} = [dat{2,3};dat{1,5}]; % ctrl
dat_switches{2,1} = [dat{1,2};dat{2,4}]; dat_switches{2,2} = [dat{1,3};dat{2,5}]; % seq

% make plot
plt_switches = struct(); plt_switches.ylim = [35,100]; plt_switches.flag = 'paper'; plt_switches.traces = 'individuals and mean'; plt_switches.xlabel = 'Trial block'; plt_switches.col = {};
plt_switches.col{1,1} = p.col.ctrl; plt_switches.col{1,2} = p.col.ctrl;
plt_switches.col{2,1} = p.col.seq; plt_switches.col{2,2} = p.col.seq;
F = learningCurves(dat_switches,plt_switches,p);

% save plot
savefig(F,[save_root_fig,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r0.fig']);
saveas(F,[save_root_png,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r0.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r0.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e1r0.txt'],'wt');
fprintf(fid,['\nBehaviouralEffectAnalysis_LearningCurve_p1e1r0\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']); fclose(fid);
fprintf(['\nBehaviouralEffectAnalysis_LearningCurve_p1e1r0\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']);



%% Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r0

% get data
dat = {};
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==7),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==0)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==8),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==0)))],j}}','Uniform',0))*100;
    end
end
dat_switches = {}; 
dat_switches{1,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
dat_switches{2,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq

% make plot
plt_switches = struct(); plt_switches.ylim = [35,100]; plt_switches.flag = 'paper'; plt_switches.traces = 'individuals and mean'; plt_switches.xlabel = 'Trial block'; plt_switches.col = {};
plt_switches.col{1,1} = p.col.ctrl; plt_switches.col{1,2} = p.col.ctrl;
plt_switches.col{2,1} = p.col.seq; plt_switches.col{2,2} = p.col.seq;
F = learningCurves(dat_switches,plt_switches,p);

% save plot
savefig(F,[save_root_fig,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r0.fig']);
saveas(F,[save_root_png,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r0.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r0.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r0.txt'],'wt');
fprintf(fid,['\nBehaviouralEffectAnalysis_LearningCurve_p1e01r0\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']); fclose(fid);
fprintf(['\nBehaviouralEffectAnalysis_LearningCurve_p1e01r0\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']);


%% --- running and e01 drill down ---

%% Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1

% get data
dat = {};
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==7),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct,{d{[intersect(find(d_info.group==8),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1)))],j}}','Uniform',0))*100;
    end
end
dat_switches = {}; 
dat_switches{1,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
dat_switches{2,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq

% make plot
plt_switches = struct(); plt_switches.ylim = [35,100]; plt_switches.flag = 'paper'; plt_switches.traces = 'individuals and mean'; plt_switches.xlabel = 'Trial block'; plt_switches.col = {};
plt_switches.col{1,1} = p.col.ctrl; plt_switches.col{1,2} = p.col.ctrl;
plt_switches.col{2,1} = p.col.seq; plt_switches.col{2,2} = p.col.seq;
F = learningCurves(dat_switches,plt_switches,p);

% save plot
savefig(F,[save_root_fig,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1.fig']);
saveas(F,[save_root_png,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1.txt'],'wt');
fprintf(fid,['\nBehaviouralEffectAnalysis_LearningCurve_p1e01r1\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']); fclose(fid);
fprintf(['\nBehaviouralEffectAnalysis_LearningCurve_p1e01r1\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']);


%% Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_stim

% get data
dat = {};
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==7),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==8),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1)))],j}}','Uniform',0))*100;
    end
end
dat_switches = {}; 
dat_switches{1,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
dat_switches{2,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq

% make plot
plt_switches = struct(); plt_switches.ylim = [35,100]; plt_switches.flag = 'paper'; plt_switches.traces = 'individuals and mean'; plt_switches.xlabel = 'Trial block'; plt_switches.col = {};
plt_switches.col{1,1} = p.col.ctrl; plt_switches.col{1,2} = p.col.ctrl;
plt_switches.col{2,1} = p.col.seq; plt_switches.col{2,2} = p.col.seq;
F = learningCurves(dat_switches,plt_switches,p);

% save plot
savefig(F,[save_root_fig,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_stim.fig']);
saveas(F,[save_root_png,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_stim.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_stim.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_stim.txt'],'wt');
fprintf(fid,['\nBehaviouralEffectAnalysis_LearningCurve_p1e01r1_stim\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']); fclose(fid);
fprintf(['\nBehaviouralEffectAnalysis_LearningCurve_p1e01r1_stim\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']);


%% Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_catch

% get data
dat = {};
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==7),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==8),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1)))],j}}','Uniform',0))*100;
    end
end
dat_switches = {}; 
dat_switches{1,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
dat_switches{2,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq

% make plot
plt_switches = struct(); plt_switches.ylim = [35,100]; plt_switches.flag = 'paper'; plt_switches.traces = 'individuals and mean'; plt_switches.xlabel = 'Trial block'; plt_switches.col = {};
plt_switches.col{1,1} = p.col.ctrl; plt_switches.col{1,2} = p.col.ctrl;
plt_switches.col{2,1} = p.col.seq; plt_switches.col{2,2} = p.col.seq;
F = learningCurves(dat_switches,plt_switches,p);

% save plot
savefig(F,[save_root_fig,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_catch.fig']);
saveas(F,[save_root_png,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_catch.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_catch.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_catch.txt'],'wt');
fprintf(fid,['\nBehaviouralEffectAnalysis_LearningCurve_p1e01r1_catch\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']); fclose(fid);
fprintf(['\nBehaviouralEffectAnalysis_LearningCurve_p1e01r1_catch\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']);


%% Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_diff

% get data
dat = {};
for j=1:d_info.numDays
    if j==1
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==7)],j}}','Uniform',0))*100 - cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==8)],j}}','Uniform',0))*100 - cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
    else
        dat{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==7),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1)))],j}}','Uniform',0))*100 - cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==7),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1)))],j}}','Uniform',0))*100;
        dat{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==8),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1)))],j}}','Uniform',0))*100 - cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==8),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1)))],j}}','Uniform',0))*100;
    end
end
dat_switches = {}; 
dat_switches{1,1} = [dat{2,2};dat{1,4}]; dat_switches{2,2} = [dat{2,3};dat{1,5}]; % ctrl
dat_switches{2,1} = [dat{1,2};dat{2,4}]; dat_switches{1,2} = [dat{1,3};dat{2,5}]; % seq

% make plot
plt_switches = struct(); plt_switches.ylim = [-50,50]; plt_switches.flag = 'paper'; plt_switches.traces = 'individuals and mean'; plt_switches.xlabel = 'Trial block'; plt_switches.col = {};
plt_switches.col{1,1} = p.col.ctrl; plt_switches.col{1,2} = p.col.ctrl;
plt_switches.col{2,1} = p.col.seq; plt_switches.col{2,2} = p.col.seq;
F = learningCurves(dat_switches,plt_switches,p);

% save plot
savefig(F,[save_root_fig,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_diff.fig']);
saveas(F,[save_root_png,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_diff.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_diff.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_diff.txt'],'wt');
fprintf(fid,['\nBehaviouralEffectAnalysis_LearningCurve_p1e01r1_diff\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']); fclose(fid);
fprintf(['\nBehaviouralEffectAnalysis_LearningCurve_p1e01r1_diff\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']);


%% --- In progress ---


%% Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_stim_vs_catch

% % get data
% dat_stim = {};
% for j=1:d_info.numDays
%     if j==1
%         dat_stim{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
%         dat_stim{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
%     else
%         dat_stim{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==7),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1)))],j}}','Uniform',0))*100;
%         dat_stim{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_stim,{d{[intersect(find(d_info.group==8),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1)))],j}}','Uniform',0))*100;
%     end
% end
% dat_switches_stim = {}; 
% dat_switches_stim{1,1} = [dat_stim{2,2};dat_stim{1,4}]; dat_switches_stim{2,2} = [dat_stim{2,3};dat_stim{1,5}]; % ctrl
% dat_switches_stim{2,1} = [dat_stim{1,2};dat_stim{2,4}]; dat_switches_stim{1,2} = [dat_stim{1,3};dat_stim{2,5}]; % seq
% dat_catch = {};
% for j=1:d_info.numDays
%     if j==1
%         dat_catch{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==7)],j}}','Uniform',0))*100;
%         dat_catch{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[find(d_info.group==8)],j}}','Uniform',0))*100;
%     else
%         dat_catch{1,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==7),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1)))],j}}','Uniform',0))*100;
%         dat_catch{2,j} = cell2mat(cellfun(@(c) c.perf.('blocks_general').correct_catch,{d{[intersect(find(d_info.group==8),intersect(find(d_info.presponsive(:,j)==1),find(d_info.running(:,j)==1)))],j}}','Uniform',0))*100;
%     end
% end
% dat_switches_catch = {}; 
% dat_switches_catch{1,1} = [dat_catch{2,2};dat_catch{1,4}]; dat_switches_catch{2,2} = [dat_catch{2,3};dat_catch{1,5}]; % ctrl
% dat_switches_catch{2,1} = [dat_catch{1,2};dat_catch{2,4}]; dat_switches_catch{1,2} = [dat_catch{1,3};dat_catch{2,5}]; % seq
% 
% % make plot
% plt_switches = struct(); plt_switches.ylim = [35,100]; plt_switches.flag = 'paper'; plt_switches.traces = 'individuals and mean'; plt_switches.xlabel = 'Trial block'; plt_switches.col = {};
% plt_switches.col{1,1} = p.col.ctrl; plt_switches.col{1,2} = p.col.ctrl;
% plt_switches.col{2,1} = p.col.seq; plt_switches.col{2,2} = p.col.seq;
% F = learningCurves(dat_switches_stim,plt_switches,p);
% 
% % make plot
% plt_switches = struct(); plt_switches.ylim = [35,100]; plt_switches.flag = 'paper'; plt_switches.traces = 'individuals and mean'; plt_switches.xlabel = 'Trial block'; plt_switches.col = {};
% plt_switches.col{1,1} = p.col.ctrl; plt_switches.col{1,2} = p.col.ctrl;
% plt_switches.col{2,1} = p.col.seq; plt_switches.col{2,2} = p.col.seq;
% F = learningCurves(dat_switches,plt_switches,p);
% 
% % save plot
% savefig(F,[save_root_fig,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_stimvscatch.fig']);
% saveas(F,[save_root_png,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_stimvscatch.png']);
% set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_stimvscatch.pdf']); set(gcf,'Color',[1,1,1])
% fid = fopen([save_root_txt,'\Fig5_BehaviouralEffectAnalysis_LearningCurve_p1e01r1_stimvscatch.txt'],'wt');
% fprintf(fid,['\nBehaviouralEffectAnalysis_LearningCurve_p1e01r1_stimvscatch\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']); fclose(fid);
% fprintf(['\nBehaviouralEffectAnalysis_LearningCurve_p1e01r1_stimvscatch\nn(d1_ctrl)=',num2str(size(dat_switches{1,1},1)),'\nn(d2_ctrl)=',num2str(size(dat_switches{1,2},1)),'\nn(d1_seq)=',num2str(size(dat_switches{2,1},1)),'\nn(d2_seq)=',num2str(size(dat_switches{2,2},1)),'\n']);
