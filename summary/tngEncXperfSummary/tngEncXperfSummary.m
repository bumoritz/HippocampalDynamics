function [tngEncXperf] = tngEncXperfSummary(d_info,d,ops,p,path)

if ~exist([path.root_summary,'plots\tngEncXperf'],'dir')
    mkdir([path.root_summary,'plots\tngEncXperf']);
end


%% Extract data

% perf.general -> performance, H, CR
performance = extractVariable(d,'perf.general.correct','array','single');
H = extractVariable(d,'perf.general.H','array','single');
CR = extractVariable(d,'perf.general.CR','array','single');

% tng_all.passed_stats.AW -> fractionOfCells
fractionOfCells.A       = extractVariable(d,'tng_all.passed_stats.AW.A_fractionOfCells','array','single');
fractionOfCells.X       = extractVariable(d,'tng_all.passed_stats.AW.X_fractionOfCells','array','single');
fractionOfCells.Aonly   = extractVariable(d,'tng_all.passed_stats.AW.Aonly_fractionOfCells','array','single');
fractionOfCells.Xonly   = extractVariable(d,'tng_all.passed_stats.AW.Xonly_fractionOfCells','array','single');
fractionOfCells.AandX   = extractVariable(d,'tng_all.passed_stats.AW.AandX_fractionOfCells','array','single');
fractionOfCells.AorX    = extractVariable(d,'tng_all.passed_stats.AW.AorX_fractionOfCells','array','single');
fractionOfCells.notAorX = 1-extractVariable(d,'tng_all.passed_stats.AW.AorX_fractionOfCells','array','single');

% build binary running vector (temporary)
running = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if isfield(d{i,1},'paq_beh') && (~isnan(nanmean(d{i,1}.paq_beh.speed)))
        running(i) = mean(d{i,1}.paq_beh.speed) > 20;
    end
end
running(16) = true; % Arwen
running(27) = false; % Stanage
running(36:42) = 0 % Python, correct after data is imported


%% Figure: running -> fractionOfCells

j = 1;

nrows = 2; ncols = 4;
F = default_figure([20,0.5,20,9.9]);

these_fields = fields(fractionOfCells);
for i=1:length(these_fields)
    this_data_running = fractionOfCells.(these_fields{i})(find(running==1),j);
    this_data_notrunning = fractionOfCells.(these_fields{i})(find(running==0),j);
end




% plt = struct(); plt.xlabel = 'Number of responsive targets'; plt.xlim = [0,240]; plt.xticks = [0,40,80,120,160,200,240];
% F = correlationPlot(responsiveTargets,performance,p,plt);
% savefig(F,[path.root_summary,'plots\respXperf\responsiveTargets_performance.fig']);
% saveas(F,[path.root_summary,'plots\respXperf\responsiveTargets_performance.png']);












