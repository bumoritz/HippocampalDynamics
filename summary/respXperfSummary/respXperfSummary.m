function respXperfSummary(d_info,d,ops,p,path)

if ~exist([path.root_summary,'plots\respXperf'],'dir')
    mkdir([path.root_summary,'plots\respXperf']);
end

%% Extract data

% perf.general.correct -> performance
temp = extractVariable(d,'perf.general.correct','array','single');
performance = createStimSquare(temp,d_info,'seq-ctrl_switch-postswitch');

% perf.general.H -> H
temp = extractVariable(d,'perf.general.H','array','single');
H = createStimSquare(temp,d_info,'seq-ctrl_switch-postswitch');

% perf.general.CR -> CR
temp = extractVariable(d,'perf.general.CR','array','single');
CR = createStimSquare(temp,d_info,'seq-ctrl_switch-postswitch');

% resp.responders_0d5z.numRespTargeted -> responsiveTargets
temp = extractVariable(d,'resp.responders_0d5z.numRespTargeted','array','single');
responsiveTargets = createStimSquare(temp,d_info,'seq-ctrl_switch-postswitch');

% resp.responders_0d5z.proportion_RespNonTargeted_RespTargeted -> offTargetsPerTarget
temp = extractVariable(d,'resp.responders_0d5z.proportion_RespNonTargeted_RespTargeted','array','single');
offTargetsPerTarget = createStimSquare(temp,d_info,'seq-ctrl_switch-postswitch');

% resp.responders_0d5z.specificity_respTargeted -> specificity
temp = extractVariable(d,'resp.responders_0d5z.specificity_respTargeted','array','single');
specificity = createStimSquare(temp,d_info,'seq-ctrl_switch-postswitch');


%% Make correlation plots - subplots

% performance ~ number of responsive targets
plt = struct(); plt.xlabel = 'Number of responsive targets'; plt.xlim = [0,240]; plt.xticks = [0,40,80,120,160,200,240];
F = correlationPlot(responsiveTargets,performance,p,plt);
savefig(F,[path.root_summary,'plots\respXperf\responsiveTargets_performance.fig']);
saveas(F,[path.root_summary,'plots\respXperf\responsiveTargets_performance.png']);

% performance ~ responsive off-targets per responsive target
plt = struct(); plt.xlabel = 'Responsive off-targets per responsive target'; plt.xlim_method = 'zeromax';
F = correlationPlot(offTargetsPerTarget,performance,p,plt);
savefig(F,[path.root_summary,'plots\respXperf\offTargetsPerTarget_performance.fig']);
saveas(F,[path.root_summary,'plots\respXperf\offTargetsPerTarget_performance.png']);

% performance ~ response specificity
plt = struct(); plt.xlabel = 'Response specificity'; plt.xlim = [0,1];
F = correlationPlot(specificity,performance,p,plt);
savefig(F,[path.root_summary,'plots\respXperf\specificity_performance.fig']);
saveas(F,[path.root_summary,'plots\respXperf\specificity_performance.png']);

% H ~ number of responsive targets
plt = struct(); plt.xlabel = 'Number of responsive targets'; plt.xlim = [0,240]; plt.xticks = [0,40,80,120,160,200,240]; plt.ylabel = 'Hit rate'; 
F = correlationPlot(responsiveTargets,H,p,plt);
savefig(F,[path.root_summary,'plots\respXperf\responsiveTargets_H.fig']);
saveas(F,[path.root_summary,'plots\respXperf\responsiveTargets_H.png']);

% CR ~ number of responsive targets
plt = struct(); plt.xlabel = 'Number of responsive targets'; plt.xlim = [0,240]; plt.xticks = [0,40,80,120,160,200,240]; plt.ylabel = 'Correct rejection rate'; 
F = correlationPlot(responsiveTargets,CR,p,plt);
savefig(F,[path.root_summary,'plots\respXperf\responsiveTargets_CR.fig']);
saveas(F,[path.root_summary,'plots\respXperf\responsiveTargets_CR.png']);


%% Make correlation plots - combined

% performance ~ number of responsive targets
plt = struct(); plt.xlabel = 'Number of responsive targets'; plt.xlim = [0,240]; plt.xticks = [0,40,80,120,160,200,240]; plt.type = 'combined';
F = correlationPlot(responsiveTargets,performance,p,plt);
savefig(F,[path.root_summary,'plots\respXperf\responsiveTargets_performance_combined.fig']);
saveas(F,[path.root_summary,'plots\respXperf\responsiveTargets_performance_combined.png']);

% performance ~ responsive off-targets per responsive target
plt = struct(); plt.xlabel = 'Responsive off-targets per responsive target'; plt.xlim_method = 'zeromax'; plt.type = 'combined';
F = correlationPlot(offTargetsPerTarget,performance,p,plt);
savefig(F,[path.root_summary,'plots\respXperf\offTargetsPerTarget_performance_combined.fig']);
saveas(F,[path.root_summary,'plots\respXperf\offTargetsPerTarget_performance_combined.png']);

% performance ~ response specificity
plt = struct(); plt.xlabel = 'Response specificity'; plt.xlim = [0,1]; plt.type = 'combined';
F = correlationPlot(specificity,performance,p,plt);
savefig(F,[path.root_summary,'plots\respXperf\specificity_performance_combined.fig']);
saveas(F,[path.root_summary,'plots\respXperf\specificity_performance_combined.png']);

% H ~ number of responsive targets
plt = struct(); plt.xlabel = 'Number of responsive targets'; plt.xlim = [0,240]; plt.xticks = [0,40,80,120,160,200,240]; plt.ylabel = 'Hit rate'; plt.type = 'combined';
F = correlationPlot(responsiveTargets,H,p,plt);
savefig(F,[path.root_summary,'plots\respXperf\responsiveTargets_H_combined.fig']);
saveas(F,[path.root_summary,'plots\respXperf\responsiveTargets_H_combined.png']);

% CR ~ number of responsive targets
plt = struct(); plt.xlabel = 'Number of responsive targets'; plt.xlim = [0,240]; plt.xticks = [0,40,80,120,160,200,240]; plt.ylabel = 'Correct rejection rate'; plt.type = 'combined';
F = correlationPlot(responsiveTargets,CR,p,plt);
savefig(F,[path.root_summary,'plots\respXperf\responsiveTargets_CR_combined.fig']);
saveas(F,[path.root_summary,'plots\respXperf\responsiveTargets_CR_combined.png']);


%% Return

disp(['--- Saved respXperf figures to ',[path.root_summary,'plots\respXperf'],'.'])
if ops.close_figures
    close all;
end
end