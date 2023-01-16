function traces_task(sca,p,info,plt)

%% Preparations

if ~exist('plt','var')
    plt = struct();
end
if ~isfield(plt,'illustrator')
    plt.illustrator = false;
end


%% Core

% task events
hold on
xline(interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),0),'Color',p.col.odour,'LineWidth',3);
xline(interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1),'Color',p.col.odour,'LineWidth',3);
xline(interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),'Color',p.col.odour,'LineWidth',3);
xline(interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2),'Color',p.col.odour,'LineWidth',3);
xline(interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),'Color',p.col.reward,'LineWidth',3);
xline(interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow),'Color',p.col.reward,'LineWidth',3);
hold off

xlim([0,length(sca.prop.t_binned)])
if isfield(plt,'ylim')
    ylim(plt.ylim)
end

if ~plt.illustrator
    set(gca,'xtick',interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),ceil(sca.prop.t_binned(1)):floor(sca.prop.t_binned(end))));
    set(gca,'xticklabel',[ceil(sca.prop.t_binned(1)):floor(sca.prop.t_binned(end))]);
    xlabel('Time (s)')
else
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    set(gca,'box','off');
    set(gca,'xcolor','none');
    set(gca,'ycolor','none');
end

end
