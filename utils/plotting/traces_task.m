function traces_task(sca,p,info,plt,paper_version)

%% Preparations

if ~exist('plt','var')
    plt = struct();
end
if ~isfield(plt,'illustrator')
    plt.illustrator = false;
end
if nargin<5
    paper_version = false;
end


%% Core

% task events
hold on
if paper_version
    this_t = p.general.t_binned; ymin = -0.1; ymax = 1.2; %ymin = -0.1; ymax = 1.2; 
    patch([interp1(this_t,1:length(this_t),0),interp1(this_t,1:length(this_t),0),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1)],...
        [ymin,ymax,ymax,ymin],mean([p.col.odour;p.col.white]),'EdgeColor','none');
    patch([interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2)],...
        [ymin,ymax,ymax,ymin],mean([p.col.odour;p.col.white]),'EdgeColor','none');
    patch([interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow)],...
        [ymin,ymax,ymax,ymin],mean([p.col.reward;p.col.white]),'EdgeColor','none');
else
    try
        xline(interp1(p.general.t_binned,1:length(p.general.t_binned),0),'Color',p.col.odour,'LineWidth',3);
        xline(interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1),'Color',p.col.odour,'LineWidth',3);
        xline(interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),'Color',p.col.odour,'LineWidth',3);
        xline(interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2),'Color',p.col.odour,'LineWidth',3);
        xline(interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),'Color',p.col.reward,'LineWidth',3);
        xline(interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow),'Color',p.col.reward,'LineWidth',3);
    catch
        xline(interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),0),'Color',p.col.odour,'LineWidth',3);
        xline(interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1),'Color',p.col.odour,'LineWidth',3);
        xline(interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),'Color',p.col.odour,'LineWidth',3);
        xline(interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2),'Color',p.col.odour,'LineWidth',3);
        xline(interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),'Color',p.col.reward,'LineWidth',3);
        xline(interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow),'Color',p.col.reward,'LineWidth',3);
    end
end
hold off

try
	xlim([0,length(p.general.t_binned)])
catch
    xlim([0,length(sca.prop.t_binned)])
end
if isfield(plt,'ylim')
    ylim(plt.ylim)
end

if ~plt.illustrator
    try
        set(gca,'xtick',interp1(p.general.t_binned,1:length(p.general.t_binned),ceil(p.general.t_binned(1)):floor(p.general.t_binned(end))));
        set(gca,'xticklabel',[ceil(p.general.t_binned(1)):floor(p.general.t_binned(end))]);
    catch
        set(gca,'xtick',interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),ceil(sca.prop.t_binned(1)):floor(sca.prop.t_binned(end))));
        set(gca,'xticklabel',[ceil(sca.prop.t_binned(1)):floor(sca.prop.t_binned(end))]);
    end
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
