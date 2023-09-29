function heatMap_task(data,these_lines,sca,p,info,plt)

if ~exist('plt','var')
    plt = struct();
end
if ~isfield(plt,'clim')
    plt.clim = [nanmin(data(:)),nanmax(data(:))];
end
if ~isfield(plt,'colormap')
    plt.colormap = parula;
end
if ~isfield(plt,'illustrator')
    plt.illustrator = false;
end
if ~isfield(plt,'lineWidth')
    plt.lineWidth = 3;
end


imagesc(data);
xlim([0,size(data,2)]+0.5);
ylim([0,size(data,1)]+0.5);

set(gca,'CLim',plt.clim);
colormap(gca,plt.colormap);

% task events
hold on
if ~isnan(these_lines)
    if length(these_lines)==1
        yline(these_lines+0.5,'Color',p.col.white,'LineWidth',plt.lineWidth);
    elseif length(these_lines)>1
        for i=1:length(these_lines)
            yline(these_lines(i)+0.5,'Color',p.col.black,'LineWidth',plt.lineWidth);
        end
    end
end
try
    xline(-0.5+interp1(p.general.t_binned,1:length(p.general.t_binned),0),'Color',p.col.odour,'LineWidth',plt.lineWidth,'alpha',1);
    %xline(-0.5+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1),'Color',p.col.odour,'LineWidth',plt.lineWidth);
    xline(-0.5+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),'Color',p.col.odour,'LineWidth',plt.lineWidth,'alpha',1);
    %xline(-0.5+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2),'Color',p.col.odour,'LineWidth',plt.lineWidth);
    xline(-0.5+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),'Color',p.col.reward,'LineWidth',plt.lineWidth,'alpha',1);
    %xline(-0.5+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow),'Color',p.col.reward,'LineWidth',plt.lineWidth);
catch
    xline(-0.5+interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),0),'Color',p.col.odour,'LineWidth',plt.lineWidth,'alpha',1);
    %xline(-0.5+interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1),'Color',p.col.odour,'LineWidth',plt.lineWidth);
    xline(-0.5+interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),'Color',p.col.odour,'LineWidth',plt.lineWidth,'alpha',1);
    %xline(-0.5+interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2),'Color',p.col.odour,'LineWidth',plt.lineWidth);
    xline(-0.5+interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),'Color',p.col.reward,'LineWidth',plt.lineWidth,'alpha',1);
    %xline(-0.5+interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow),'Color',p.col.reward,'LineWidth',plt.lineWidth);
end
hold off

if ~plt.illustrator
    try        
        set(gca,'xtick',-0.5+interp1(p.general.t_binned,1:length(p.general.t_binned),ceil(p.general.t_binned(1)):floor(p.general.t_binned(end))));
        set(gca,'xticklabel',[ceil(p.general.t_binned(1)):floor(p.general.t_binned(end))]);        
    catch
        set(gca,'xtick',-0.5+interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),ceil(sca.prop.t_binned(1)):floor(sca.prop.t_binned(end))));
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

if isfield(plt,'yticks') && isfield(plt,'yticklabels')
    set(gca,'box','on');
    try
        yticks(plt.yticks)
        yticklabels(plt.yticklabels)
    catch
    end
    set(gca,'xcolor',p.col.black);
    set(gca,'ycolor',p.col.black);
end

if isfield(plt,'ylim')
    ylim(plt.ylim);
end
if isfield(plt,'ylabel')
    ylabel(plt.ylabel);
end

end
