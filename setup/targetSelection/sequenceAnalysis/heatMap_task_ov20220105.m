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

imagesc(data);

% 0.0756
% %nanmax(nanmean(data,1))
% set(gca,'CLim',[0,0.0756]); %set(gca,'CLim',[0,nanmax(data(:))/2]); %

set(gca,'CLim',plt.clim);
colormap(gca,plt.colormap);

% task events
hold on
if ~isnan(these_lines)
    if length(these_lines)==1
        yline(these_lines+0.5,'Color',p.col.white,'LineWidth',3);
    elseif length(these_lines)>1
        for i=1:length(these_lines)
            yline(these_lines(i)+0.5,'Color',p.col.black,'LineWidth',3);
        end
    end
end
xline(-0.5+interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),0),'Color',p.col.odour,'LineWidth',3);
xline(-0.5+interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1),'Color',p.col.odour,'LineWidth',3);
xline(-0.5+interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),'Color',p.col.odour,'LineWidth',3);
xline(-0.5+interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2),'Color',p.col.odour,'LineWidth',3);
xline(-0.5+interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),'Color',p.col.reward,'LineWidth',3);
xline(-0.5+interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow),'Color',p.col.reward,'LineWidth',3);
hold off

if ~plt.illustrator
    set(gca,'xtick',-0.5+interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),ceil(sca.prop.t_binned(1)):floor(sca.prop.t_binned(end))));
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

if isfield(plt,'yticks') && isfield(plt,'yticklabels')
    set(gca,'box','on');
    yticks(plt.yticks)
    yticklabels(plt.yticklabels)
    set(gca,'xcolor',p.col.black);
    set(gca,'ycolor',p.col.black);
end
end
