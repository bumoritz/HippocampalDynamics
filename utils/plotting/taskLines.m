function taskLines(p,info,timeWindow,plotType,do_labels,do_x,do_y)


if nargin < 3
    timeWindow='full'; % 'full','AW'
end
if nargin < 4
    plotType='traces'; % 'traces','heatmap','decoding'
end
if nargin < 5
    do_labels=true;
end
if nargin < 6
    do_x=1;
end
if nargin < 7
    do_y=0;
end

if strcmp(plotType,'traces')
    offset = 0;
elseif strcmp(plotType,'heatmap')
    offset = -0.5;
elseif strcmp(plotType,'decoding')
    offset = -15.5;
end

if do_x
    if strcmp(timeWindow,'full')
        xline(offset+interp1(p.general.t_binned,1:length(p.general.t_binned),0),'Color',p.col.odour,'LineWidth',3);
        xline(offset+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1),'Color',p.col.odour,'LineWidth',3);
        xline(offset+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),'Color',p.col.odour,'LineWidth',3);
        xline(offset+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2),'Color',p.col.odour,'LineWidth',3);
        xline(offset+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),'Color',p.col.reward,'LineWidth',3);
        xline(offset+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow),'Color',p.col.reward,'LineWidth',3);
    elseif strcmp(timeWindow,'AW')
        xline(offset,'Color',p.col.odour,'LineWidth',3);
        xline(-interp1(p.general.t_binned,1:length(p.general.t_binned),0)+offset+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1),'Color',p.col.odour,'LineWidth',3);
        xline(-interp1(p.general.t_binned,1:length(p.general.t_binned),0)+offset+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),'Color',p.col.odour,'LineWidth',3);
        xlim([0,-interp1(p.general.t_binned,1:length(p.general.t_binned),0)+offset+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap)])
    end
    if do_labels
        xlabel('Time (s)')
        if strcmp(plotType,'decoding')
            lower = 0.5;
            upper = 52.5;
            xlim([lower,upper])
            set(gca,'xtick',-0.5+interp1(p.general.t_binned(16:67),1:52,ceil(p.general.t_binned(16)):floor(p.general.t_binned(67))));
            set(gca,'xticklabel',[ceil(p.general.t_binned(16)):floor(p.general.t_binned(67))]);
        else
            lower = 0;
            upper = length(p.general.t_binned);
            xlim([lower,upper])
            set(gca,'xtick',interp1(p.general.t_binned,1:length(p.general.t_binned),ceil(p.general.t_binned(1)):floor(p.general.t_binned(end))));
            set(gca,'xticklabel',[ceil(p.general.t_binned(1)):floor(p.general.t_binned(end))]);
        end
    end
end

if do_y
    if strcmp(timeWindow,'full')
        yline(offset+interp1(p.general.t_binned,1:length(p.general.t_binned),0),'Color',p.col.odour,'LineWidth',3);
        yline(offset+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1),'Color',p.col.odour,'LineWidth',3);
        yline(offset+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),'Color',p.col.odour,'LineWidth',3);
        yline(offset+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2),'Color',p.col.odour,'LineWidth',3);
        yline(offset+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),'Color',p.col.reward,'LineWidth',3);
        yline(offset+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow),'Color',p.col.reward,'LineWidth',3);
    elseif strcmp(timeWindow,'AW')
        yline(offset,'Color',p.col.odour,'LineWidth',3);
        yline(-interp1(p.general.t_binned,1:length(p.general.t_binned),0)+offset+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1),'Color',p.col.odour,'LineWidth',3);
        yline(-interp1(p.general.t_binned,1:length(p.general.t_binned),0)+offset+interp1(p.general.t_binned,1:length(p.general.t_binned),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),'Color',p.col.odour,'LineWidth',3);
    end
end

end