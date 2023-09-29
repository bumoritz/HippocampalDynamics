function taskLines(p,info,timeWindow,plotType,do_labels,do_x,do_y,this_t,paper)

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
if nargin < 8
    this_t=p.general.t_binned; % p.general.t_unbinned
end
if nargin < 9
    paper = false;
end

if strcmp(plotType,'traces')
    offset = 0;
elseif strcmp(plotType,'heatmap')
    offset = -0.5;
elseif strcmp(plotType,'decoding')
    offset = -15.5;
end

if paper
%     patch([offset+interp1(this_t,1:length(this_t),0),offset+interp1(this_t,1:length(this_t),0)],[ylim],'Color',p.col.odour)
% xline(offset+interp1(this_t,1:length(this_t),0),'Color',p.col.odour,'LineWidth',3);
% xline(offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1),'Color',p.col.odour,'LineWidth',3);
% xline(offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),'Color',p.col.odour,'LineWidth',3);
% xline(offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2),'Color',p.col.odour,'LineWidth',3);
% xline(offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),'Color',p.col.reward,'LineWidth',3);
% xline(offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow),'Color',p.col.reward,'LineWidth',3);        
else
    if do_x
        if strcmp(timeWindow,'full')
                ymin = -1000; ymax = 1000;
                patch([interp1(this_t,1:length(this_t),0),interp1(this_t,1:length(this_t),0),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1)],...
                    [ymin,ymax,ymax,ymin],mean([p.col.odour;p.col.white]),'EdgeColor','none');
                patch([interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2)],...
                    [ymin,ymax,ymax,ymin],mean([p.col.odour;p.col.white]),'EdgeColor','none');
                patch([interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow),interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow)],...
                    [ymin,ymax,ymax,ymin],mean([p.col.reward;p.col.white]),'EdgeColor','none');
%             xline(offset+interp1(this_t,1:length(this_t),0),'Color',p.col.odour,'LineWidth',3);
%             xline(offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1),'Color',p.col.odour,'LineWidth',3);
%             xline(offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),'Color',p.col.odour,'LineWidth',3);
%             xline(offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2),'Color',p.col.odour,'LineWidth',3);
%             xline(offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),'Color',p.col.reward,'LineWidth',3);
%             xline(offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow),'Color',p.col.reward,'LineWidth',3);
        elseif strcmp(timeWindow,'AW')
            xline(offset,'Color',p.col.odour,'LineWidth',3);
            xline(-interp1(this_t,1:length(this_t),0)+offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1),'Color',p.col.odour,'LineWidth',3);
            xline(-interp1(this_t,1:length(this_t),0)+offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),'Color',p.col.odour,'LineWidth',3);
            xlim([0,-interp1(this_t,1:length(this_t),0)+offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap)])
        end
        if do_labels
            xlabel('Time (s)')
            if strcmp(plotType,'decoding')
                lower = 0.5;
                upper = 52.5;
                xlim([lower,upper])
                set(gca,'xtick',-0.5+interp1(this_t(16:67),1:52,ceil(this_t(16)):floor(this_t(67))));
                set(gca,'xticklabel',[ceil(this_t(16)):floor(this_t(67))]);
            else
                lower = 0;
                upper = length(this_t);
                xlim([lower,upper])
                set(gca,'xtick',interp1(this_t,1:length(this_t),ceil(this_t(1)):floor(this_t(end))));
                set(gca,'xticklabel',[ceil(this_t(1)):floor(this_t(end))]);
            end
        end
    end

    if do_y
        if strcmp(timeWindow,'full')
            yline(offset+interp1(this_t,1:length(this_t),0),'Color',p.col.odour,'LineWidth',3);
            yline(offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1),'Color',p.col.odour,'LineWidth',3);
            yline(offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),'Color',p.col.odour,'LineWidth',3);
            yline(offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2),'Color',p.col.odour,'LineWidth',3);
            yline(offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay),'Color',p.col.reward,'LineWidth',3);
            yline(offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap+info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow),'Color',p.col.reward,'LineWidth',3);
        elseif strcmp(timeWindow,'AW')
            yline(offset,'Color',p.col.odour,'LineWidth',3);
            yline(-interp1(this_t,1:length(this_t),0)+offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1),'Color',p.col.odour,'LineWidth',3);
            yline(-interp1(this_t,1:length(this_t),0)+offset+interp1(this_t,1:length(this_t),info.task.trialStructure.tOdour1+info.task.trialStructure.tGap),'Color',p.col.odour,'LineWidth',3);
        end
    end
end
end