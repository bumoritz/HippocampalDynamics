%function singleCell_waterfall(sca,idx)
idx = idcs_sorted_Aseq(45);
sca.prop.metric = 'dFF';

default_figure();

subplot(1,2,1)
data_A = squeeze(sca.traces_A(idx,:,:))';
waterfall(data_A);
xlabel('Time')
ylabel('Trial')
if strcmp(sca.prop.metric,'dFF')
    zlabel('dFF')
end

subplot(1,2,2)
data_X = squeeze(sca.traces_X(idx,:,:))';
waterfall(data_X);
xlabel('Time')
ylabel('Trial')
if strcmp(sca.prop.metric,'dFF')
    zlabel('dFF')
end


% xlabels = -1:1:(time(end)-1);
% xticks = floor((xlabels+1)/BinDurationAfterBinning)+2;
% set(gca, 'XTick', xticks, 'XTickLabel', xlabels);
% xlabel('Time (s)')
% ylabel('Trial Number')
% line([first_odour first_odour],[0 size(data,1)],'Color','red','LineWidth',2);
% line([second_odour second_odour],[0 size(data,1)],'Color','red','LineWidth',2);
% line([first_odour+round(1/BinDurationAfterBinning) first_odour+round(1/BinDurationAfterBinning)],[0 size(data,1)],'Color','red','LineWidth',2,'LineStyle','--');
% line([second_odour+round(1/BinDurationAfterBinning) second_odour+round(1/BinDurationAfterBinning)],[0 size(data,1)],'Color','red','LineWidth',2,'LineStyle','--');
% title([Relative ' trials'])