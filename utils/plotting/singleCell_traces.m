%function singleCell_traces(sca,idx)
idx = idcs_unsorted_Aseq(50);
%sca.prop.metric = 'dFF';

default_figure();

for i=1:2
    if i==1
        data = squeeze(sca.traces_A(idx,:,:))';
    elseif i==2
        data = squeeze(sca.traces_X(idx,:,:))';
    end
    
    subplot(1,2,i)
    
    for j=1:size(data,1)
        plot(sca.prop.t_binned,data(j,:))
        hold on
    end
    
    hold off
    
    xlabel('Time (s)')
    ylabel('Trial')  
end

% colorbar(winter)
% ylabel('Trial')