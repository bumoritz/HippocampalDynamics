%% Plot sequence cells by block, bin and animal - all 5 days
% 
% default_figure(); 
% 
% subplot(1,2,1); hold on;
% this_data = seqCellsByTime_num;
% for n=1:numBins
%     temp0 = flipud(eval('winter'));
%     [~,temp1] = discretize([1,numBins],size(temp0,1));
%     temp2 = discretize(1:numBins,temp1);
%     shadedErrorBar(1:numBlocks,nanmean(this_data(n,:,:),3),nansem(this_data(n,:,:),3),'lineProps',temp0(temp2(n),:));
% end
% temp = cumsum(blocks(1:end-1))+0.5;
% for i=1:length(temp)
%     if i==1 || i==3
%         xline(temp(i),'k-');
%     else
%         xline(temp(i),'k:');
%     end
% end
% xlim([0,numBlocks])
% xlabel('trial block')
% ylabel('Number of sequence cells')
% 
% subplot(1,2,2); hold on;
% this_data = seqCellsByTime_prop;
% for n=1:numBins
%     temp0 = flipud(eval('winter'));
%     [~,temp1] = discretize([1,numBins],size(temp0,1));
%     temp2 = discretize(1:numBins,temp1);
%     shadedErrorBar(1:numBlocks,nanmean(this_data(n,:,:),3)*100,nansem(this_data(n,:,:),3)*100,'lineProps',temp0(temp2(n),:));
% end
% temp = cumsum(blocks(1:end-1))+0.5;
% for i=1:length(temp)
%     if i==1 || i==3
%         xline(temp(i),'k-');
%     else
%         xline(temp(i),'k:');
%     end
% end
% xlim([0,numBlocks])
% xlabel('Trial block')
% ytickformat('percentage')
% ylabel('Proportion of sequence cells (rel. to all cells)')
% 


%% Plot correlation with performance - animal-by-animal

% default_figure();
% 
% idx = 3;
% 
% subplot(1,2,1); hold on;
% this_data = seqCellsByTime_num;
% for n=1:numBins
%     temp0 = flipud(eval('winter'));
%     [~,temp1] = discretize([1,numBins],size(temp0,1));
%     temp2 = discretize(1:numBins,temp1);
%     plot(1:numBlocks,this_data(n,:,idx),'Color',temp0(temp2(n),:))
% end
% temp = cumsum(blocks(1:end-1))+0.5;
% for i=1:length(temp)
%     if i==1 || i==3
%         xline(temp(i),'k-');
%     else
%         xline(temp(i),'k:');
%     end
% end
% xlim([0,numBlocks])
% xlabel('trial block')
% ylabel('Number of sequence cells')
% 
% subplot(1,2,2); hold on;
% this_data = seqCellsByTime_prop;
% for n=1:numBins
%     temp0 = flipud(eval('winter'));
%     [~,temp1] = discretize([1,numBins],size(temp0,1));
%     temp2 = discretize(1:numBins,temp1);
%     plot(1:numBlocks,this_data(n,:,idx)*100,'Color',temp0(temp2(n),:));
% end
% temp = cumsum(blocks(1:end-1))+0.5;
% for i=1:length(temp)
%     if i==1 || i==3
%         xline(temp(i),'k-');
%     else
%         xline(temp(i),'k:');
%     end
% end
% xlim([0,numBlocks])
% xlabel('Trial block')
% ytickformat('percentage')
% ylabel('Proportion of sequence cells (rel. to all cells)')
% 
% suptitle(['Sequences across days (',d_info.animals{idx},', green = early sequence cells, blue = late sequence cells)'])

%% Plot sequence cells by block, bin and animal - switch-averaged

default_figure(); 

subplot(1,2,1); hold on;
this_data = seqCellsByTime_switch_num;
for n=1:numBins
    temp0 = flipud(eval('winter'));
    [~,temp1] = discretize([1,numBins],size(temp0,1));
    temp2 = discretize(1:numBins,temp1);
    shadedErrorBar(1:numBlocks_switch,nanmean(this_data(n,:,:),3),nansem(this_data(n,:,:),3),'lineProps',temp0(temp2(n),:));
end
xline(blocks_switch(1)+0.5,'k-');
xlim([0,numBlocks_switch])
xlabel('trial block')
ylabel('Number of sequence cells')

subplot(1,2,2); hold on;
this_data = seqCellsByTime_switch_prop;
for n=1:numBins
    temp0 = flipud(eval('winter'));
    [~,temp1] = discretize([1,numBins],size(temp0,1));
    temp2 = discretize(1:numBins,temp1);
    shadedErrorBar(1:numBlocks_switch,nanmean(this_data(n,:,:),3)*100,nansem(this_data(n,:,:),3)*100,'lineProps',temp0(temp2(n),:));
end
xline(blocks_switch(1)+0.5,'k-');
xlim([0,numBlocks_switch])
xlabel('Trial block')
ytickformat('percentage')
ylabel('Proportion of sequence cells (rel. to all cells)')

suptitle('Sequences across switch (green = early sequence cells, blue = late sequence cells)')


%% Plot sequence cells by block and bin - all 5 days - animal-by-animal

for idx=14%1:d_info.numAnimals
    if ~isempty(d{idx,1})

        default_figure();
        
        subplot(2,3,1); hold on;
        this_data = correct;
        plot(1:numBlocks,this_data(idx,:)*100,'Color','k','LineWidth',3)
        temp = cumsum(blocks(1:end-1))+0.5;
        for i=1:length(temp)
            if i==1 || i==3
                xline(temp(i),'k-');
            else
                xline(temp(i),'k:');
            end
        end
        xlim([0,numBlocks])
        xlabel('Trial block')
        ytickformat('percentage')
        ylabel('Performance (%correct)')

        subplot(2,3,2); hold on;
        this_data = seqCellsByTime_num;
        for n=1:numBins
            temp0 = flipud(eval('winter'));
            [~,temp1] = discretize([1,numBins],size(temp0,1));
            temp2 = discretize(1:numBins,temp1);
            plot(1:numBlocks,this_data(n,:,idx),'Color',temp0(temp2(n),:))
        end
        temp = cumsum(blocks(1:end-1))+0.5;
        for i=1:length(temp)
            if i==1 || i==3
                xline(temp(i),'k-');
            else
                xline(temp(i),'k:');
            end
        end
        xlim([0,numBlocks])
        xlabel('Trial block')
        ylabel('Number of sequence cells')

        subplot(2,3,3); hold on;
        this_data = seqCellsByTime_prop;
        for n=1:numBins
            temp0 = flipud(eval('winter'));
            [~,temp1] = discretize([1,numBins],size(temp0,1));
            temp2 = discretize(1:numBins,temp1);
            plot(1:numBlocks,this_data(n,:,idx)*100,'Color',temp0(temp2(n),:));
        end
        temp = cumsum(blocks(1:end-1))+0.5;
        for i=1:length(temp)
            if i==1 || i==3
                xline(temp(i),'k-');
            else
                xline(temp(i),'k:');
            end
        end
        xlim([0,numBlocks])
        xlabel('Trial block')
        ytickformat('percentage')
        ylabel('Proportion of sequence cells (rel. to all cells)')

        subplot(2,3,5); hold on;
        this_data_x = seqCellsByTime_num;
        this_data_y = correct;
        for n=1:numBins
            temp0 = flipud(eval('winter'));
            [~,temp1] = discretize([1,numBins],size(temp0,1));
            temp2 = discretize(1:numBins,temp1);
            scatter(this_data_x(n,:,idx),this_data_y(idx,:)*100,'MarkerEdgeColor',temp0(temp2(n),:));
            fitLine(this_data_x(n,:,idx)',this_data_y(idx,:)'*100,temp0(temp2(n),:));
        end
        xlabel('Number of sequence cells')
        ylim([0,100])
        ytickformat('percentage')
        ylabel('Performance (%correct)')
        
        subplot(2,3,6); hold on;
        this_data_x = seqCellsByTime_prop;
        this_data_y = correct;
        for n=1:numBins
            temp0 = flipud(eval('winter'));
            [~,temp1] = discretize([1,numBins],size(temp0,1));
            temp2 = discretize(1:numBins,temp1);
            scatter(this_data_x(n,:,idx)*100,this_data_y(idx,:)*100,'MarkerEdgeColor',temp0(temp2(n),:));
            fitLine(this_data_x(n,:,idx)'*100,this_data_y(idx,:)'*100,temp0(temp2(n),:));
        end
        xtickformat('percentage')
        xlabel('Proportion of sequence cells (rel. to all cells)')
        ylim([0,100])
        ytickformat('percentage')
        ylabel('Performance (%correct)')
        
        suptitle(['Sequences across days (',d_info.animals{idx},', green = early sequence cells, blue = late sequence cells)'])
    end
end

