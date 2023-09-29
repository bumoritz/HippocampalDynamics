function F = noStimResponseFigure(s2p_meta,iscell,trg,resp,p,info)

F = default_figure();


%% Example first cluster responses 

% subplot(2,4,1)
% clusterResponseMap(s2p_meta,iscell,trg,resp,resp.noStim.clusterIdx(1,1),1,1); 
% if strcmp(resp.stimType,'seq')
%     title('First cluster (A)')
% elseif strcmp(resp.stimType,'ctrl')
%     title(['Example first cluster (',num2str(resp.noStim.clusterIdx(1,1)),')'])
% end
% 
% subplot(2,4,2)
% clusterResponseMap(s2p_meta,iscell,trg,resp,resp.noStim.clusterIdx(1,2),1,1); 
% if strcmp(resp.stimType,'seq')
%     title('First cluster (X)')
% elseif strcmp(resp.stimType,'ctrl')
%     title(['Example first cluster (',num2str(resp.noStim.clusterIdx(1,2)),')'])
% end


%% Photostimulation responders

these_labels = categorical({char('pFDR'),char('pFDR, 0.5 z'),char('pFDR, 1 z')});
these_labels = reordercats(these_labels,{char('pFDR'),char('pFDR, 0.5 z'),char('pFDR, 1 z')});

subplot(2,4,3)
hold on
bar(these_labels,[resp.noStim.responders_ext.numRespTargeted,resp.noStim.responders_0d5z_ext.numRespTargeted,resp.noStim.responders_1z_ext.numRespTargeted],'FaceColor',p.col.photostim)
yline(resp.noStim.responders_ext.numIdentifiedTargeted,'LineStyle',':','LineWidth',2);
ylim([0,resp.noStim.responders_ext.numAllTargeted])
yticks([0:40:240])
ylabel('Number of responsive targeted cells')
title('Photostimulation target responses')

subplot(2,4,4)
hold on
bar(these_labels,[resp.noStim.responders_ext.proportion_RespNonTargeted_RespTargeted,resp.noStim.responders_0d5z_ext.proportion_RespNonTargeted_RespTargeted,resp.noStim.responders_1z_ext.proportion_RespNonTargeted_RespTargeted],'FaceColor',p.col.darkGray)
yline(0,'LineStyle','-');
ylabel('Responsive non-targeted cells per targeted cell')
title('Photostimulation non-target responses')


%% Photostimulation amplitude

subplot(2,4,5)

temp = resp.noStim.avgAct_net(:);
hold on
histogram(temp(find(resp.noStim.stats_sig_pos(:)==1 & resp.noStim.responders_main.targeted(:)==1)),'Normalization','probability','BinWidth',0.25,'FaceColor',p.col.photostim)
histogram(temp(find(resp.noStim.stats_sig_pos(:)==1 & resp.noStim.responders_main.targeted(:)==0)),'Normalization','probability','BinWidth',0.25,'FaceColor',p.col.darkGray)
xline(p.resp.mainAmplitudeCritierion,'LineStyle',':','LineWidth',2);
hold off

xlabel('Response amplitude\newline(\Delta_{resp-base} z-scored \DeltaF/F)')
ylabel('Proportion')
legend('Responsive targeted cells','Responsive non-targeted cells')
title('Photostimulation amplitude')


%% Photostimulation stability

% subplot(2,4,6)
% 
% shadedErrorBar(1:size(resp.trgRespAmps_cells_bw,2),nanmean(resp.trgRespAmps_cells_bw,1),nansem(resp.trgRespAmps_cells_bw,1),'lineProps',p.col.black);
% 
% xlim([0,info.task.numBlocks])
% ylim([0,inf])
% xlabel('Trial block')
% ylabel('Mean response amplitude of targeted cells\newline(\Delta_{resp-base} z-scored \DeltaF/F)')
% title('Photostimulation stability')


%% Photostimulation resolution

subplot(2,4,7)

binsize = 2.5;
hold on
for i=1:length(p.resp.amplitudeCritiera)
    if isnan(p.resp.amplitudeCritiera(i))
        fieldname = char('responders');
    else
        fieldname = char(['responders_',num2str(p.resp.amplitudeCritiera(i)),'z']);
        fieldname(find(fieldname=='.'))='d';
    end
    
    temp = resp.noStim.(fieldname).respAll(:);
    temp2 = discretize(trg.dist_closestLaser(:),[0:binsize:50]);
    responseProbability = nan(nanmax(temp2),1);
    for j=1:nanmax(temp2)
        responseProbability(j) = nanmean(temp(find(temp2==j)));
    end
    plot(responseProbability*100,'LineWidth',2)
end

xticks([0:4:length([0:binsize:50])]+1)
xticklabels({'0','10','20','30','40','50'})
xlim([1,length([0:binsize:50])])
ylim([0,100])
ytickformat('percentage');
xlabel('Distance from closest laser beamlet (\mum)\newline(2.5\mum bin size)')
ylabel('Response probability')
legend('pFDR','pFDR, 0.5 z','pFDR, 1 z')
title('Photostimulation resolution')


%% Photostimulation followers

subplot(2,4,8)

binsize = 10;
hold on

temp = resp.noStim.stats_sig_pos(:);
temp2 = discretize(trg.dist_closestLaser(:),[0:binsize:400]);
responseProbability = nan(nanmax(temp2),1);
for j=1:nanmax(temp2)
    responseProbability(j) = nanmean(temp(find(temp2==j)));
end
plot(responseProbability*100,'LineWidth',2,'Color','r')
this_upper = responseProbability(3)*100;

temp = resp.noStim.stats_sig_neg(:);
temp2 = discretize(trg.dist_closestLaser(:),[0:binsize:400]);
responseProbability = nan(nanmax(temp2),1);
for j=1:nanmax(temp2)
    responseProbability(j) = nanmean(temp(find(temp2==j)));
end
plot(responseProbability*100,'LineWidth',2,'Color','b')

xticks([0:(length([0:binsize:400]))/8:length([0:binsize:400])]+1)
xticklabels({'0','50','100','150','200','250','300','350','400'})
xlim([1,length([0:binsize:400])+1])
ylim([0,this_upper])
ytickformat('percentage');
xlabel('Distance from closest laser beamlet (\mum)\newline(25\mum bin size)')
ylabel('Response probability')
legend('positive','negative')
title('Photostimulation followers')


%% Return

suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType])
drawnow;
end
