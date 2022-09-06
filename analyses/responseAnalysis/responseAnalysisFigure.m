function F = responseAnalysisFigure(s2p_meta,iscell,trg,resp,p,info)

F = default_figure([20,0.5,20,9.9]);


%% Example cluster responses

subplot(2,4,1)
clusterResponseMap(s2p_meta,iscell,trg,resp,randi(size(trg.idcs_targetedCells,2)),1); 
title('Example photostimulation cluster')


%% Photostimulation responders

these_labels = categorical({char('pFDR'),char('pFDR, 0.5 z'),char('pFDR, 1 z')});
these_labels = reordercats(these_labels,{char('pFDR'),char('pFDR, 0.5 z'),char('pFDR, 1 z')});

subplot(2,4,2)
hold on
bar(these_labels,[resp.responders_ext.numRespTargeted,resp.responders_0d5z_ext.numRespTargeted,resp.responders_1z_ext.numRespTargeted],'FaceColor',p.col.photostim)
yline(resp.responders_ext.numIdentifiedTargeted,'LineStyle',':','LineWidth',2);
ylim([0,resp.responders_ext.numAllTargeted])
yticks([0:40:240])
ylabel('Number of responsive targeted cells')
title('Photostimulation target responses')

subplot(2,4,3)
hold on
bar(these_labels,[resp.responders_ext.proportion_RespNonTargeted_RespTargeted,resp.responders_0d5z_ext.proportion_RespNonTargeted_RespTargeted,resp.responders_1z_ext.proportion_RespNonTargeted_RespTargeted],'FaceColor',p.col.darkGray)
yline(0,'LineStyle','-');
ylabel('Responsive non-targeted cells per targeted cell')
title('Photostimulation non-target responses')

subplot(2,4,4)
hold on
v=bar(these_labels,[resp.responders_ext.specificity_respTargeted,resp.responders_ext.specificity_respIpsiGroup-resp.responders_ext.specificity_respTargeted,resp.responders_ext.specificity_respContraGroup,resp.responders_ext.specificity_respNoGroupCells;...
    resp.responders_0d5z_ext.specificity_respTargeted,resp.responders_0d5z_ext.specificity_respIpsiGroup-resp.responders_0d5z_ext.specificity_respTargeted,resp.responders_0d5z_ext.specificity_respContraGroup,resp.responders_0d5z_ext.specificity_respNoGroupCells;...
    resp.responders_1z_ext.specificity_respTargeted,resp.responders_1z_ext.specificity_respIpsiGroup-resp.responders_1z_ext.specificity_respTargeted,resp.responders_1z_ext.specificity_respContraGroup,resp.responders_1z_ext.specificity_respNoGroupCells]);
v(1).FaceColor=p.col.photostim; v(2).FaceColor=mean([p.col.photostim;p.col.gray],1); v(3).FaceColor=mean([p.col.photostim;p.col.darkGray],1); v(4).FaceColor=p.col.darkGray;

ylim([0,1])
ylabel('Proportion of all responsive cells')
legend('targeted cells','same stim group (non-targeted)','other stim group','no stim group')
title('Photostimulation specificity')


%% Photostimulation amplitude

subplot(2,4,5)

temp = resp.avgAct_net(:);
hold on
histogram(temp(find(resp.stats_sig_pos(:)==1 & resp.responders_main.targeted(:)==1)),'Normalization','probability','BinWidth',0.25,'FaceColor',p.col.photostim)
histogram(temp(find(resp.stats_sig_pos(:)==1 & resp.responders_main.targeted(:)==0)),'Normalization','probability','BinWidth',0.25,'FaceColor',p.col.darkGray)
xline(p.resp.mainAmplitudeCritierion,'LineStyle',':','LineWidth',2);
hold off

xlabel('Response amplitude\newline(\Delta_{resp-base} z-scored \DeltaF/F)')
ylabel('Proportion')
legend('Responsive targeted cells','Responsive non-targeted cells')
title('Photostimulation amplitude')


%% Photostimulation stability

subplot(2,4,6)

shadedErrorBar(1:size(resp.trgRespAmps_cells_bw,2),nanmean(resp.trgRespAmps_cells_bw,1),nansem(resp.trgRespAmps_cells_bw,1),'lineProps',p.col.black);

xlim([0,info.task.numBlocks])
ylim([0,inf])
xlabel('Trial block')
ylabel('Mean response amplitude of targeted cells\newline(\Delta_{resp-base} z-scored \DeltaF/F)')
title('Photostimulation stability')


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
    
    temp = resp.(fieldname).respAll(:);
    temp2 = discretize(trg.dist_closestLaser(:),[0:binsize:50]);
    responseProbability = nan(nanmax(temp2),1);
    for j=1:nanmax(temp2)
        responseProbability(j) = nanmean(temp(find(temp2==j)));
    end
    plot(responseProbability*100,'LineWidth',2)
end

xticks([0:4:length([0:binsize:50])-1])
xticklabels({'0','10','20','30','40','50'})
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

temp = resp.stats_sig_pos(:);
temp2 = discretize(trg.dist_closestLaser(:),[0:binsize:400]);
responseProbability = nan(nanmax(temp2),1);
for j=1:nanmax(temp2)
    responseProbability(j) = nanmean(temp(find(temp2==j)));
end
plot(responseProbability*100,'LineWidth',2,'Color','r')
this_upper = responseProbability(3)*100;

temp = resp.stats_sig_neg(:);
temp2 = discretize(trg.dist_closestLaser(:),[0:binsize:400]);
responseProbability = nan(nanmax(temp2),1);
for j=1:nanmax(temp2)
    responseProbability(j) = nanmean(temp(find(temp2==j)));
end
plot(responseProbability*100,'LineWidth',2,'Color','b')

xticks([0:(length([0:binsize:400])-1)/8:length([0:binsize:400])-1])
xticklabels({'0','50','100','150','200','250','300','350','400'})
xlim([0,length([0:binsize:400])-1])
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
