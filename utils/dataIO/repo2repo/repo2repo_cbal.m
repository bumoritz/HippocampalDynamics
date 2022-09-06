function path = repo2repo_cbal(info,ops,p,path)

%% Load trg files

disp('--- Loading trg and meta files...')

trg1 = load([path.root_repoX,info.animal,'\',info.animal,'_',info.date(1,:),'\',info.animal,'_',info.date(1,:),'_','trg_online.mat']);
trg1 = trg1.trg_online;
trg2 = load([path.root_repoX,info.animal,'\',info.animal,'_',info.date(2,:),'\',info.animal,'_',info.date(2,:),'_','trg_online.mat']);
trg2 = trg2.trg_online;
if exist([path.filepart_out,'meta.mat'])==2
    load([path.filepart_out,'meta.mat']);
end


%% Kruskal-Wallis tests

% calculate p values - group comparisons and cluster comparisons
pval=NaN(18,2);
for j=1:18
    if j<18
        pval(j,1) = kruskalwallis([trg1.features.cellw1(:,j),trg1.features.cellw2(:,j),trg2.features.cellw1(:,j),trg2.features.cellw2(:,j)],[],'off');
    end
    pval(j,2) = kruskalwallis([trg1.features.clusterw1(:,j),trg1.features.clusterw2(:,j),trg2.features.clusterw1(:,j),trg2.features.clusterw2(:,j)],[],'off');
end


%% Plot - Counterbalancing of target groups

nrows = 4;
ncols = 5;
these_labels = {trg1.p.groupNames(1),trg1.p.groupNames(2),trg2.p.groupNames(1),trg2.p.groupNames(2)};
plot_positions = [1,6,2,7,12,17,3,8,13,18,11,4,9,14,5,10,15];
exc = [1,1,0,0,0,0,0,2,2,0,0,0,0,0,0,0,0];

F1 = default_figure([20,0.5,20,9.9]);
for n=1:17
    subplot(nrows,ncols,plot_positions(n));
    if exc(n)==2
        v = violinplot([trg1.features.cellw1(:,n),trg1.features.cellw2(:,n),trg2.features.cellw1(:,n),trg2.features.cellw2(:,n)]*100,these_labels);
        ytickformat('percentage');
    else
        v = violinplot([trg1.features.cellw1(:,n),trg1.features.cellw2(:,n),trg2.features.cellw1(:,n),trg2.features.cellw2(:,n)],these_labels);
    end
    if strcmp(trg1.p.stimType,'seq')
        v(1).ViolinColor = p.col.seq; v(2).ViolinColor = p.col.seq; v(3).ViolinColor = p.col.ctrl; v(4).ViolinColor = p.col.ctrl; 
    elseif strcmp(trg2.p.stimType,'seq')
        v(1).ViolinColor = p.col.ctrl; v(2).ViolinColor = p.col.ctrl; v(3).ViolinColor = p.col.seq; v(4).ViolinColor = p.col.seq; 
    end
    v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(3).BoxColor = 'k'; v(4).BoxColor = 'k';
    ylabel(trg1.features.names(n));
    title(['kwtest: p=',num2str(pval(n,1),2)],'FontSize',10,'FontWeight','normal');
    
    if exc(n)==1
        ylim([0,floor(trg1.info.FOVsize)*trg1.info.numPixels/trg1.info.FOVsize]);
        yticks([0:(floor(trg1.info.FOVsize)*trg1.info.numPixels/trg1.info.FOVsize)/4:floor(trg1.info.FOVsize)*trg1.info.numPixels/trg1.info.FOVsize])
        yticklabels(strsplit(num2str([0:floor(trg1.info.FOVsize)/4:floor(trg1.info.FOVsize)])))
    end 
end
suptitle('Counterbalancing of target groups');


%% Plot - Counterbalancing of target clusters

nrows = 4;
ncols = 5;
these_labels = categorical({char(trg1.p.groupNames(1)),char(trg1.p.groupNames(2)),char(trg2.p.groupNames(1)),char(trg2.p.groupNames(2))});
these_labels = reordercats(these_labels,{char(trg1.p.groupNames(1)),char(trg1.p.groupNames(2)),char(trg2.p.groupNames(1)),char(trg2.p.groupNames(2))});
plot_positions = [1,6,2,7,12,17,3,8,13,18,11,4,9,14,5,10,15,16];
exc = [1,1,0,3,3,3,0,2,2,4,0,0,0,0,3,3,3,0];

F2 = default_figure([20,0.5,20,9.9]);
for n=1:18
    subplot(nrows,ncols,plot_positions(n));
    hold on
    if exc(n)==2
        v = bar(these_labels,[mean(trg1.features.clusterw1(:,n)),mean(trg1.features.clusterw2(:,n)),mean(trg2.features.clusterw1(:,n)),mean(trg2.features.clusterw2(:,n))]*100);
        for i=1:trg1.p.numClustersPerGroup
            plot(these_labels(1:2),[trg1.features.clusterw1(i,n),trg1.features.clusterw2(i,n)]*100,'-k')
        end
        for i=1:trg2.p.numClustersPerGroup
            plot(these_labels(3:4),[trg2.features.clusterw1(i,n),trg2.features.clusterw2(i,n)]*100,'-k')
        end
        ytickformat('percentage');
    else
        v = bar(these_labels,[mean(trg1.features.clusterw1(:,n)),mean(trg1.features.clusterw2(:,n)),mean(trg2.features.clusterw1(:,n)),mean(trg2.features.clusterw2(:,n))]);
        for i=1:trg1.p.numClustersPerGroup
            plot(these_labels(1:2),[trg1.features.clusterw1(i,n),trg1.features.clusterw2(i,n)],'-k')
        end
        for i=1:trg2.p.numClustersPerGroup
            plot(these_labels(3:4),[trg2.features.clusterw1(i,n),trg2.features.clusterw2(i,n)],'-k')
        end
    end
    hold off
    v.FaceColor = 'flat'; 
    if strcmp(trg1.p.stimType,'seq')
        v.CData(1,:) = p.col.seq; v.CData(2,:) = p.col.seq; v.CData(3,:) = p.col.ctrl; v.CData(4,:) = p.col.ctrl;
    elseif strcmp(trg2.p.stimType,'seq')
        v.CData(1,:) = p.col.ctrl; v.CData(2,:) = p.col.ctrl; v.CData(3,:) = p.col.seq; v.CData(4,:) = p.col.seq;
    end
    ylabel(trg1.features.names(n));
    title(['kwtest: p=',num2str(pval(n,2),2)],'FontSize',10,'FontWeight','normal');
    
    if exc(n)==1
        ylim([0,floor(trg1.info.FOVsize)*trg1.info.numPixels/trg1.info.FOVsize]);
        yticks([0:(floor(trg1.info.FOVsize)*trg1.info.numPixels/trg1.info.FOVsize)/4:floor(trg1.info.FOVsize)*trg1.info.numPixels/trg1.info.FOVsize])
        yticklabels(strsplit(num2str([0:floor(trg1.info.FOVsize)/4:floor(trg1.info.FOVsize)])))
    elseif exc(n)==4
        ylim([floor(min(min([trg1.features.clusterw1(:,n),trg1.features.clusterw2(:,n),trg2.features.clusterw1(:,n),trg2.features.clusterw2(:,n)]))*10)/10,ceil(max(max([trg1.features.clusterw1(:,n),trg1.features.clusterw2(:,n),trg2.features.clusterw1(:,n),trg2.features.clusterw2(:,n)]))*10)/10]);
    elseif exc(n)==3
        ylim([floor(min(min([trg1.features.clusterw1(:,n),trg1.features.clusterw2(:,n),trg2.features.clusterw1(:,n),trg2.features.clusterw2(:,n)]))*100)/100,ceil(max(max([trg1.features.clusterw1(:,n),trg1.features.clusterw2(:,n),trg2.features.clusterw1(:,n),trg2.features.clusterw2(:,n)]))*100)/100]);
    elseif exc(n)==2
        ylim([floor(min(min([trg1.features.clusterw1(:,n),trg1.features.clusterw2(:,n),trg2.features.clusterw1(:,n),trg2.features.clusterw2(:,n)]*100))),ceil(max(max([trg1.features.clusterw1(:,n),trg1.features.clusterw2(:,n),trg2.features.clusterw1(:,n),trg2.features.clusterw2(:,n)]*100)))]);
    else
        ylim([floor(min(min([trg1.features.clusterw1(:,n),trg1.features.clusterw2(:,n),trg2.features.clusterw1(:,n),trg2.features.clusterw2(:,n)]))),ceil(max(max([trg1.features.clusterw1(:,n),trg1.features.clusterw2(:,n),trg2.features.clusterw1(:,n),trg2.features.clusterw2(:,n)])))]);
    end    
end
suptitle('Counterbalancing of target clusters');


%% Save

savefig(F1,[path.filepart_out,'plots\',info.animal,'_cbal_cellw.fig']);
saveas(F1,[path.filepart_out,'plots\',info.animal,'_cbal_cellw.png']);
savefig(F2,[path.filepart_out,'plots\',info.animal,'_cbal_clusterw.fig']);
saveas(F2,[path.filepart_out,'plots\',info.animal,'_cbal_clusterw.png']);
disp(['--- Saved cbal figures to ',path.filepart_out,'plots.'])
if ~ops.cbal.showFigures
    close all;
end

meta.cbal.pval = pval;
meta = orderfields(meta);
save([path.filepart_out,'meta.mat'],'meta','-v7.3');
disp(['--- Overwrote meta file in repo at ',[path.filepart_out,'meta.mat'],'.'])
    
end