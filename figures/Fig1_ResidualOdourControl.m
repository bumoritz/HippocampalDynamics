%% Select data

clc; clear; 

p = get_p;

info.animal     = ["Guzman";"Guzman";"AlCapone";"AlCapone";"Scarface";"Scarface";"Gambino";"Gambino";"Escobar";"Escobar";"Baskin";"Baskin";...
    "Costello";"Gotti";"Gotti";"Guzman";"Guzman";"AlCapone";"AlCapone";"Baskin";"Baskin";"DocAntle";"DocAntle";...
    "Guzman";"Guzman";"AlCapone";"AlCapone";"Gotti";"Gotti";"Escobar";"Escobar";"Gambino";"Gambino";"Glover";"Glover"];
info.date       = ['20200805';'20200807';'20200805';'20200807';'20200815';'20200817';'20200818';'20200821';'20200825';'20200831';'20200921';'20201002';...
    '20200807';'20200814';'20200815';'20200815';'20200816';'20200813';'20200814';'20200910';'20200911';'20200917';'20200920';...
    '20200820';'20200824';'20200820';'20200821';'20200825';'20200901';'20200909';'20200910';'20200911';'20200912';'20200925';'20200930'];
info.lowConc    = [0.001,0,0.001,0,0.001,0,0.001,0,0.001,0,0.001,0,...
    0.001,0.001,0,0.001,0,0.001,0,0.001,0,0.001,0,...
    0.001,0,0.001,0,0.001,0,0.001,0,0.001,0,0.001,0];
info.set        = [1,1,1,1,1,1,1,1,1,1,1,1,...
    2,2,2,2,2,2,2,2,2,2,2,...
    3,3,3,3,3,3,3,3,3,3,3,3];

info.dataRoot       = 'E:\SniffinHippo\RepoO\';
info.analysisRoot   = 'E:\SniffinHippo\RepoO\';

path.root_summary = 'C:\SniffinHippo\Summary\';
save_root_fig = [path.root_summary,'figures\Fig1_fig\'];
save_root_png = [path.root_summary,'figures\Fig1_png\'];
save_root_pdf = [path.root_summary,'figures\Fig1_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig1_txt\'];


%% Get data from every session

conc = [];
perf = [];
for sess=1:length(info.animal)

    %% Load data

    trl = load([info.dataRoot,'ResidualOdorControl\',char(info.animal(sess)),'_',char(info.date(sess,:)),'.mat']);
    trl = trl.trl;

    %% Get session information
    numBlocks = 5;
    type1(sess,1) = mean(trl.performance.blocks_type1.correct(1:numBlocks));
    type2(sess,1) = mean(trl.performance.blocks_type2.correct(1:numBlocks));
    type3(sess,1) = mean(trl.performance.blocks_type3.correct(1:numBlocks));
    type4(sess,1) = mean(trl.performance.blocks_type4.correct(1:numBlocks));
    type5(sess,1) = mean(trl.performance.blocks_type5.correct(1:numBlocks));
    type6(sess,1) = mean(trl.performance.blocks_type6.correct(1:numBlocks));
    type7(sess,1) = mean(trl.performance.blocks_type7.correct(1:numBlocks));
    type8(sess,1) = mean(trl.performance.blocks_type8.correct(1:numBlocks));
%     type1(sess,1) = nanmean(trl.performance.blocks_type1.correct(1:numBlocks));
%     type2(sess,1) = nanmean(trl.performance.blocks_type2.correct(1:numBlocks));
%     type3(sess,1) = nanmean(trl.performance.blocks_type3.correct(1:numBlocks));
%     type4(sess,1) = nanmean(trl.performance.blocks_type4.correct(1:numBlocks));
%     type5(sess,1) = nanmean(trl.performance.blocks_type5.correct(1:numBlocks));
%     type6(sess,1) = nanmean(trl.performance.blocks_type6.correct(1:numBlocks));
%     type7(sess,1) = nanmean(trl.performance.blocks_type7.correct(1:numBlocks));
%     type8(sess,1) = nanmean(trl.performance.blocks_type8.correct(1:numBlocks));
end

high_A = (type1+type3)/2;
high_X = (type2+type4)/2;
high = (high_A+high_X)/2;
low_A = (type5+type7)/2;
low_X = (type6+type8)/2;
low = (low_A+low_X)/2;
all_A = [high_A,low_A];
all_X = [high_X,low_X];
all = [high,low];


%% Final-type plot

figure;

subplot(2,3,1)
hold on
data = all(info.set==1&info.lowConc==0.001,:)*100;
bar([1,2],mean(data,1),'FaceColor',p.col.odour)
for i=1:size(data,1)
    plot([1,2],data(i,:),'o-k')
end
yline(50,'Color',p.col.darkGray,'LineStyle',':');
hold off
xlim([0,3])
xticks([1,2])
xticklabels({'100%','0.1%'})
ylim([40,100])
yticks([50,75,100])
ytickformat('percentage')
ylabel('Performance')
title('Set i - residual odour control')

subplot(2,3,2)
hold on
data = all(info.set==2&info.lowConc==0.001,:)*100;
bar([1,2],mean(data,1),'FaceColor',p.col.odour)
for i=1:size(data,1)
    plot([1,2],data(i,:),'o-k')
end
yline(50,'Color',p.col.darkGray,'LineStyle',':');
hold off
xlim([0,3])
xticks([1,2])
xticklabels({'100%','0.1%'})
ylim([40,100])
yticks([50,75,100])
ytickformat('percentage')
ylabel('Performance')
title('Set ii - residual odour control')

subplot(2,3,3)
hold on
data = all(info.set==3&info.lowConc==0.001,:)*100;
bar([1,2],mean(data,1),'FaceColor',p.col.odour)
for i=1:size(data,1)
    plot([1,2],data(i,:),'o-k')
end
yline(50,'Color',p.col.darkGray,'LineStyle',':');
hold off
xlim([0,3])
xticks([1,2])
xticklabels({'100%','0.1%'})
ylim([40,100])
yticks([50,75,100])
ytickformat('percentage')
ylabel('Performance')
title('Set iii - residual odour control')

subplot(2,3,4)
hold on
data = all(info.set==1&info.lowConc==0,:)*100;
bar([1,2],mean(data,1),'FaceColor',p.col.odour)
for i=1:size(data,1)
    plot([1,2],data(i,:),'o-k')
end
yline(50,'Color',p.col.darkGray,'LineStyle',':');
hold off
xlim([0,3])
xticks([1,2])
xticklabels({'100%','0%'})
ylim([40,100])
yticks([50,75,100])
ytickformat('percentage')
ylabel('Performance')
title('Set i - no odour control')

subplot(2,3,5)
hold on
data = all(info.set==2&info.lowConc==0,:)*100;
bar([1,2],mean(data,1),'FaceColor',p.col.odour)
for i=1:size(data,1)
    plot([1,2],data(i,:),'o-k')
end
yline(50,'Color',p.col.darkGray,'LineStyle',':');
hold off
xlim([0,3])
xticks([1,2])
xticklabels({'100%','0%'})
ylim([40,100])
yticks([50,75,100])
ytickformat('percentage')
ylabel('Performance')
title('Set ii - no odour control')

subplot(2,3,6)
hold on
data = all(info.set==3&info.lowConc==0,:)*100;
bar([1,2],mean(data,1),'FaceColor',p.col.odour)
for i=1:size(data,1)
    plot([1,2],data(i,:),'o-k')
end
yline(50,'Color',p.col.darkGray,'LineStyle',':');
hold off
xlim([0,3])
xticks([1,2])
xticklabels({'100%','0%'})
ylim([40,100])
yticks([50,75,100])
ytickformat('percentage')
ylabel('Performance')
title('Set iii - no odour control')


%% Fig1_ResidualOdourControl_i

these_labels = categorical({'Normal','Low','normal','None'});
these_labels = reordercats(these_labels,{'Normal','Low','normal','None'});
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;

this_data_1 = all(info.set==1&info.lowConc==0.001,:)*100;
this_data_2 = all(info.set==1&info.lowConc==0,:)*100;
v = bar(these_labels,[nanmean(this_data_1,1),nanmean(this_data_2,1)]);
%v.FaceColor = 'flat'; v.CData(1,:) = p.col.odour; v.CData(2,:) = mean([p.col.odour;p.col.white]); v.CData(3,:) = p.col.odour; v.CData(4,:) = p.col.gray; v.EdgeColor = 'none'; 
v.FaceColor = 'flat'; v.CData(1,:) = p.col.odour; v.CData(2,:) = p.col.gray; v.CData(3,:) = p.col.odour; v.CData(4,:) = p.col.darkGray; v.EdgeColor = 'none'; 
yline(50,':');
for i=1:size(this_data_1,1)
    plot(these_labels(1:2),this_data_1(i,:),'-k','LineWidth',1)
end
for i=1:size(this_data_2,1)
    plot(these_labels(3:4),this_data_2(i,:),'-k','LineWidth',1)
end
xlabel({'Concentration of 1st odor '})
ylim([35,100])
yticks([50,100])
ylabel({'Performance'})

savefig(F,[save_root_fig,'\Fig1_ResidualOdourControl_i.fig']);
saveas(F,[save_root_png,'\Fig1_ResidualOdourControl_i.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_ResidualOdourControl_i.pdf']); set(gcf,'Color',[1,1,1])

[temp1,~,temp2] = signrank(this_data_1(:,2)-50)
[temp1,~,temp2] = signrank(this_data_2(:,2)-50)
[temp1,~,temp2] = signrank(this_data_1(:,2)-50,this_data_2(:,2)-50)


%% Fig1_ResidualOdourControl_ii

these_labels = categorical({'Normal','Low','normal','None'});
these_labels = reordercats(these_labels,{'Normal','Low','normal','None'});
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;

this_data_1 = all(info.set==2&info.lowConc==0.001,:)*100;
this_data_2 = all(info.set==2&info.lowConc==0,:)*100;
v = bar(these_labels,[nanmean(this_data_1,1),nanmean(this_data_2,1)]);
%v.FaceColor = 'flat'; v.CData(1,:) = p.col.odour; v.CData(2,:) = mean([p.col.odour;p.col.white]); v.CData(3,:) = p.col.odour; v.CData(4,:) = p.col.gray; v.EdgeColor = 'none'; 
v.FaceColor = 'flat'; v.CData(1,:) = p.col.odour; v.CData(2,:) = p.col.gray; v.CData(3,:) = p.col.odour; v.CData(4,:) = p.col.darkGray; v.EdgeColor = 'none'; 
yline(50,':');
for i=1:size(this_data_1,1)
    plot(these_labels(1:2),this_data_1(i,:),'-k','LineWidth',1)
end
for i=1:size(this_data_2,1)
    plot(these_labels(3:4),this_data_2(i,:),'-k','LineWidth',1)
end
xlabel({'Concentration of 1st odor '})
ylim([35,100])
yticks([50,100])
ylabel({'Performance'})

savefig(F,[save_root_fig,'\Fig1_ResidualOdourControl_ii.fig']);
saveas(F,[save_root_png,'\Fig1_ResidualOdourControl_ii.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_ResidualOdourControl_ii.pdf']); set(gcf,'Color',[1,1,1])

[temp1,~,temp2] = signrank(this_data_1(:,2)-50)
[temp1,~,temp2] = signrank(this_data_2(:,2)-50)
[temp1,~,temp2] = signrank(this_data_1(2:end,2)-50,this_data_2(:,2)-50)


%% Fig1_ResidualOdourControl_iii

these_labels = categorical({'Normal','Low','normal','None'});
these_labels = reordercats(these_labels,{'Normal','Low','normal','None'});
F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(34)]); hold on;

this_data_1 = all(info.set==3&info.lowConc==0.001,:)*100;
this_data_2 = all(info.set==3&info.lowConc==0,:)*100;
v = bar(these_labels,[nanmean(this_data_1,1),nanmean(this_data_2,1)]);
%v.FaceColor = 'flat'; v.CData(1,:) = p.col.odour; v.CData(2,:) = mean([p.col.odour;p.col.white]); v.CData(3,:) = p.col.odour; v.CData(4,:) = p.col.gray; v.EdgeColor = 'none'; 
v.FaceColor = 'flat'; v.CData(1,:) = p.col.odour; v.CData(2,:) = p.col.gray; v.CData(3,:) = p.col.odour; v.CData(4,:) = p.col.darkGray; v.EdgeColor = 'none'; 
yline(50,':');
for i=1:size(this_data_1,1)
    plot(these_labels(1:2),this_data_1(i,:),'-k','LineWidth',1)
end
for i=1:size(this_data_2,1)
    plot(these_labels(3:4),this_data_2(i,:),'-k','LineWidth',1)
end
xlabel({'Concentration of 1st odor '})
ylim([35,100])
yticks([50,100])
ylabel({'Performance'})

savefig(F,[save_root_fig,'\Fig1_ResidualOdourControl_iii.fig']);
saveas(F,[save_root_png,'\Fig1_ResidualOdourControl_iii.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig1_ResidualOdourControl_iii.pdf']); set(gcf,'Color',[1,1,1])

[temp1,~,temp2] = signrank(this_data_1(:,2)-50)
[temp1,~,temp2] = signrank(this_data_2(:,2)-50)
[temp1,~,temp2] = signrank(this_data_1(:,2)-50,this_data_2(:,2)-50)










