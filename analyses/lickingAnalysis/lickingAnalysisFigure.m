function F = lickingAnalysisFigure(lick,perf,p,info)

F = default_figure([20,0.5,20,9.9]);


%% Response probability

these_labels = categorical({char('A-B'),char('X-Y'),char('A-Y'),char('X-B')});
these_labels = reordercats(these_labels,{char('A-B'),char('X-Y'),char('A-Y'),char('X-B')});

subplot(2,5,1)
hold on
v = bar(these_labels,[perf.type1.correct,perf.type2.correct,perf.type3.incorrect,perf.type4.incorrect]*100);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.AB; v.CData(2,:) = p.col.XY; v.CData(3,:) = p.col.AY; v.CData(4,:) = p.col.XB; set(gca, 'box','off');
yline(lick.engagement*100,'LineStyle',':','LineWidth',2,'Color',p.col.darkGray);
yline(p.lick.crit_engagement*100,'LineStyle',':','LineWidth',2,'Color',p.col.photostim);
ylim([0,100])
ytickformat('percentage');
ylabel('P(licks in response window)')
title('Response probability')


%% Response probability as a function of trial block

subplot(2,5,6)
hold on
plot(smoothdata(perf.blocks_type1.correct,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',2,'Color',p.col.AB)
plot(smoothdata(perf.blocks_type2.correct,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',2,'Color',p.col.XY)
plot(smoothdata(perf.blocks_type3.incorrect,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',2,'Color',p.col.AY)
plot(smoothdata(perf.blocks_type4.incorrect,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',2,'Color',p.col.XB)
ylim([0,100])
ytickformat('percentage');
xlabel('Trial block')
ylabel('P(licks in response window)')
title('Response probability')


%% Response bias

these_labels = categorical({char('general'),char('1st odour'),char('2nd odour')});
these_labels = reordercats(these_labels,{char('general'),char('1st odour'),char('2nd odour')});

subplot(2,5,2)
hold on
v = bar(these_labels,[lick.rbias.general,lick.rbias.odourA,lick.rbias.odourB]);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.darkGray; v.CData(2,:) = p.col.AB; v.CData(3,:) = p.col.XB; set(gca, 'box','off');
yline(p.lick.crit_bias,'LineStyle',':','LineWidth',2,'Color',p.col.photostim);
yline(-p.lick.crit_bias,'LineStyle',':','LineWidth',2,'Color',p.col.photostim);
ylim([-1,1])
ylabel('Response bias index')
title('Response bias')


%% Response bias as a function of trial block

subplot(2,5,7)
hold on
yline(0,'LineStyle',':','LineWidth',2,'Color',p.col.black);
plot(smoothdata(lick.rbias_blocks.general,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',2,'Color',p.col.darkGray)
plot(smoothdata(lick.rbias_blocks.odourA,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',2,'Color',p.col.AB)
plot(smoothdata(lick.rbias_blocks.odourB,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',2,'Color',p.col.XB)
ylim([-1,1])
xlabel('Trial block')
ylabel('Response bias index')
title('Response bias')


%% Lick probability

these_labels = categorical({char('A-B'),char('X-Y'),char('A-Y'),char('X-B')});
these_labels = reordercats(these_labels,{char('A-B'),char('X-Y'),char('A-Y'),char('X-B')});

subplot(2,5,3)
v = bar(these_labels,[lick.lprob.AB,lick.lprob.XY,lick.lprob.AY,lick.lprob.XB]*100);
v.FaceColor = 'flat'; v.CData(1,:) = p.col.AB; v.CData(2,:) = p.col.XY; v.CData(3,:) = p.col.AY; v.CData(4,:) = p.col.XB; set(gca, 'box','off');
ylim([0,100])
ytickformat('percentage');
ylabel('P(licks after 2nd odour onset)')
title('Lick probability')


%% Lick probability as a function of trial block

subplot(2,5,8)
hold on
plot(smoothdata(lick.lprob_blocks.AB,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',2,'Color',p.col.AB)
plot(smoothdata(lick.lprob_blocks.XY,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',2,'Color',p.col.XY)
plot(smoothdata(lick.lprob_blocks.AY,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',2,'Color',p.col.AY)
plot(smoothdata(lick.lprob_blocks.XB,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',2,'Color',p.col.XB)
ylim([0,100])
ytickformat('percentage');
xlabel('Trial block')
ylabel('P(licks after 2nd odour onset)')
title('Lick probability')


%% Lick latency

these_labels = {char('A-B'),char('X-Y'),char('A-Y'),char('X-B')};
this_data = nan(info.task.numTrials/4,4);
this_data(1:length(lick.latency.AB_H),1) = lick.latency.AB_H;
this_data(1:length(lick.latency.XY_H),2) = lick.latency.XY_H;
this_data(1:length(lick.latency.AY_FA),3) = lick.latency.AY_FA;
this_data(1:length(lick.latency.XB_FA),4) = lick.latency.XB_FA;

subplot(2,5,4)
hold on
yline(info.task.trialStructure.tOdour2,'LineStyle',':','LineWidth',2,'Color',p.col.odour);
yline(info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay,'LineStyle',':','LineWidth',2,'Color',p.col.reward);
yline(info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow,'LineStyle',':','LineWidth',2,'Color',p.col.reward);
v = violinplot([this_data(:,1),this_data(:,2),this_data(:,3),this_data(:,4)],these_labels);
v(1).ViolinColor = p.col.AB; v(2).ViolinColor = p.col.XY; v(3).ViolinColor = p.col.AY; v(4).ViolinColor = p.col.XB; 
v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(3).BoxColor = 'k'; v(4).BoxColor = 'k';
ylim([0,info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow])
ylabel('Time since 2nd odour onset (s)')
title('Lick latency')


%% Lick latency as a function of trial block

subplot(2,5,9)
hold on
yline(info.task.trialStructure.tOdour2,'LineStyle',':','LineWidth',2,'Color',p.col.odour);
yline(info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay,'LineStyle',':','LineWidth',2,'Color',p.col.reward);
yline(info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow,'LineStyle',':','LineWidth',2,'Color',p.col.reward);
plot(smoothdata(lick.latency_blocks.AB_H,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',2,'Color',p.col.AB)
plot(smoothdata(lick.latency_blocks.XY_H,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',2,'Color',p.col.XY)
plot(smoothdata(lick.latency_blocks.AY_FA,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',2,'Color',p.col.AY)
plot(smoothdata(lick.latency_blocks.XB_FA,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',2,'Color',p.col.XB)
ylim([0,info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow])
xlabel('Trial block')
ylabel('Time since 2nd odour onset (s)')
title('Lick latency')


%% Lick rigour

these_labels = {char('A-B'),char('X-Y'),char('A-Y'),char('X-B')};
this_data = nan(info.task.numTrials/4,4);
this_data(1:length(lick.rigour.AB_H),1) = lick.rigour.AB_H;
this_data(1:length(lick.rigour.XY_H),2) = lick.rigour.XY_H;
this_data(1:length(lick.rigour.AY_FA),3) = lick.rigour.AY_FA;
this_data(1:length(lick.rigour.XB_FA),4) = lick.rigour.XB_FA;

subplot(2,5,5)
v = violinplot([this_data(:,1),this_data(:,2),this_data(:,3),this_data(:,4)],these_labels);
v(1).ViolinColor = p.col.AB; v(2).ViolinColor = p.col.XY; v(3).ViolinColor = p.col.AY; v(4).ViolinColor = p.col.XB; 
v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(3).BoxColor = 'k'; v(4).BoxColor = 'k';
ylim([0,inf])
ylabel(['Number of licks from 2nd odour onset',newline,'(until reward delivery)'])
%ylabel(['Number of licks since 2nd odour onset'])
title('Lick rigour')


%% Lick rigour as a function of trial block

subplot(2,5,10)
hold on
plot(smoothdata(lick.rigour_blocks.AB_H,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',2,'Color',p.col.AB)
plot(smoothdata(lick.rigour_blocks.XY_H,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',2,'Color',p.col.XY)
plot(smoothdata(lick.rigour_blocks.AY_FA,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',2,'Color',p.col.AY)
plot(smoothdata(lick.rigour_blocks.XB_FA,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',2,'Color',p.col.XB)
ylim([0,inf])
xlabel('Trial block')
ylabel(['Number of licks from 2nd odour onset',newline,'(until reward delivery)'])
%ylabel(['Number of licks since 2nd odour onset'])
title('Lick rigour')


%% Return

if info.stimSession
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType])
else
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-','nostim']);
end

drawnow;
end

