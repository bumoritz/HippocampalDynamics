function F = lickingAnalysisFigure_stim(lick,perf,p,info)

F = default_figure([20,0.5,20,9.9]);


%% Response probability

these_labels = categorical({char('A-B'),char('X-Y'),char('A-Y'),char('X-B')});
these_labels = reordercats(these_labels,{char('A-B'),char('X-Y'),char('A-Y'),char('X-B')});

subplot(2,5,1)
hold on
v = bar(these_labels,[perf.type1.correct_stim,perf.type1.correct_catch;perf.type2.correct_stim,perf.type2.correct_catch;perf.type3.incorrect_stim,perf.type3.incorrect_catch;perf.type4.incorrect_stim,perf.type4.incorrect_catch]*100);
v(1).FaceColor = 'flat'; v(2).FaceColor = 'flat'; set(gca,'box','off');
v(1).CData(1,:) = p.col.AB; v(2).CData(1,:) = mean([p.col.AB;p.col.gray],1); 
v(1).CData(2,:) = p.col.XY; v(2).CData(2,:) = mean([p.col.XY;p.col.gray],1); 
v(1).CData(3,:) = p.col.AY; v(2).CData(3,:) = mean([p.col.AY;p.col.gray],1); 
v(1).CData(4,:) = p.col.XB; v(2).CData(4,:) = mean([p.col.XB;p.col.gray],1); 
yline(lick.engagement*100,'LineStyle',':','LineWidth',2,'Color',p.col.darkGray);
yline(p.lick.crit_engagement*100,'LineStyle',':','LineWidth',2,'Color',p.col.photostim);
ylim([0,100])
ytickformat('percentage');
ylabel('P(licks in response window)')
title('Response probability')


%% Response probability as a function of trial block

subplot(2,5,6)
hold on
plot(smoothdata(perf.blocks_type1.correct_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineStyle','--','LineWidth',1,'Color',p.col.AB)
plot(smoothdata(perf.blocks_type2.correct_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineStyle','--','LineWidth',1,'Color',p.col.XY)
plot(smoothdata(perf.blocks_type3.incorrect_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineStyle','--','LineWidth',1,'Color',p.col.AY)
plot(smoothdata(perf.blocks_type4.incorrect_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineStyle','--','LineWidth',1,'Color',p.col.XB)
plot(smoothdata(perf.blocks_type1.correct_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',1,'Color',p.col.AB)
plot(smoothdata(perf.blocks_type2.correct_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',1,'Color',p.col.XY)
plot(smoothdata(perf.blocks_type3.incorrect_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',1,'Color',p.col.AY)
plot(smoothdata(perf.blocks_type4.incorrect_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',1,'Color',p.col.XB)
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
v = bar(these_labels,[lick.rbias.general_stim,lick.rbias.general_catch;lick.rbias.odourA_stim,lick.rbias.odourA_catch;lick.rbias.odourB_stim,lick.rbias.odourB_catch]);
v(1).FaceColor = 'flat'; v(2).FaceColor = 'flat'; set(gca,'box','off');
v(1).CData(1,:) = p.col.darkGray; v(2).CData(1,:) = mean([p.col.darkGray;p.col.gray],1); 
v(1).CData(2,:) = p.col.AB; v(2).CData(2,:) = mean([p.col.AB;p.col.gray],1); 
v(1).CData(3,:) = p.col.XB; v(2).CData(3,:) = mean([p.col.XB;p.col.gray],1); 
yline(p.lick.crit_bias,'LineStyle',':','LineWidth',2,'Color',p.col.photostim);
yline(-p.lick.crit_bias,'LineStyle',':','LineWidth',2,'Color',p.col.photostim);
ylim([-1,1])
ylabel('Response bias index')
title('Response bias')


%% Response bias as a function of trial block

subplot(2,5,7)
hold on
yline(0,'LineStyle',':','LineWidth',2,'Color',p.col.black);
plot(smoothdata(lick.rbias_blocks.general_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineStyle','--','LineWidth',1,'Color',p.col.darkGray)
plot(smoothdata(lick.rbias_blocks.odourA_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineStyle','--','LineWidth',1,'Color',p.col.AB)
plot(smoothdata(lick.rbias_blocks.odourB_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineStyle','--','LineWidth',1,'Color',p.col.XB)
plot(smoothdata(lick.rbias_blocks.general_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',1,'Color',p.col.darkGray)
plot(smoothdata(lick.rbias_blocks.odourA_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',1,'Color',p.col.AB)
plot(smoothdata(lick.rbias_blocks.odourB_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',1,'Color',p.col.XB)
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
v = bar(these_labels,[lick.lprob.AB_stim,lick.lprob.AB_catch;lick.lprob.XY_stim,lick.lprob.XY_catch;lick.lprob.AY_stim,lick.lprob.AY_catch;lick.lprob.XB_stim,lick.lprob.XB_catch]*100);
v(1).FaceColor = 'flat'; v(2).FaceColor = 'flat'; set(gca,'box','off');
v(1).CData(1,:) = p.col.AB; v(2).CData(1,:) = mean([p.col.AB;p.col.gray],1); 
v(1).CData(2,:) = p.col.XY; v(2).CData(2,:) = mean([p.col.XY;p.col.gray],1); 
v(1).CData(3,:) = p.col.AY; v(2).CData(3,:) = mean([p.col.AY;p.col.gray],1); 
v(1).CData(4,:) = p.col.XB; v(2).CData(4,:) = mean([p.col.XB;p.col.gray],1); 
ylim([0,100])
ytickformat('percentage');
ylabel('P(licks after 2nd odour onset)')
title('Lick probability')


%% Lick probability as a function of trial block

subplot(2,5,8)
hold on
plot(smoothdata(lick.lprob_blocks.AB_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineStyle','--','LineWidth',1,'Color',p.col.AB)
plot(smoothdata(lick.lprob_blocks.XY_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineStyle','--','LineWidth',1,'Color',p.col.XY)
plot(smoothdata(lick.lprob_blocks.AY_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineStyle','--','LineWidth',1,'Color',p.col.AY)
plot(smoothdata(lick.lprob_blocks.XB_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineStyle','--','LineWidth',1,'Color',p.col.XB)
plot(smoothdata(lick.lprob_blocks.AB_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',1,'Color',p.col.AB)
plot(smoothdata(lick.lprob_blocks.XY_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',1,'Color',p.col.XY)
plot(smoothdata(lick.lprob_blocks.AY_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',1,'Color',p.col.AY)
plot(smoothdata(lick.lprob_blocks.XB_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan')*100,'LineWidth',1,'Color',p.col.XB)
ylim([0,100])
ytickformat('percentage');
xlabel('Trial block')
ylabel('P(licks after 2nd odour onset)')
title('Lick probability')


%% Lick latency

these_labels = {char('      A-B'),char(''),char(''),char('       X-Y'),char(''),char(''),char('       A-Y'),char(''),char(''),char('       X-B'),char('')};
this_data = nan(info.task.numTrials/4,8);
this_data(1:length(lick.latency.AB_H_stim),1) = lick.latency.AB_H_stim;
this_data(1:length(lick.latency.AB_H_catch),2) = lick.latency.AB_H_catch;
this_data(1:length(lick.latency.XY_H_stim),3) = lick.latency.XY_H_stim;
this_data(1:length(lick.latency.XY_H_catch),4) = lick.latency.XY_H_catch;
this_data(1:length(lick.latency.AY_FA_stim),5) = lick.latency.AY_FA_stim;
this_data(1:length(lick.latency.AY_FA_catch),6) = lick.latency.AY_FA_catch;
this_data(1:length(lick.latency.XB_FA_stim),7) = lick.latency.XB_FA_stim;
this_data(1:length(lick.latency.XB_FA_catch),8) = lick.latency.XB_FA_catch;

subplot(2,5,4)
hold on
yline(info.task.trialStructure.tOdour2,'LineStyle',':','LineWidth',2,'Color',p.col.odour);
yline(info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay,'LineStyle',':','LineWidth',2,'Color',p.col.reward);
yline(info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow,'LineStyle',':','LineWidth',2,'Color',p.col.reward);
v = violinplot([this_data(:,1),this_data(:,2),nan(info.task.numTrials/4,1),this_data(:,3),this_data(:,4),nan(info.task.numTrials/4,1),this_data(:,5),this_data(:,6),nan(info.task.numTrials/4,1),this_data(:,7),this_data(:,8)],these_labels);
v(1).ViolinColor = p.col.AB; v(2).ViolinColor = mean([p.col.AB;p.col.gray],1); v(4).ViolinColor = p.col.XY; v(5).ViolinColor = mean([p.col.XY;p.col.gray],1); v(7).ViolinColor = p.col.AY; v(8).ViolinColor = mean([p.col.AY;p.col.gray],1); v(10).ViolinColor = p.col.XB; v(11).ViolinColor = mean([p.col.XB;p.col.gray],1);
v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(4).BoxColor = 'k'; v(5).BoxColor = 'k'; v(7).BoxColor = 'k'; v(8).BoxColor = 'k'; v(10).BoxColor = 'k'; v(11).BoxColor = 'k';
xlim([0,length(v)+1])
ylim([0,info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow])
ylabel('Time since 2nd odour onset (s)')
title('Lick latency')


%% Lick latency as a function of trial block

subplot(2,5,9)
hold on
yline(info.task.trialStructure.tOdour2,'LineStyle',':','LineWidth',2,'Color',p.col.odour);
yline(info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay,'LineStyle',':','LineWidth',2,'Color',p.col.reward);
yline(info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow,'LineStyle',':','LineWidth',2,'Color',p.col.reward);
plot(smoothdata(lick.latency_blocks.AB_H_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineStyle','--','LineWidth',1,'Color',p.col.AB)
plot(smoothdata(lick.latency_blocks.XY_H_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineStyle','--','LineWidth',1,'Color',p.col.XY)
plot(smoothdata(lick.latency_blocks.AY_FA_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineStyle','--','LineWidth',1,'Color',p.col.AY)
plot(smoothdata(lick.latency_blocks.XB_FA_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineStyle','--','LineWidth',1,'Color',p.col.XB)
plot(smoothdata(lick.latency_blocks.AB_H_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',1,'Color',p.col.AB)
plot(smoothdata(lick.latency_blocks.XY_H_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',1,'Color',p.col.XY)
plot(smoothdata(lick.latency_blocks.AY_FA_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',1,'Color',p.col.AY)
plot(smoothdata(lick.latency_blocks.XB_FA_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',1,'Color',p.col.XB)
ylim([0,info.task.trialStructure.tOdour2+info.task.trialStructure.tRespDelay+info.task.trialStructure.tRespWindow])
xlabel('Trial block')
ylabel('Time since 2nd odour onset (s)')
title('Lick latency')


%% Lick rigour

these_labels = {char('      A-B'),char(''),char(''),char('       X-Y'),char(''),char(''),char('       A-Y'),char(''),char(''),char('       X-B'),char('')};
this_data = nan(info.task.numTrials/4,8);
this_data(1:length(lick.rigour.AB_H_stim),1) = lick.rigour.AB_H_stim;
this_data(1:length(lick.rigour.AB_H_catch),2) = lick.rigour.AB_H_catch;
this_data(1:length(lick.rigour.XY_H_stim),3) = lick.rigour.XY_H_stim;
this_data(1:length(lick.rigour.XY_H_catch),4) = lick.rigour.XY_H_catch;
this_data(1:length(lick.rigour.AY_FA_stim),5) = lick.rigour.AY_FA_stim;
this_data(1:length(lick.rigour.AY_FA_catch),6) = lick.rigour.AY_FA_catch;
this_data(1:length(lick.rigour.XB_FA_stim),7) = lick.rigour.XB_FA_stim;
this_data(1:length(lick.rigour.XB_FA_catch),8) = lick.rigour.XB_FA_catch;

subplot(2,5,5)
v = violinplot([this_data(:,1),this_data(:,2),nan(info.task.numTrials/4,1),this_data(:,3),this_data(:,4),nan(info.task.numTrials/4,1),this_data(:,5),this_data(:,6),nan(info.task.numTrials/4,1),this_data(:,7),this_data(:,8)],these_labels);
v(1).ViolinColor = p.col.AB; v(2).ViolinColor = mean([p.col.AB;p.col.gray],1); v(4).ViolinColor = p.col.XY; v(5).ViolinColor = mean([p.col.XY;p.col.gray],1); v(7).ViolinColor = p.col.AY; v(8).ViolinColor = mean([p.col.AY;p.col.gray],1); v(10).ViolinColor = p.col.XB; v(11).ViolinColor = mean([p.col.XB;p.col.gray],1);
v(1).BoxColor = 'k'; v(2).BoxColor = 'k'; v(4).BoxColor = 'k'; v(5).BoxColor = 'k'; v(7).BoxColor = 'k'; v(8).BoxColor = 'k'; v(10).BoxColor = 'k'; v(11).BoxColor = 'k';
xlim([0,length(v)+1])
ylim([0,inf])
ylabel(['Number of licks from 2nd odour onset',newline,'(until reward delivery)'])
%ylabel(['Number of licks since 2nd odour onset'])
title('Lick rigour')


%% Lick rigour as a function of trial block

subplot(2,5,10)
hold on
plot(smoothdata(lick.rigour_blocks.AB_H_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineStyle','--','LineWidth',1,'Color',p.col.AB)
plot(smoothdata(lick.rigour_blocks.XY_H_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineStyle','--','LineWidth',1,'Color',p.col.XY)
plot(smoothdata(lick.rigour_blocks.AY_FA_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineStyle','--','LineWidth',1,'Color',p.col.AY)
plot(smoothdata(lick.rigour_blocks.XB_FA_catch,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineStyle','--','LineWidth',1,'Color',p.col.XB)
plot(smoothdata(lick.rigour_blocks.AB_H_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',1,'Color',p.col.AB)
plot(smoothdata(lick.rigour_blocks.XY_H_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',1,'Color',p.col.XY)
plot(smoothdata(lick.rigour_blocks.AY_FA_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',1,'Color',p.col.AY)
plot(smoothdata(lick.rigour_blocks.XB_FA_stim,2,'gaussian',p.lick.smoothing_sd*5,'includenan'),'LineWidth',1,'Color',p.col.XB)
ylim([0,inf])
xlabel('Trial block')
ylabel(['Number of licks from 2nd odour onset',newline,'(until reward delivery)'])
%ylabel(['Number of licks since 2nd odour onset'])
title('Lick rigour')


%% Return

suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', stim vs catch'])

drawnow;
end

