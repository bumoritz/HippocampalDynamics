function behConfoundAnalysis(info,ops,p,path,task,paq_beh)

%% Split running data into trials

[prop,~,~] = preprocessActivityMeasure([],p.bcon,p,paq_beh,[],true);
[events_binned] = binPaqEvents(paq_beh,task,p,prop);
trials_all = createTrialsStruct(task,1:prop.numTrials_incl);

try
    licking = reshape(events_binned.lick(prop.trial_frames_binned'),size(prop.trial_frames_binned,2),size(prop.trial_frames_binned,1));
    distance = reshape(events_binned.distance(prop.trial_frames_binned'),size(prop.trial_frames_binned,2),size(prop.trial_frames_binned,1));
    velocity = reshape(events_binned.velocity(prop.trial_frames_binned'),size(prop.trial_frames_binned,2),size(prop.trial_frames_binned,1));
    acceleration = reshape(events_binned.acceleration(prop.trial_frames_binned'),size(prop.trial_frames_binned,2),size(prop.trial_frames_binned,1));
catch
    licking = reshape(events_binned.lick(prop.trial_frames_binned(:,1:end-1)'),size(prop.trial_frames_binned,2)-1,size(prop.trial_frames_binned,1));
    distance = reshape(events_binned.distance(prop.trial_frames_binned(:,1:end-1)'),size(prop.trial_frames_binned,2)-1,size(prop.trial_frames_binned,1));
    velocity = reshape(events_binned.velocity(prop.trial_frames_binned(:,1:end-1)'),size(prop.trial_frames_binned,2)-1,size(prop.trial_frames_binned,1));
    acceleration = reshape(events_binned.acceleration(prop.trial_frames_binned(:,1:end-1)'),size(prop.trial_frames_binned,2)-1,size(prop.trial_frames_binned,1));
end
% odourA = reshape(events_binned.odourA(prop.trial_frames_binned'),size(prop.trial_frames_binned,2),size(prop.trial_frames_binned,1));


%% Create bcon struct

bcon.binwise_full.licking = licking;
bcon.binwise_AW.licking = licking(:,p.general.bins_analysisWindow);
bcon.trialwise_AW.licking = nanmean(bcon.binwise_AW.licking,2);

bcon.binwise_full.distance = distance;
bcon.binwise_full.velocity = velocity;
bcon.binwise_full.acceleration = acceleration;

bcon.binwise_AW.distance = distance(:,p.general.bins_analysisWindow);
bcon.binwise_AW.velocity = velocity(:,p.general.bins_analysisWindow);
bcon.binwise_AW.acceleration = acceleration(:,p.general.bins_analysisWindow);

bcon.trialwise_full.distance = nanmean(bcon.binwise_full.distance,2);
bcon.trialwise_full.velocity = nanmean(bcon.binwise_full.velocity,2);
bcon.trialwise_full.acceleration = nanmean(bcon.binwise_full.acceleration,2);

bcon.trialwise_AW.distance = nanmean(bcon.binwise_AW.distance,2);
bcon.trialwise_AW.velocity = nanmean(bcon.binwise_AW.velocity,2);
bcon.trialwise_AW.acceleration = nanmean(bcon.binwise_AW.acceleration,2);

% fit line to velocity in AW
bcon.binwise_AW.velocity_lin_m = nan(size(bcon.binwise_AW.velocity,1),1);
bcon.binwise_AW.velocity_lin_b = nan(size(bcon.binwise_AW.velocity,1),1);
for i=1:size(bcon.binwise_AW.velocity,1)
    temp = fit(p.general.t_binned(p.general.bins_analysisWindow)',bcon.binwise_AW.velocity(i,:)','poly1');
    bcon.binwise_AW.velocity_lin_m(i) = temp.p1;
    bcon.binwise_AW.velocity_lin_b(i) = temp.p2;
end

% blockwise
bcon.blockwise.licking_AW = nan(size(velocity,1)/100,1);
bcon.blockwise.velocity_full = nan(size(velocity,1)/100,1);
bcon.blockwise.distance_AW = nan(size(velocity,1)/100,1);
bcon.blockwise.velocity_AW = nan(size(velocity,1)/100,1);
bcon.blockwise.velocity_AW_m = nan(size(velocity,1)/100,1);
bcon.blockwise.velocity_AW_b = nan(size(velocity,1)/100,1);
bcon.blockwise.acceleration_AW = nan(size(velocity,1)/100,1);
for i=1:size(velocity,1)/100
    bcon.blockwise.licking_AW(i) = nanmean(bcon.trialwise_AW.licking((i-1)*100+1:i*100));
    bcon.blockwise.velocity_full(i) = nanmean(bcon.trialwise_full.velocity((i-1)*100+1:i*100));
    bcon.blockwise.velocity_AW(i) = nanmean(bcon.trialwise_AW.velocity((i-1)*100+1:i*100));
    bcon.blockwise.velocity_AW_m(i) = nanmean(bcon.binwise_AW.velocity_lin_m((i-1)*100+1:i*100));
    bcon.blockwise.velocity_AW_b(i) = nanmean(bcon.binwise_AW.velocity_lin_b((i-1)*100+1:i*100));
    bcon.blockwise.distance_AW(i) = nanmean(bcon.trialwise_AW.distance((i-1)*100+1:i*100));
    bcon.blockwise.acceleration_AW(i) = nanmean(bcon.trialwise_AW.acceleration((i-1)*100+1:i*100));
end

% save
bcon = orderfields(bcon);
save([path.filepart_out,'bcon.mat'],'bcon','-v7.3');
disp(['--- Saved bcon file as ',[path.filepart_out,'bcon.mat'],'.'])


%% Trial-by-trial figure

% nrows = 2; ncols = 5; m=0;
% F = default_figure();
% plt = struct(); plt.illustrator = false; plt.colormap = 'copper'; plt.ylim = [1,length(trials_all.stimuli.AB)]; plt.ylabel = 'Trial'; plt.prop = prop; plt.clim = [0,200]; %[nanmin(velocity(:)),nanmax(velocity(:))]; % flipud(eval('copper'))
% 
% m = m+1; subplot(nrows,ncols,m); hold on;
% heatMap_task(velocity(trials_all.stimuli.AB,:),[],plt,p,info,plt); xlim([0,p.general.numBins]); 
% title('AB')
% 
% m = m+1; subplot(nrows,ncols,m); hold on;
% heatMap_task(velocity(trials_all.stimuli.AY,:),[],plt,p,info,plt); xlim([0,p.general.numBins])
% title('AY')
% 
% m = m+1; subplot(nrows,ncols,m); hold on;
% heatMap_task(velocity(trials_all.stimuli.XY,:),[],plt,p,info,plt); xlim([0,p.general.numBins])
% title('XY')
% 
% m = m+1; subplot(nrows,ncols,m); hold on;
% heatMap_task(velocity(trials_all.stimuli.XB,:),[],plt,p,info,plt); xlim([0,p.general.numBins])
% title('XB')
% 
% m = m+1; subplot(nrows,ncols,m); hold on;
% colormap(gca,plt.colormap); set(gca,'CLim',plt.clim); h=colorbar; h.Label.String = 'Velocity (cm/s)'; set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]); set(gca,'box','off');
% 
% plt.colormap = redblue; plt.clim = [-100,100]; %[-nanmax(abs(acceleration(:))),nanmax(abs(acceleration(:)))];
% 
% m = m+1; subplot(nrows,ncols,m); hold on;
% heatMap_task(acceleration(trials_all.stimuli.AB,:),[],plt,p,info,plt); xlim([0,p.general.numBins])
% title('AB')
% 
% m = m+1; subplot(nrows,ncols,m); hold on;
% heatMap_task(acceleration(trials_all.stimuli.AY,:),[],plt,p,info,plt); xlim([0,p.general.numBins])
% title('AY')
% 
% m = m+1; subplot(nrows,ncols,m); hold on;
% heatMap_task(acceleration(trials_all.stimuli.XY,:),[],plt,p,info,plt); xlim([0,p.general.numBins])
% title('XY')
% 
% m = m+1; subplot(nrows,ncols,m); hold on;
% heatMap_task(acceleration(trials_all.stimuli.XB,:),[],plt,p,info,plt); xlim([0,p.general.numBins])
% title('XB')
% 
% m = m+1; subplot(nrows,ncols,m); hold on;
% colormap(gca,plt.colormap); set(gca,'CLim',plt.clim); h=colorbar; h.Label.String = 'Acceleration (cm/s^2)'; set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]); set(gca,'box','off');
% 
% if info.stimSession
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType])
% else
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-','nostim']);
% end
% 
% %F = lickingAnalysisFigure(lick,perf,p,info);
% savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','runningBehaviour.fig']);
% saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','runningBehaviour.png']);
% disp(['--- Saved running behaviour figure to ',path.filepart_outX,'plots.'])
% drawnow;
% 
% 
% %% Summary figure
% 
% nrows = 1; ncols = 2; m=0;
% F = default_figure([-20,0.5,10,5]);
% 
% m = m+1; subplot(nrows,ncols,m); hold on;
% shadedErrorBar(1:p.general.numBins,nanmean(velocity(trials_all.stimuli.XB,:),1),nansem(velocity(trials_all.stimuli.XB,:),1),'lineProps',p.col.XB)
% shadedErrorBar(1:p.general.numBins,nanmean(velocity(trials_all.stimuli.AY,:),1),nansem(velocity(trials_all.stimuli.AY,:),1),'lineProps',p.col.AY)
% shadedErrorBar(1:p.general.numBins,nanmean(velocity(trials_all.stimuli.XY,:),1),nansem(velocity(trials_all.stimuli.XY,:),1),'lineProps',p.col.XY)
% shadedErrorBar(1:p.general.numBins,nanmean(velocity(trials_all.stimuli.AB,:),1),nansem(velocity(trials_all.stimuli.AB,:),1),'lineProps',p.col.AB)
% % shadedErrorBar(1:p.general.numBins,nanmean(velocity,1),nansem(velocity,1),'lineProps',p.col.black)
% taskLines(p,info);
% ylim([0,100])
% ylabel('Velocity (cm/s)')
% 
% m = m+1; subplot(nrows,ncols,m); hold on;
% shadedErrorBar(1:p.general.numBins,nanmean(acceleration(trials_all.stimuli.XB,:),1),nansem(acceleration(trials_all.stimuli.XB,:),1),'lineProps',p.col.XB)
% shadedErrorBar(1:p.general.numBins,nanmean(acceleration(trials_all.stimuli.AY,:),1),nansem(acceleration(trials_all.stimuli.AY,:),1),'lineProps',p.col.AY)
% shadedErrorBar(1:p.general.numBins,nanmean(acceleration(trials_all.stimuli.XY,:),1),nansem(acceleration(trials_all.stimuli.XY,:),1),'lineProps',p.col.XY)
% shadedErrorBar(1:p.general.numBins,nanmean(acceleration(trials_all.stimuli.AB,:),1),nansem(acceleration(trials_all.stimuli.AB,:),1),'lineProps',p.col.AB)
% % shadedErrorBar(1:p.general.numBins,nanmean(velocity,1),nansem(velocity,1),'lineProps',p.col.black)
% taskLines(p,info);
% ylim([-100,100])
% ylabel('Acceleration (cm/s^2)')
% 
% 
% if info.stimSession
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType])
% else
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-','nostim']);
% end
% savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','runningBehaviour_cmpr.fig']);
% saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','runningBehaviour_cmpr.png']);
% disp(['--- Saved running behaviour figure to ',path.filepart_outX,'plots.'])
% drawnow;
% 
% 
% %% Summary figure - 100t blocks
% 
% nrows = 1; ncols = 3; m=0;
% F = default_figure([-20,0.5,15,5]);
% these_rgbs = discretisedColourMap('winter',false,size(velocity,1)/100);
% 
% m = m+1; subplot(nrows,ncols,m); hold on;
% for i=1:size(velocity,1)/100
%     shadedErrorBar(1:p.general.numBins,nanmean(velocity((i-1)*100+1:i*100,:),1),nansem(velocity((i-1)*100+1:i*100,:),1),'lineProps',these_rgbs(i,:))
% end
% taskLines(p,info);
% ylim([0,120])
% ylabel('Velocity (cm/s)')
% 
% m = m+1; subplot(nrows,ncols,m); hold on;
% for i=1:size(acceleration,1)/100
%     shadedErrorBar(1:p.general.numBins,nanmean(acceleration((i-1)*100+1:i*100,:),1),nansem(acceleration((i-1)*100+1:i*100,:),1),'lineProps',these_rgbs(i,:))
% end
% taskLines(p,info);
% ylim([-100,100])
% ylabel('Acceleration (cm/s^2)')
% 
% m = m+1; subplot(nrows,ncols,m); hold on;
% colormap(gca,winter); set(gca,'CLim',[1,size(velocity,1)/100]); h=colorbar; h.Ticks=[1:size(velocity,1)/100]; h.Label.String = 'Trial block (100 trials)'; set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]); set(gca,'box','off');
% 
% if info.stimSession
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType])
% else
%     suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-','nostim']);
% end
% savefig(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','runningBehaviour_cmpr_100t.fig']);
% saveas(F,[path.filepart_outX,'plots\',info.animal,'_',info.date,'_','runningBehaviour_cmpr_100t.png']);
% disp(['--- Saved running behaviour figure to ',path.filepart_outX,'plots.'])
% drawnow;
% 
% 
% %% Summary figure - 100t blocks
% 
% F = default_figure([-20,0.5,5,5]);
% these_rgbs = discretisedColourMap('winter',false,size(velocity,1)/100);
% 
% hold on;
% for i=1:size(velocity,1)/100
%     shadedErrorBar(1:p.general.numBins,nanmean(licking((i-1)*100+1:i*100,:),1),nansem(licking((i-1)*100+1:i*100,:),1),'lineProps',these_rgbs(i,:))
% end
% taskLines(p,info);
% %ylim([0,120])
% ylabel('Lick probability')
% 
% F = default_figure([-20,0.5,5,5]);
% these_rgbs = discretisedColourMap('winter',false,size(velocity,1)/100);
% 
% hold on;
% for i=1:size(velocity,1)/100
%     shadedErrorBar(1:p.general.numBins,nanmean(velocity((i-1)*100+1:i*100,:),1),nansem(velocity((i-1)*100+1:i*100,:),1),'lineProps',these_rgbs(i,:))
% end
% taskLines(p,info);
% ylim([0,120])
% ylabel('Velocity (cm/s)')
% 
% F = default_figure([-20,0.5,5,5]);
% these_rgbs = discretisedColourMap('winter',false,size(velocity,1)/100);
% 
% hold on;
% yline(0,':')
% for i=1:size(velocity,1)/100
%     shadedErrorBar(1:p.general.numBins,nanmean(velocity((i-1)*100+1:i*100,:),1)-nanmean(velocity((1-1)*100+1:1*100,:),1),nansem(velocity((i-1)*100+1:i*100,:),1),'lineProps',these_rgbs(i,:))
% end
% taskLines(p,info);
% ylim([-80,80])
% ylabel('Delta velocity (cm/s)')
% 
% 
% F = default_figure([-20,0.5,5,5]);
% these_rgbs = discretisedColourMap('winter',false,size(velocity,1)/100);
% 
% hold on;
% yline(0,':')
% for i=1:size(velocity,1)/100
%     shadedErrorBar(1:p.general.numBins,nanmean(velocity((i-1)*100+1:i*100,:),1)-nanmean(velocity,1),nansem(velocity((i-1)*100+1:i*100,:),1),'lineProps',these_rgbs(i,:))
% end
% taskLines(p,info);
% ylim([-80,80])
% ylabel('Deviation from average velocity profile (cm/s)')


%% Return

if ops.close_figures
    close all;
end
end





