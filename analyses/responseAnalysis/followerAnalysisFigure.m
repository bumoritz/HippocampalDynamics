function F = followerAnalysisFigure(s2p_meta,iscell,trg,resp,flw,p,info)

F = default_figure();


%% Photostimulation followers

subplot(2,4,1)

binsize = 10;
hold on

temp = flw.stats_sig_pos(:);
temp2 = discretize(trg.dist_closestLaser(:),[0:binsize:1000]);
responseProbability = nan(nanmax(temp2),1);
for j=1:nanmax(temp2)
    responseProbability(j) = nanmean(temp(find(temp2==j)));
end
plot(responseProbability*100,'LineWidth',2,'Color','r')

temp = flw.stats_sig_neg(:);
temp2 = discretize(trg.dist_closestLaser(:),[0:binsize:1000]);
responseProbability = nan(nanmax(temp2),1);
for j=1:nanmax(temp2)
    responseProbability(j) = nanmean(temp(find(temp2==j)));
end
plot(responseProbability*100,'LineWidth',2,'Color','b')

% xticks([0:(length([0:binsize:1000]))/10:length([0:binsize:1000])]+1)
% xticklabels({'0','100','200','300','400','500','600','700','800','900','1000'})
% xlim([1,length([0:binsize:400])+1])
ytickformat('percentage');
xlabel('Distance from closest laser beamlet (\mum)\newline(25\mum bin size)')
ylabel('Response probability (stats only)')
legend('positive','negative')
title('Photostimulation followers')


%% Amplitude

nrows = 3; ncols = 3; m=0;
default_figure();

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgTrgAmp(:),flw.ampByResp.avgPosAmp(:),'.','MarkerEdgeColor','r');
[corr_r,corr_p] = fitLine(flw.ampByResp.avgTrgAmp(:),flw.ampByResp.avgPosAmp(:),'k','Pearson');
xlabel('Amplitude of targeted cells')
ylabel('Amplitude of excited followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgTrgAmp(:),flw.ampByResp.avgNeuAmp(:),'.','MarkerEdgeColor',p.col.gray)
[corr_r,corr_p] = fitLine(flw.ampByResp.avgTrgAmp(:),flw.ampByResp.avgNeuAmp(:),'k','Pearson');
xlabel('Amplitude of targeted cells')
ylabel('Amplitude of non-followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgTrgAmp(:),flw.ampByResp.avgNegAmp(:),'.','MarkerEdgeColor','b')
[corr_r,corr_p] = fitLine(flw.ampByResp.avgTrgAmp(:),flw.ampByResp.avgNegAmp(:),'k','Pearson');
xlabel('Amplitude of targeted cells')
ylabel('Amplitude of inhibited followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgPosAmp(:),'.','MarkerEdgeColor','r');
[corr_r,corr_p] = fitLine(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgPosAmp(:),'k','Pearson');
xlabel('Amplitude of responsive cells')
ylabel('Amplitude of excited followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgNeuAmp(:),'.','MarkerEdgeColor',p.col.gray)
[corr_r,corr_p] = fitLine(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgNeuAmp(:),'k','Pearson');
xlabel('Amplitude of responsive cells')
ylabel('Amplitude of non-followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgNegAmp(:),'.','MarkerEdgeColor','b')
[corr_r,corr_p] = fitLine(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgNegAmp(:),'k','Pearson');
xlabel('Amplitude of responsive cells')
ylabel('Amplitude of inhibited followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgRespTrgAmp(:),flw.ampByResp.avgPosAmp(:),'.','MarkerEdgeColor','r');
[corr_r,corr_p] = fitLine(flw.ampByResp.avgRespTrgAmp(:),flw.ampByResp.avgPosAmp(:),'k','Pearson');
xlabel('Amplitude of responsive targeted cells')
ylabel('Amplitude of excited followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgRespTrgAmp(:),flw.ampByResp.avgNeuAmp(:),'.','MarkerEdgeColor',p.col.gray)
[corr_r,corr_p] = fitLine(flw.ampByResp.avgRespTrgAmp(:),flw.ampByResp.avgNeuAmp(:),'k','Pearson');
xlabel('Amplitude of responsive targeted cells')
ylabel('Amplitude of non-followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgRespTrgAmp(:),flw.ampByResp.avgNegAmp(:),'.','MarkerEdgeColor','b')
[corr_r,corr_p] = fitLine(flw.ampByResp.avgRespTrgAmp(:),flw.ampByResp.avgNegAmp(:),'k','Pearson');
xlabel('Amplitude of responsive targeted cells')
ylabel('Amplitude of inhibited followers')

suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', each dot is a trial x cluster combination'])


%% Amplitude

nrows = 3; ncols = 3; m=0;
default_figure();

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgTrgAmp(:),flw.ampByResp.avgPosAmp(:),'.','MarkerEdgeColor','r');
[corr_r,corr_p] = fitLine(flw.ampByResp.avgTrgAmp(:),flw.ampByResp.avgPosAmp(:),'k','Pearson');
plot(p.resp.ampBins_x,flw.ampByResp.avgPosAmp_by_avgTrgAmpBin,'g-','LineWidth',2)
xlabel('Amplitude of targeted cells')
ylabel('Amplitude of excited followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgTrgAmp(:),flw.ampByResp.avgNeuAmp(:),'.','MarkerEdgeColor',p.col.gray)
[corr_r,corr_p] = fitLine(flw.ampByResp.avgTrgAmp(:),flw.ampByResp.avgNeuAmp(:),'k','Pearson');
plot(p.resp.ampBins_x,flw.ampByResp.avgNeuAmp_by_avgTrgAmpBin,'g-','LineWidth',2)
xlabel('Amplitude of targeted cells')
ylabel('Amplitude of non-followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgTrgAmp(:),flw.ampByResp.avgNegAmp(:),'.','MarkerEdgeColor','b')
[corr_r,corr_p] = fitLine(flw.ampByResp.avgTrgAmp(:),flw.ampByResp.avgNegAmp(:),'k','Pearson');
plot(p.resp.ampBins_x,flw.ampByResp.avgNegAmp_by_avgTrgAmpBin,'g-','LineWidth',2)
xlabel('Amplitude of targeted cells')
ylabel('Amplitude of inhibited followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgPosAmp(:),'.','MarkerEdgeColor','r');
[corr_r,corr_p] = fitLine(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgPosAmp(:),'k','Pearson');
plot(p.resp.ampBins_x,flw.ampByResp.avgPosAmp_by_avgRespAmpBin,'g-','LineWidth',2)
xlabel('Amplitude of responsive cells')
ylabel('Amplitude of excited followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgNeuAmp(:),'.','MarkerEdgeColor',p.col.gray)
[corr_r,corr_p] = fitLine(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgNeuAmp(:),'k','Pearson');
plot(p.resp.ampBins_x,flw.ampByResp.avgNeuAmp_by_avgRespAmpBin,'g-','LineWidth',2)
xlabel('Amplitude of responsive cells')
ylabel('Amplitude of non-followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgNegAmp(:),'.','MarkerEdgeColor','b')
[corr_r,corr_p] = fitLine(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgNegAmp(:),'k','Pearson');
plot(p.resp.ampBins_x,flw.ampByResp.avgNegAmp_by_avgRespAmpBin,'g-','LineWidth',2)
xlabel('Amplitude of responsive cells')
ylabel('Amplitude of inhibited followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgRespTrgAmp(:),flw.ampByResp.avgPosAmp(:),'.','MarkerEdgeColor','r');
[corr_r,corr_p] = fitLine(flw.ampByResp.avgRespTrgAmp(:),flw.ampByResp.avgPosAmp(:),'k','Pearson');
plot(p.resp.ampBins_x,flw.ampByResp.avgPosAmp_by_avgRespTrgAmpBin,'g-','LineWidth',2)
xlabel('Amplitude of responsive targeted cells')
ylabel('Amplitude of excited followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgRespTrgAmp(:),flw.ampByResp.avgNeuAmp(:),'.','MarkerEdgeColor',p.col.gray)
[corr_r,corr_p] = fitLine(flw.ampByResp.avgRespTrgAmp(:),flw.ampByResp.avgNeuAmp(:),'k','Pearson');
plot(p.resp.ampBins_x,flw.ampByResp.avgNeuAmp_by_avgRespTrgAmpBin,'g-','LineWidth',2)
xlabel('Amplitude of responsive targeted cells')
ylabel('Amplitude of non-followers')

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgRespTrgAmp(:),flw.ampByResp.avgNegAmp(:),'.','MarkerEdgeColor','b')
[corr_r,corr_p] = fitLine(flw.ampByResp.avgRespTrgAmp(:),flw.ampByResp.avgNegAmp(:),'k','Pearson');
plot(p.resp.ampBins_x,flw.ampByResp.avgNegAmp_by_avgRespTrgAmpBin,'g-','LineWidth',2)
xlabel('Amplitude of responsive targeted cells')
ylabel('Amplitude of inhibited followers')

suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', each dot is a trial x cluster combination'])


%% Amplitude

nrows = 1; ncols = 3; m=0;
default_figure();

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgPosAmp(:),'.','MarkerEdgeColor','r');
[corr_r,corr_p] = fitLine(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgPosAmp(:),'k','Pearson');
%plot(p.resp.ampBins_x,flw.ampByResp.avgPosAmp_by_avgRespAmpBin,'g-','LineWidth',2)
xlabel('Amplitude of responsive cells')
ylabel('Amplitude of excited followers')
title(['rho=',num2str(corr_r,2),', p=',num2str(corr_p,2)])

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgNeuAmp(:),'.','MarkerEdgeColor',p.col.gray)
[corr_r,corr_p] = fitLine(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgNeuAmp(:),'k','Pearson');
%plot(p.resp.ampBins_x,flw.ampByResp.avgNeuAmp_by_avgRespAmpBin,'g-','LineWidth',2)
xlabel('Amplitude of responsive cells')
ylabel('Amplitude of non-followers')
title(['rho=',num2str(corr_r,2),', p=',num2str(corr_p,2)])

m = m+1; subplot(nrows,ncols,m); hold on;
scatter(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgNegAmp(:),'.','MarkerEdgeColor','b')
[corr_r,corr_p] = fitLine(flw.ampByResp.avgRespAmp(:),flw.ampByResp.avgNegAmp(:),'k','Pearson');
%plot(p.resp.ampBins_x,flw.ampByResp.avgNegAmp_by_avgRespAmpBin,'g-','LineWidth',2)
xlabel('Amplitude of responsive cells')
ylabel('Amplitude of inhibited followers')
title(['rho=',num2str(corr_r,2),', p=',num2str(corr_p,2)])

suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', each dot is a trial x cluster combination'])


%% Return

suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType])
drawnow;
end
