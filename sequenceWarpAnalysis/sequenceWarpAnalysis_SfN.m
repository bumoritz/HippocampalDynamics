%% sequenceWarpAnalysis_SfN

warp1 = load('C:\Users\Moritz\Documents\Illustrator\SfN2022\Data\Elrond_20210601_warp.mat');
warp1 = warp1.warp;
perf1 = load('C:\Users\Moritz\Documents\Illustrator\SfN2022\Data\Elrond_20210601_perf.mat');
perf1 = perf1.perf;
paq_beh1 = load('C:\Users\Moritz\Documents\Illustrator\SfN2022\Data\Elrond_20210601_paq_beh.mat');
paq_beh1 = paq_beh1.paq_beh;
bcon1 = load('C:\Users\Moritz\Documents\Illustrator\SfN2022\Data\Elrond_20210601_bcon.mat');
bcon1 = bcon1.bcon;

warp2 = load('C:\Users\Moritz\Documents\Illustrator\SfN2022\Data\Elrond_20210602_warp.mat');
warp2 = warp2.warp;
perf2 = load('C:\Users\Moritz\Documents\Illustrator\SfN2022\Data\Elrond_20210602_perf.mat');
perf2 = perf2.perf;
paq_beh2 = load('C:\Users\Moritz\Documents\Illustrator\SfN2022\Data\Elrond_20210602_paq_beh.mat');
paq_beh2 = paq_beh2.paq_beh;
bcon2 = load('C:\Users\Moritz\Documents\Illustrator\SfN2022\Data\Elrond_20210602_bcon.mat');
bcon2 = bcon2.bcon;


%%

nrows = 4; ncols = 1; m=0;
F = default_figure([-20,0.5,5,5]); %default_figure([-20,0.5,18,1.5]);

m = m+1; subplot(nrows,ncols,m); hold on;
xline(8.5,'k:')
yline(50,'k:')
temp = [movmean(perf1.blocks_general.correct,5),movmean(perf2.blocks_general.correct,5)];
plot(temp(3:5:end)*100,'k-')
xticks([1:13])
xticklabels({})
ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel('Performance')

m = m+1; subplot(nrows,ncols,m); hold on;
xline(8.5,'k:')
these_edges1 = 1:length(paq_beh1.speed)/9:length(paq_beh1.speed); 
these_edges2 = 1:length(paq_beh2.speed)/6:length(paq_beh2.speed);
temp = nan(1,13);
for k=1:8
    temp(k) = nanmean(paq_beh1.speed(these_edges1(k):these_edges1(k+1)));
end
for k=1:5
    temp(8+k) = nanmean(paq_beh2.speed(these_edges2(k):these_edges2(k+1)));
end 
plot(temp,'k-')
xticks([1:13])
xticklabels({})
ylim([0,100])
yticks([0,50,100])
ylabel(['Mean velocity',newline,'(cm/s)'])

m = m+1; subplot(nrows,ncols,m); hold on;
xline(8.5,'k:')
temp = nan(1,13);
for k=1:8
    temp(k) = warp1{k}.exp1.a;
end
for k=1:5
    temp(8+k) = warp2{k}.exp1.a;
end 
plot(temp,'k-')
xticks([1:13])
xticklabels({})
ylim([0,0.1])
yticks([0,0.1])
ylabel(['Coefficient a',newline,'(intercept)'])

m = m+1; subplot(nrows,ncols,m); hold on;
xline(8.5,'k:')
yline(0,'k:')
temp = nan(1,13);
for k=1:8
    temp(k) = warp1{k}.exp1.b;
end
for k=1:5
    temp(8+k) = warp2{k}.exp1.b;
end 
plot(temp,'k-')
xticks([1:13])
%xticklabels({'1',''}) 
ylim([-0.5,0.5])
yticks([-0.5,0,0.5])
xlabel('Trial block (in blocks of 100 trials)')
ylabel(['Coefficient b',newline,'(decay rate, s^{-1})'])

% savefig(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Examples\Elrond (switch imaging, runner)\summary.fig']);
% saveas(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Examples\Elrond (switch imaging, runner)\summary.png']);


%%

nrows = 4; ncols = 1; m=0;
F = default_figure([-20,0.5,5,5]); %default_figure([-20,0.5,18,1.5]);

m = m+1; subplot(nrows,ncols,m); hold on;
xline(8.5,'k:')
yline(50,'k:')
temp = [movmean(perf1.blocks_general.correct,5),movmean(perf2.blocks_general.correct,5)];
plot(temp(3:5:end)*100,'k-')
xticks([1:13])
xticklabels({})
ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel('Performance')

m = m+1; subplot(nrows,ncols,m); hold on;
xline(8.5,'k:')
plot([bcon1.blockwise.velocity_AW;bcon2.blockwise.velocity_AW],'k-')
xticks([1:13])
xticklabels({})
ylabel(['Velocity in AW'])

m = m+1; subplot(nrows,ncols,m); hold on;
xline(8.5,'k:')
plot([bcon1.blockwise.velocity_AW_m;bcon2.blockwise.velocity_AW_m],'k-')
xticks([1:13])
xticklabels({})
ylabel(['Velocity slope in AW'])

m = m+1; subplot(nrows,ncols,m); hold on;
xline(8.5,'k:')
yline(0,'k:')
temp = nan(1,13);
for k=1:8
    temp(k) = warp1{k}.exp1.b;
end
for k=1:5
    temp(8+k) = warp2{k}.exp1.b;
end 
plot(temp,'k-')
xticks([1:13])
%xticklabels({'1',''}) 
ylim([-0.5,0.5])
yticks([-0.5,0,0.5])
xlabel('Trial block (in blocks of 100 trials)')
ylabel(['Coefficient b',newline,'(decay rate, s^{-1})'])


%%

nrows = 4; ncols = 1; m=0;
F = default_figure([-20,0.5,5,5]); %default_figure([-20,0.5,18,1.5]);

m = m+1; subplot(nrows,ncols,m); hold on;
xline(8.5,'k:')
yline(50,'k:')
temp = [movmean(perf1.blocks_general.correct,5),movmean(perf2.blocks_general.correct,5)];
plot(temp(3:5:end)*100,'k-')
xticks([1:13])
xticklabels({})
ylim([35,100])
yticks([50,100])
ytickformat('percentage')
ylabel('Performance')

m = m+1; subplot(nrows,ncols,m); hold on;
xline(8.5,'k:')
plot([bcon1.blockwise.acceleration_AW;bcon2.blockwise.acceleration_AW],'k-')
xticks([1:13])
xticklabels({})
ylabel(['Acceleration in AW'])

m = m+1; subplot(nrows,ncols,m); hold on;
xline(8.5,'k:')
plot([bcon1.blockwise.licking_AW;bcon2.blockwise.licking_AW],'k-')
xticks([1:13])
xticklabels({})
ylabel(['Licking in AW'])

m = m+1; subplot(nrows,ncols,m); hold on;
xline(8.5,'k:')
yline(0,'k:')
temp = nan(1,13);
for k=1:8
    temp(k) = warp1{k}.exp1.b;
end
for k=1:5
    temp(8+k) = warp2{k}.exp1.b;
end 
plot(temp,'k-')
xticks([1:13])
%xticklabels({'1',''}) 
ylim([-0.5,0.5])
yticks([-0.5,0,0.5])
xlabel('Trial block (in blocks of 100 trials)')
ylabel(['Coefficient b',newline,'(decay rate, s^{-1})'])



