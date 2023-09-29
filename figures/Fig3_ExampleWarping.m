%% Fig3_ExampleWarping

p = get_p; info = get_info;

path.root_summary = 'C:\SniffinHippo\Summary\';
save_root_fig = [path.root_summary,'figures\Fig3_fig\'];
save_root_png = [path.root_summary,'figures\Fig3_png\'];
save_root_pdf = [path.root_summary,'figures\Fig3_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig3_txt\'];

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


%% Fig3_ExampleWarping_Distance

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;
these_rgbs = discretisedColourMap('winter',false,13);

hold on;
taskLines(p,info);
yline(0,'k-');
for i=1:8
    plot(1:p.general.numBins,nanmean(bcon1.binwise_full.distance((i-1)*100+1:i*100,:),1),'Color',these_rgbs(i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(bcon1.binwise_full.distance((i-1)*100+1:i*100,:),1),nansem(bcon1.binwise_full.distance((i-1)*100+1:i*100,:),1),'lineProps',these_rgbs(i,:))
end
for i=1:5
    plot(1:p.general.numBins,nanmean(bcon2.binwise_full.distance((i-1)*100+1:i*100,:),1),'Color',these_rgbs(8+i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(bcon2.binwise_full.distance((i-1)*100+1:i*100,:),1),nansem(bcon2.binwise_full.distance((i-1)*100+1:i*100,:),1),'lineProps',these_rgbs(8+i,:))
end
ylim([-400,800])
yticks([-400:400:800])
ylabel('Distance traveled (cm)')

savefig(F,[save_root_fig,'\Fig3_ExampleWarping_Distance.fig']);
saveas(F,[save_root_png,'\Fig3_ExampleWarping_Distance.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_ExampleWarping_Distance.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_ExampleWarping_Velocity

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;
these_rgbs = discretisedColourMap('winter',false,13);

hold on;
taskLines(p,info);
for i=1:8
    plot(1:p.general.numBins,nanmean(bcon1.binwise_full.velocity((i-1)*100+1:i*100,:),1),'Color',these_rgbs(i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(bcon1.binwise_full.velocity((i-1)*100+1:i*100,:),1),nansem(bcon1.binwise_full.velocity((i-1)*100+1:i*100,:),1),'lineProps',these_rgbs(i,:))
end
for i=1:5
    plot(1:p.general.numBins,nanmean(bcon2.binwise_full.velocity((i-1)*100+1:i*100,:),1),'Color',these_rgbs(8+i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(bcon2.binwise_full.velocity((i-1)*100+1:i*100,:),1),nansem(bcon2.binwise_full.velocity((i-1)*100+1:i*100,:),1),'lineProps',these_rgbs(8+i,:))
end
ylim([0,120])
yticks([0:40:120])
ylabel('Velocity (cm/s)')

savefig(F,[save_root_fig,'\Fig3_ExampleWarping_Velocity.fig']);
saveas(F,[save_root_png,'\Fig3_ExampleWarping_Velocity.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_ExampleWarping_Velocity.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_ExampleWarping_Acceleration

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;
these_rgbs = discretisedColourMap('winter',false,13);

hold on;
taskLines(p,info);
yline(0,'k-');
for i=1:8
    plot(1:p.general.numBins,nanmean(bcon1.binwise_full.acceleration((i-1)*100+1:i*100,:),1),'Color',these_rgbs(i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(bcon1.binwise_full.acceleration((i-1)*100+1:i*100,:),1),nansem(bcon1.binwise_full.acceleration((i-1)*100+1:i*100,:),1),'lineProps',these_rgbs(i,:))
end
for i=1:5
    plot(1:p.general.numBins,nanmean(bcon2.binwise_full.acceleration((i-1)*100+1:i*100,:),1),'Color',these_rgbs(8+i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(bcon2.binwise_full.acceleration((i-1)*100+1:i*100,:),1),nansem(bcon2.binwise_full.acceleration((i-1)*100+1:i*100,:),1),'lineProps',these_rgbs(8+i,:))
end
ylim([-80,80])
yticks([-80:40:80])
ylabel('Acceleration (cm/s2)')

savefig(F,[save_root_fig,'\Fig3_ExampleWarping_Acceleration.fig']);
saveas(F,[save_root_png,'\Fig3_ExampleWarping_Acceleration.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_ExampleWarping_Acceleration.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_ExampleWarping_Licking

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;
these_rgbs = discretisedColourMap('winter',false,13);

hold on;
taskLines(p,info);
for i=1:8
    plot(1:p.general.numBins,nanmean(bcon1.binwise_full.licking((i-1)*100+1:i*100,:),1)*100,'Color',these_rgbs(i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(bcon1.binwise_full.licking((i-1)*100+1:i*100,:),1)*100,nansem(bcon1.binwise_full.licking((i-1)*100+1:i*100,:),1)*100,'lineProps',these_rgbs(i,:))
end
for i=1:5
    plot(1:p.general.numBins,nanmean(bcon2.binwise_full.licking((i-1)*100+1:i*100,:),1)*100,'Color',these_rgbs(8+i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(bcon2.binwise_full.licking((i-1)*100+1:i*100,:),1)*100,nansem(bcon2.binwise_full.licking((i-1)*100+1:i*100,:),1)*100,'lineProps',these_rgbs(8+i,:))
end
ytickformat('percentage')
ylim([0,100])
yticks([0,50,100])
ylabel('Lick probability')

savefig(F,[save_root_fig,'\Fig3_ExampleWarping_Licking.fig']);
saveas(F,[save_root_png,'\Fig3_ExampleWarping_Licking.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_ExampleWarping_Licking.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_ExampleWarping_Colorbar

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
cmap = winter(13);
colormap(cmap);
h=colorbar;

savefig(F,[save_root_fig,'\Fig3_ExampleWarping_Colorbar.fig']);
saveas(F,[save_root_png,'\Fig3_ExampleWarping_Colorbar.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_ExampleWarping_Colorbar.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_ExampleWarping_Warping

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
these_rgbs = discretisedColourMap('winter',false,13);

this_x = 0:0.01:5;
for i=1:8
    this_y = warp1{i}.exp1.a*exp(warp1{i}.exp1.b*this_x);
    h=plot(this_x,this_y*100,'Color',these_rgbs(i,:));
end
for i=1:5
    this_y = warp2{i}.exp1.a*exp(warp2{i}.exp1.b*this_x);
    h=plot(this_x,this_y*100,'Color',these_rgbs(8+i,:));
end
ytickformat('percentage')
xlim([0,5])
xticks([0,5])
xlabel('Firing field peak (s)')
ylim([0,10])
yticks([0,5,10])
ylabel('Proportion of sequence cells')

savefig(F,[save_root_fig,'\Fig3_ExampleWarping_Warping.fig']);
saveas(F,[save_root_png,'\Fig3_ExampleWarping_Warping.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_ExampleWarping_Warping.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_ExampleWarping_BconSummary

nrows = 4; ncols = 1; m=0;
F = paper_figure([0,0.5,mm2inch(34),mm2inch(1*34)]); hold on;

m = m+1; subplot(nrows,ncols,m); hold on;
% xline(8.5,'k-')
temp = [movmean(perf1.blocks_general.correct,5),movmean(perf2.blocks_general.correct,5)];
plot(temp(3:5:end)*100,'k-')
xlim([0,13]); xticks([1:13]); xticklabels({});
%ylim([35,100])
yticks([50,100])
ytickformat('percentage')
%ylabel('Performance')

m = m+1; subplot(nrows,ncols,m); hold on;
% xline(8.5,'k-')
temp = nan(1,13);
for k=1:8
    temp(k) = warp1{k}.exp1.b;
end
for k=1:5
    temp(8+k) = warp2{k}.exp1.b;
end 
plot(temp,'k-')
xlim([0,13]); xticks([1:13]); xticklabels({});
ylim([-0.5,0.5]); yticks([-0.5,0.5]);
%ylabel(['Coefficient'])

m = m+1; subplot(nrows,ncols,m); hold on;
% xline(8.5,'k-')
temp = [bcon1.trialwise_AW.velocity;bcon2.trialwise_AW.velocity];
these_means = nan(13,1); these_sems = nan(13,1);
for i=1:13
    these_means(i) = nanmean(temp((i-1)*100+1:i*100,:));
    these_sems(i) = nansem(temp((i-1)*100+1:i*100,:)');
end
%shadedErrorBar(1:13,these_means,these_sems,'lineProps',p.col.black)
plot(1:13,these_means,'k-')
xlim([0,13]); xticks([1:13]); xticklabels({});
ylim([40,80]); yticks([40,80]);
%ylabel(['Velocity'])

m = m+1; subplot(nrows,ncols,m); hold on;
% xline(8.5,'k-')
temp = [bcon1.trialwise_AW.acceleration;bcon2.trialwise_AW.acceleration];
these_means = nan(13,1); these_sems = nan(13,1);
for i=1:13
    these_means(i) = nanmean(temp((i-1)*100+1:i*100,:));
    these_sems(i) = nansem(temp((i-1)*100+1:i*100,:)');
end
%shadedErrorBar(1:13,these_means,these_sems,'lineProps',p.col.black)
plot(1:13,these_means,'k-')
xlim([0,13]); xticks([1:13]); xticklabels({});
ylim([-12,2]); yticks([-12,2]);
%ylabel(['Acceleration'])

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
%xlabel('Block of 100 trials)')

savefig(F,[save_root_fig,'\Fig3_ExampleWarping_BconSummary.fig']);
saveas(F,[save_root_png,'\Fig3_ExampleWarping_BconSummary.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_ExampleWarping_BconSummary.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_ExampleWarping_BconSummary2

nrows = 2; ncols = 1; m=0;
F = paper_figure([0,0.5,mm2inch(34),mm2inch(1*34)]); hold on;

m = m+1; subplot(nrows,ncols,m); hold on;
% xline(8.5,'k-')
temp = [movmean(perf1.blocks_general.correct,5),movmean(perf2.blocks_general.correct,5)];
plot(temp(3:5:end)*100,'k-')
xlim([0,13]); xticks([1:13]); xticklabels({});
%ylim([35,100])
yticks([50,100])
ytickformat('percentage')
%ylabel('Performance')

m = m+1; subplot(nrows,ncols,m); hold on;
% xline(8.5,'k-')
temp = nan(1,13);
for k=1:8
    temp(k) = warp1{k}.exp1.b;
end
for k=1:5
    temp(8+k) = warp2{k}.exp1.b;
end 
plot(temp,'k-')
xlim([0,13]); xticks([1:13]); xticklabels({});
ylim([-0.5,0.5]); yticks([-0.5,0.5]);
%ylabel(['Coefficient'])

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials)')

savefig(F,[save_root_fig,'\Fig3_ExampleWarping_BconSummary2.fig']);
saveas(F,[save_root_png,'\Fig3_ExampleWarping_BconSummary2.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_ExampleWarping_BconSummary2.pdf']); set(gcf,'Color',[1,1,1])













