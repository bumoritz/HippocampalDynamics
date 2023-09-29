%% Fig2_SpatialAnalysis

% import data using Summary_Master with ops.do_sequenceCellSummary = true;
% d_info = selectDataset(d_info,'-g02345678-d1-e01-r01-p01-l01-i00',sheet,path,ops);
% run all data gathering steps from Fig2_SequenceAnalysis


%% Extract data

% peakAmplitude, meanAmplitude, meanAmplitudeActive

corr_speed_timeAtPeak.activeTrials_non = nan(d_info.numAnimals,1);
corr_speed_timeAtPeak.activeTrials_pos = nan(d_info.numAnimals,1);
corr_speed_timeAtPeak.activeTrials_neg = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
	if ~isempty(idcs.iscells{i})
        temp = [d{i,1}.pca_all.corr_speed_timeAtPeak_activeTrials_A(idcs.passed.Aonly{i},:);d{i,1}.pca_all.corr_speed_timeAtPeak_activeTrials_X(idcs.passed.Xonly{i},:)];
        corr_speed_timeAtPeak.activeTrials_non(i) = nansum(temp(:,2)>=0.05)/size(temp,1);
        corr_speed_timeAtPeak.activeTrials_pos(i) = nansum((temp(:,2)<0.05)&(temp(:,1)>0))/size(temp,1);
        corr_speed_timeAtPeak.activeTrials_neg(i) = nansum((temp(:,2)<0.05)&(temp(:,1)<0))/size(temp,1);
    end
end


%% Fig2_SpatialAnalysis_corr_pos

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_nonrunners = corr_speed_timeAtPeak.activeTrials_pos(d_info.running(:,1)==0);
this_data_runners = corr_speed_timeAtPeak.activeTrials_pos(d_info.running(:,1)==1);
this_speed_nonrunners = d_info.speed_sess(d_info.running(:,1)==0,1);
this_speed_runners = d_info.speed_sess(d_info.running(:,1)==1,1);

xline(10,'k:');
scatter(this_speed_nonrunners,this_data_nonrunners,'.','SizeData',100,'MarkerEdgeColor',p.col.nonrunner);
scatter(this_speed_runners,this_data_runners,'.','SizeData',100,'MarkerEdgeColor',p.col.runner);
%[this_corr_comb_r,this_corr_comb_p] = fitLine([this_speed_runners;this_speed_nonrunners],[this_data_runners;this_data_nonrunners],p.col.black)
[this_corr_r,this_corr_p] = fitLine(this_speed_runners,this_data_runners,p.col.black)
[pval,~,stats]=ranksum(this_data_nonrunners,this_data_runners)

xlim([-10,100])
xticks([0,50,100])
xlabel({'Running speed (cm/s)'})
% ylim([0,50])
% yticks([0,25,50])
% ytickformat('percentage')
ylabel(['Proportion of cells with positive correlation']);

savefig(F,[save_root_fig,'\Fig2_SpatialAnalysis_corr_pos.fig']);
saveas(F,[save_root_png,'\Fig2_SpatialAnalysis_corr_pos.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_SpatialAnalysis_corr_pos.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_SpatialAnalysis_corr_pos.txt'],'wt');


%% Fig2_SpatialAnalysis_corr_neg

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_nonrunners = corr_speed_timeAtPeak.activeTrials_neg(d_info.running(:,1)==0);
this_data_runners = corr_speed_timeAtPeak.activeTrials_neg(d_info.running(:,1)==1);
this_speed_nonrunners = d_info.speed_sess(d_info.running(:,1)==0,1);
this_speed_runners = d_info.speed_sess(d_info.running(:,1)==1,1);

xline(10,'k:');
scatter(this_speed_nonrunners,this_data_nonrunners,'.','SizeData',100,'MarkerEdgeColor',p.col.nonrunner);
scatter(this_speed_runners,this_data_runners,'.','SizeData',100,'MarkerEdgeColor',p.col.runner);
%[this_corr_comb_r,this_corr_comb_p] = fitLine([this_speed_runners;this_speed_nonrunners],[this_data_runners;this_data_nonrunners],p.col.black)
[this_corr_r,this_corr_p] = fitLine(this_speed_runners,this_data_runners,p.col.black)
[pval,~,stats]=ranksum(this_data_nonrunners,this_data_runners)

xlim([-10,100])
xticks([0,50,100])
xlabel({'Running speed (cm/s)'})
% ylim([0,50])
% yticks([0,25,50])
% ytickformat('percentage')
ylabel(['Proportion of cells with negative correlation']);

savefig(F,[save_root_fig,'\Fig2_SpatialAnalysis_corr_neg.fig']);
saveas(F,[save_root_png,'\Fig2_SpatialAnalysis_corr_neg.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_SpatialAnalysis_corr_neg.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_SpatialAnalysis_corr_neg.txt'],'wt');


%% Fig2_SpatialAnalysis_corr_non

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_nonrunners = corr_speed_timeAtPeak.activeTrials_non(d_info.running(:,1)==0);
this_data_runners = corr_speed_timeAtPeak.activeTrials_non(d_info.running(:,1)==1);
this_speed_nonrunners = d_info.speed_sess(d_info.running(:,1)==0,1);
this_speed_runners = d_info.speed_sess(d_info.running(:,1)==1,1);

xline(10,'k:');
scatter(this_speed_nonrunners,this_data_nonrunners,'.','SizeData',100,'MarkerEdgeColor',p.col.nonrunner);
scatter(this_speed_runners,this_data_runners,'.','SizeData',100,'MarkerEdgeColor',p.col.runner);
%[this_corr_comb_r,this_corr_comb_p] = fitLine([this_speed_runners;this_speed_nonrunners],[this_data_runners;this_data_nonrunners],p.col.black)
[this_corr_r,this_corr_p] = fitLine(this_speed_runners,this_data_runners,p.col.black)
[pval,~,stats]=ranksum(this_data_nonrunners,this_data_runners)

xlim([-10,100])
xticks([0,50,100])
xlabel({'Running speed (cm/s)'})
% ylim([0,50])
% yticks([0,25,50])
% ytickformat('percentage')
ylabel(['Proportion of cells with no correlation']);

savefig(F,[save_root_fig,'\Fig2_SpatialAnalysis_corr_non.fig']);
saveas(F,[save_root_png,'\Fig2_SpatialAnalysis_corr_non.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_SpatialAnalysis_corr_non.pdf']); set(gcf,'Color',[1,1,1])
fid = fopen([save_root_txt,'\Fig2_SpatialAnalysis_corr_non.txt'],'wt');



%% ----- performance corr -----

this_data = corr_speed_timeAtPeak.activeTrials_neg;
this_perf = correct;

figure; hold on
scatter(this_data,this_perf)
[this_corr_r,this_corr_p] = fitLine(this_data,this_perf,p.col.black)

xlabel('Prop cells w neg corr (all)')
ylabel('Performance')


this_data = corr_speed_timeAtPeak.activeTrials_neg(d_info.running(:,1)==1);
this_perf = correct(d_info.running(:,1)==1);

figure; hold on
scatter(this_data,this_perf)
[this_corr_r,this_corr_p] = fitLine(this_data,this_perf,p.col.black)

xlabel('Prop cells w neg corr (R)')
ylabel('Performance')
