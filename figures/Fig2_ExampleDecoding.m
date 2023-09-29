%% Preparations

p = get_p; info = get_info;
path.root_summary = 'C:\SniffinHippo\Summary\';
save_root_fig = [path.root_summary,'figures\Fig2_fig\'];
save_root_png = [path.root_summary,'figures\Fig2_png\'];
save_root_pdf = [path.root_summary,'figures\Fig2_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig2_txt\'];

dec_NR = load('F:\SniffinHippo\AnalysisX\BullyBoy\BullyBoy_20220214\BullyBoy_20220214_dec_all.mat');
dec_NR = dec_NR.dec_all;
dec_R = load('F:\SniffinHippo\AnalysisX\Ao\Ao_20220402\Ao_20220402_dec_all.mat');
dec_R = dec_R.dec_all;


%% Fig2_ExampleDecoding_Nonrunner

dec = dec_NR;
this_data_A = nanmean(dec.trainingSet_seq.posterior(:,:,find(dec.trainingSet_seq.trialType==1)),3);
this_data_X = nanmean(dec.trainingSet_seq.posterior(:,:,find(dec.trainingSet_seq.trialType==2)),3);
this_data_ipsi = nanmean(cat(3,this_data_A(1:25,:),this_data_X(26:50,:)),3);
this_data_contra = nanmean(cat(3,this_data_A(26:50,:),this_data_X(1:25,:)),3);

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.5*34)]); hold on;

%subplot(1,2,1); hold on;
imagesc(this_data_ipsi(1:end-1,1:end-1),[0,0.07]);
daspect([1,1,1])
colormap('hot');
xlabel('Time (s)')
ylabel('Decoded time (s)')

xlim([0,25])
xticks([0:5:25])
xticklabels({'0','','','','','5'})
ylim([0,25])
yticks([0:5:25])
yticklabels({'0','','','','','5'})

% save plot
savefig(F,[save_root_fig,'\Fig2_ExampleDecoding_Nonrunner.fig']);
saveas(F,[save_root_png,'\Fig2_ExampleDecoding_Nonrunner.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleDecoding_Nonrunner.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_ExampleDecoding_Quantification_Nonrunner

dec = dec_NR;

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(0.7*34)]); hold on;

this_data = dec.analysis.trainingSet_seq.typeDecodingCorrect;
shadedErrorBar(1:dec.bins.numBins,nanmean(this_data,1)*100,nansem(this_data,1)*100,'lineProps',p.col.black)
xlim([0,25]); xticks([0:5:25]); xticklabels({'0','','','','','5'})
xlabel('Time (s)')
ytickformat('percentage')
ylim([50,100]); yticks([50,100]); yticklabels([5,1]);
ylabel('1st Odor decoding accuracy')

% save plot
savefig(F,[save_root_fig,'\Fig2_ExampleDecoding_Quantification_Nonrunner.fig']);
saveas(F,[save_root_png,'\Fig2_ExampleDecoding_Quantification_Nonrunner.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleDecoding_Quantification_Nonrunner.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_ExampleDecoding_Runner

dec = dec_R;
this_data_A = nanmean(dec.trainingSet_seq.posterior(:,:,find(dec.trainingSet_seq.trialType==1)),3);
this_data_X = nanmean(dec.trainingSet_seq.posterior(:,:,find(dec.trainingSet_seq.trialType==2)),3);
this_data_ipsi = nanmean(cat(3,this_data_A(1:25,:),this_data_X(26:50,:)),3);
this_data_contra = nanmean(cat(3,this_data_A(26:50,:),this_data_X(1:25,:)),3);

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(1.5*34)]); hold on;

%subplot(1,2,1); hold on;
imagesc(this_data_ipsi(1:end-1,1:end-1),[0,0.07]);
daspect([1,1,1])
colormap('hot');
xlabel('Time (s)')
ylabel('Decoded time (s)')

xlim([0,25])
xticks([0:5:25])
xticklabels({'0','','','','','5'})
ylim([0,25])
yticks([0:5:25])
yticklabels({'0','','','','','5'})

% save plot
savefig(F,[save_root_fig,'\Fig2_ExampleDecoding_Runner.fig']);
saveas(F,[save_root_png,'\Fig2_ExampleDecoding_Runner.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleDecoding_Runner.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_ExampleDecoding_Quantification_Runner

dec = dec_R;

F = paper_figure([0,0.5,mm2inch(1.5*34),mm2inch(0.7*34)]); hold on;

this_data = dec.analysis.trainingSet_seq.typeDecodingCorrect;
shadedErrorBar(1:dec.bins.numBins,nanmean(this_data,1)*100,nansem(this_data,1)*100,'lineProps',p.col.black)
xlim([0,25]); xticks([0:5:25]); xticklabels({'0','','','','','5'})
xlabel('Time (s)')
ytickformat('percentage')
ylim([50,100]); yticks([50,100]); yticklabels([5,1]);
ylabel('1st Odor decoding accuracy')

% save plot
savefig(F,[save_root_fig,'\Fig2_ExampleDecoding_Quantification_Runner.fig']);
saveas(F,[save_root_png,'\Fig2_ExampleDecoding_Quantification_Runner.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleDecoding_Quantification_Runner.pdf']); set(gcf,'Color',[1,1,1])


%% Fig2_ExampleDecoding_Colorbar

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;
imagesc(this_data_ipsi,[0,0.07])
colormap('hot')
h=colorbar;

savefig(F,[save_root_fig,'\Fig2_ExampleDecoding_Colorbar.fig']);
saveas(F,[save_root_png,'\Fig2_ExampleDecoding_Colorbar.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig2_ExampleDecoding_Colorbar.pdf']); set(gcf,'Color',[1,1,1])
