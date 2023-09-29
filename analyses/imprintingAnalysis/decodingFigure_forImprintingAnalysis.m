%function decodingFigure
dec = dec_all; 

%% Average posterior snake plot - training set

this_set = 'completeSet';
this_trialSet = 'trialsForTrainingSet'; %'trialsForTestSet_incorrM';
this_upper = 0.2;

default_figure();
these_conditions = {'A','X'};
for c=1:length(these_conditions)
    these_trials = intersect(dec.cv.(this_trialSet).combined,find(dec.(this_set).trialType==c));
    this_data = nanmean(dec.(this_set).posterior(:,:,these_trials),3);

    subplot(1,2,c)
    %imagesc(this_data,[0,prctile(this_decs.(['avgPosterior_type',num2str(c)'])(:),99.5)]);
    imagesc(this_data,[0,this_upper]);
    colormap('hot');
    %h = colorbar;
    taskLines(p,info,'full','decoding',1,1,0)
    xlim([0.5,26.5])

    n=0;
    these_ticks = nan(1,dec.supClasses.numSupClasses);
    for i=1:dec.supClasses.numSupClasses
        temp = length(dec.supClasses.bins{i});
        these_ticks(i) = n+temp/2+0.5;
        n=n+temp;
        yline(n+0.5,'LineStyle','-','Color',p.col.gray);
    end
    yticks(these_ticks);
    yticklabels(dec.supClasses.labels_extd);

    temp = dec.(this_set).Y(min(find(dec.(this_set).trialType==c)),:);
    line([1-0.5,26+0.5],[temp(1)-0.5,temp(26)+0.5],'LineStyle','-','Color','g');
    %line([27-0.5,52+0.5],[temp(27)-0.5,temp(52)+0.5],'LineStyle','-','Color','g');

    xlabel('Time (s)')
    ylabel('Decoded time (s) and trial type')
    title(['Average Bayesian decoding - ',these_conditions{c},' - ',num2str(length(these_trials)),' trials'])
end
suptitle([this_set,', ',this_trialSet,', clim [0,',num2str(this_upper),']']);


%% Average posterior snake plot - test set

this_set = 'completeSet';
this_trialSet = 'trialsForTestSet'; %'trialsForTestSet_incorrM';
this_upper = 0.2;

default_figure();
these_conditions = {'A','X'};
for c=1:length(these_conditions)
    these_trials = intersect(dec.cv.(this_trialSet).combined,find(dec.(this_set).trialType==c));
    this_data = nanmean(dec.(this_set).posterior(:,:,these_trials),3);

    subplot(1,2,c)
    %imagesc(this_data,[0,prctile(this_decs.(['avgPosterior_type',num2str(c)'])(:),99.5)]);
    imagesc(this_data,[0,this_upper]);
    colormap('hot');
    %h = colorbar;
    taskLines(p,info,'full','decoding',1,1,0)
    xlim([0.5,26.5])

    n=0;
    these_ticks = nan(1,dec.supClasses.numSupClasses);
    for i=1:dec.supClasses.numSupClasses
        temp = length(dec.supClasses.bins{i});
        these_ticks(i) = n+temp/2+0.5;
        n=n+temp;
        yline(n+0.5,'LineStyle','-','Color',p.col.gray);
    end
    yticks(these_ticks);
    yticklabels(dec.supClasses.labels_extd);

    temp = dec.(this_set).Y(min(find(dec.(this_set).trialType==c)),:);
    line([1-0.5,26+0.5],[temp(1)-0.5,temp(26)+0.5],'LineStyle','-','Color','g');
    %line([27-0.5,52+0.5],[temp(27)-0.5,temp(52)+0.5],'LineStyle','-','Color','g');

    xlabel('Time (s)')
    ylabel('Decoded time (s) and trial type')
    title(['Average Bayesian decoding - ',these_conditions{c},' - ',num2str(length(these_trials)),' trials'])
end
suptitle([this_set,', ',this_trialSet,', clim [0,',num2str(this_upper),']']);





%% Average posterior snake plot - training set     - replotting

this_set = 'trainingSet';
this_upper = 0.2;

default_figure();
these_conditions = {'A','X'};
for c=1:length(these_conditions)
    these_trials =  find(dec.(this_set).trialType==c);
    this_data = nanmean(dec.(this_set).posterior(:,:,these_trials),3);

    subplot(1,2,c)
    %imagesc(this_data,[0,prctile(this_decs.(['avgPosterior_type',num2str(c)'])(:),99.5)]);
    imagesc(this_data,[0,this_upper]);
    colormap('hot');
    %h = colorbar;
    taskLines(p,info,'full','decoding',1,1,0)
    xlim([0.5,26.5])

    n=0;
    these_ticks = nan(1,dec.supClasses.numSupClasses);
    for i=1:dec.supClasses.numSupClasses
        temp = length(dec.supClasses.bins{i});
        these_ticks(i) = n+temp/2+0.5;
        n=n+temp;
        yline(n+0.5,'LineStyle','-','Color',p.col.gray);
    end
    yticks(these_ticks);
    yticklabels(dec.supClasses.labels_extd);

    temp = dec.(this_set).Y(min(find(dec.(this_set).trialType==c)),:);
    line([1-0.5,26+0.5],[temp(1)-0.5,temp(26)+0.5],'LineStyle','-','Color','g');
    %line([27-0.5,52+0.5],[temp(27)-0.5,temp(52)+0.5],'LineStyle','-','Color','g');

    xlabel('Time (s)')
    ylabel('Decoded time (s) and trial type')
    title(['Average Bayesian decoding - ',these_conditions{c}])
end
suptitle([this_set,', clim [0,',num2str(this_upper),']']);















%% Single-trial posterior snake plot

this_set = 'completeSet';
this_trialSet = 'trialsForTestSet_incorrM';

for i=1:1

    this_trial = dec.cv.(this_trialSet).combined(i);
    this_data = posterior_within(:,:,this_trial); 
    % this_data = dec.(this_set).posterior(:,:,this_trial);

    default_figure();
    imagesc(this_data,[0,1]);
    colormap('hot');
    h = colorbar;
    taskLines(p,info,'full','decoding',1,1,0)

    n=0;
    these_ticks = nan(1,dec.supClasses.numSupClasses);
    for i=1:dec.supClasses.numSupClasses
        temp = length(dec.supClasses.bins{i});
        these_ticks(i) = n+temp/2+0.5;
        n=n+temp;
        yline(n+0.5,'LineStyle','-','Color',p.col.gray);
    end
    yticks(these_ticks);
    yticklabels(dec.supClasses.labels_extd);

    % correct labels
    temp = dec.(this_set).Y(this_trial,:);
    line([1-0.5,26+0.5],[temp(1)-0.5,temp(26)+0.5],'LineStyle','-','Color','g')
    line([27-0.5,52+0.5],[temp(27)-0.5,temp(52)+0.5],'LineStyle','-','Color','g')

    xlabel('Time (s)')
    ylabel('Decoded time (s) and trial type')
    
    title(['Single-trial Bayesian decoding, trial=',num2str(this_trial),', type=',num2str(dec.(this_set).trialType(this_trial)),', ',this_trialSet])
end


%% Properties

nrows = 2; ncols = 4;
F = default_figure();

subplot(nrows,ncols,1)
hold on
shadedErrorBar(1:dec_all.bins.numBins,nanmean(dec_all.analysis.testSet_incorrM.timeDecodingError_s,1),nanstd(dec_all.analysis.testSet_incorrM.timeDecodingError_s,[],1),...
    'lineProps','r');
shadedErrorBar(1:dec_all.bins.numBins,nanmean(dec_all.analysis.testSet_corrM.timeDecodingError_s,1),nanstd(dec_all.analysis.testSet_corrM.timeDecodingError_s,[],1),...
    'lineProps','g');
xline(26.5)
xlabel('Bin')
ylabel('Decoding error (s)')
title('Time decoding error')

subplot(nrows,ncols,2)
hold on
shadedErrorBar(1:dec_all.bins.numBins,nanmean(dec_all.analysis.testSet_incorrM.timeDecodingError_withinCat_s,1),nanstd(dec_all.analysis.testSet_incorrM.timeDecodingError_withinCat_s,[],1),...
    'lineProps','r');
shadedErrorBar(1:dec_all.bins.numBins,nanmean(dec_all.analysis.testSet_corrM.timeDecodingError_withinCat_s,1),nanstd(dec_all.analysis.testSet_corrM.timeDecodingError_withinCat_s,[],1),...
    'lineProps','g');
xline(26.5)
xlabel('Bin')
ylabel('Decoding error (s)')
title('Time decoding error (within category)')

subplot(nrows,ncols,3)
hold on
shadedErrorBar(1:dec_all.bins.numBins,nanmean(dec_all.analysis.testSet_incorrM.typeDecodingCorrect,1),nanstd(dec_all.analysis.testSet_incorrM.typeDecodingCorrect,[],1),...
    'lineProps','r');
shadedErrorBar(1:dec_all.bins.numBins,nanmean(dec_all.analysis.testSet_corrM.typeDecodingCorrect,1),nanstd(dec_all.analysis.testSet_corrM.typeDecodingCorrect,[],1),...
    'lineProps','g');
xline(26.5)
xlabel('Bin')
ylabel('Accuracy (P_{correct})')
title('Type decoding accuracy')


subplot(nrows,ncols,4)
hold on
shadedErrorBar(1:dec_all.bins.numBins,nanmean(dec_all.analysis.testSet_incorrM.typeDecodingCorrect_withinWindow,1),nanstd(dec_all.analysis.testSet_incorrM.typeDecodingCorrect_withinWindow,[],1),...
    'lineProps','r');
shadedErrorBar(1:dec_all.bins.numBins,nanmean(dec_all.analysis.testSet_corrM.typeDecodingCorrect_withinWindow,1),nanstd(dec_all.analysis.testSet_corrM.typeDecodingCorrect_withinWindow,[],1),...
    'lineProps','g');
xline(26.5)
xlabel('Bin')
ylabel('Accuracy (P_{correct})')
title('Type decoding accuracy (within window)')
