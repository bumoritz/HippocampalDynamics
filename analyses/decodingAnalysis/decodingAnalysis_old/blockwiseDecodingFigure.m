%% Average posterior snake plot

this_set = 'trainingSet';
this_trialSet = 'trialsForCompleteSet';
this_maxBin = 26;

n = 0; nrows = 2; ncols = 2; %length(dec_100t); 
F = default_figure([20,0.5,20,9.9]);


for i=1:length(dec_100t)
    for c=1:2
        
        % get data
        if c==1 % AB or AY
            these_trials = [find(dec_100t{i}.(this_set).trialType==1);find(dec_100t{i}.(this_set).trialType==3)];
        elseif c==2 % XY or XB
            these_trials = [find(dec_100t{i}.(this_set).trialType==2);find(dec_100t{i}.(this_set).trialType==4)];
        end
        this_data = nanmean(dec_100t{i}.(this_set).posterior(:,:,these_trials),3);
        
        % make subplot 
        n = n+1; subplot(nrows,ncols,n); hold on;
        imagesc(this_data,[0,0.1]);
        colormap('hot');
        %h = colorbar;
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

        temp = dec.(this_set).Y(min(find(dec.(this_set).trialType==c)),:);
        line([1-0.5,26+0.5],[temp(1)-0.5,temp(26)+0.5],'LineStyle','-','Color','g');
        line([27-0.5,52+0.5],[temp(27)-0.5,temp(52)+0.5],'LineStyle','-','Color','g');

        xlabel('Time (s)')
        ylabel('Decoded time (s) and trial type')
        title(['Average Bayesian decoding - ',these_conditions{c},' - ',num2str(length(these_trials)),' trials'])

    end
end








%%

this_set = 'trainingSet';
this_trialSet = 'trialsForTestSet_incorrM';

these_conditions = {'AB','XY','AY','XB'};
for c=1:4
    these_trials = intersect(dec.cv.(this_trialSet).combined,find(dec.(this_set).trialType==c));
    this_data = nanmean(dec.(this_set).posterior(:,:,these_trials),3);

    subplot(2,2,c)
    %imagesc(this_data,[0,prctile(this_decs.(['avgPosterior_type',num2str(c)'])(:),99.5)]);
    imagesc(this_data,[0,0.3]);
    colormap('hot');
    %h = colorbar;
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

    temp = dec.(this_set).Y(min(find(dec.(this_set).trialType==c)),:);
    line([1-0.5,26+0.5],[temp(1)-0.5,temp(26)+0.5],'LineStyle','-','Color','g');
    line([27-0.5,52+0.5],[temp(27)-0.5,temp(52)+0.5],'LineStyle','-','Color','g');

    xlabel('Time (s)')
    ylabel('Decoded time (s) and trial type')
    title(['Average Bayesian decoding - ',these_conditions{c},' - ',num2str(length(these_trials)),' trials'])
end
suptitle([this_set,', ',this_trialSet]);


