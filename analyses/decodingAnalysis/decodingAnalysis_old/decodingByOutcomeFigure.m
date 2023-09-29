% function decodingByOutcomeFigure

%% Figure layout

nrows = 2; ncols = 4;
F = default_figure([20,0.5,20,9.9]);

this_set = 'completeSet';
this_trialSet_corre = 'trialsForTestSet_incorrM';

default_figure([20,0.5,20,9.9]);
these_conditions = {'AB','XY','AY','XB'};
k=0;
for i=1:4
    for o=1:2
        k=k+1;
        if o==1
            these_trials = intersect(trials_all.outcome.correct,find(dec.(this_set).trialType==i));
            this_outcome = 'correct';
        elseif o==2
            these_trials = intersect(trials_all.outcome.incorrect,find(dec.(this_set).trialType==i));
            this_outcome = 'incorrect';
        end
        this_data = nanmean(dec.(this_set).posterior(:,:,these_trials),3);

        subplot(nrows,ncols,k)
        %imagesc(this_data,[0,prctile(this_decs.(['avgPosterior_type',num2str(c)'])(:),99.5)]);
        imagesc(this_data,[0,0.1]); % [0,0.5]
        colormap('hot');
        %h = colorbar;
        taskLines(p,info,'full','decoding',1,1,0);

        n=0;
        these_ticks = nan(1,dec.supClasses.numSupClasses);
        for j=1:dec.supClasses.numSupClasses
            temp = length(dec.supClasses.bins{j});
            these_ticks(j) = n+temp/2+0.5;
            n=n+temp;
            yline(n+0.5,'LineStyle','-','Color',p.col.gray);
        end
        yticks(these_ticks);
        yticklabels(dec.supClasses.labels_extd);

%         temp = dec.(this_set).Y(min(find(dec.(this_set).trialType==i)),:);
%         line([1-0.5,26+0.5],[temp(1)-0.5,temp(26)+0.5],'LineStyle','-','Color','g');
%         line([27-0.5,52+0.5],[temp(27)-0.5,temp(52)+0.5],'LineStyle','-','Color','g');
%         
        xlim([0,41])
        xlabel('Time (s)')
        if k==1 || k==5
            ylabel('Decoded time (s) and trial type')
        end
        title([these_conditions{i},', ',this_outcome,' (',num2str(length(these_trials)),' trials)'])
    end
end
if info.stimSession
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-',info.stimType,', Bayesian decoding']);
else
    suptitle([info.animal,'-',info.date,'-d',num2str(info.expDay),'-','nostim',', Bayesian decoding']);
end
drawnow;




%% Quick improvisation for single trial

these_conditions = {'AB','XY','AY','XB'};
k=0;
for i=1:1
    for o=1:1
        k=k+1;
        this_data = nanmean(dec.(this_set).posterior(:,:,4),3);

        subplot(nrows,ncols,k)
        %imagesc(this_data,[0,prctile(this_decs.(['avgPosterior_type',num2str(c)'])(:),99.5)]);
        imagesc(this_data,[0,0.5]);
        colormap('hot');
        %h = colorbar;
        taskLines(p,info,'full','decoding',1,1,0);

        n=0;
        these_ticks = nan(1,dec.supClasses.numSupClasses);
        for j=1:dec.supClasses.numSupClasses
            temp = length(dec.supClasses.bins{j});
            these_ticks(j) = n+temp/2+0.5;
            n=n+temp;
            yline(n+0.5,'LineStyle','-','Color',p.col.gray);
        end
        yticks(these_ticks);
        yticklabels(dec.supClasses.labels_extd);

%         temp = dec.(this_set).Y(min(find(dec.(this_set).trialType==i)),:);
%         line([1-0.5,26+0.5],[temp(1)-0.5,temp(26)+0.5],'LineStyle','-','Color','g');
%         line([27-0.5,52+0.5],[temp(27)-0.5,temp(52)+0.5],'LineStyle','-','Color','g');
%         
        xlim([0,41])
        xlabel('Time (s)')
        if k==1 || k==5
            ylabel('Decoded time (s) and trial type')
        end
        title([these_conditions{i},', ',this_outcome,' (',num2str(length(these_trials)),' trials)'])
    end
end
drawnow;
