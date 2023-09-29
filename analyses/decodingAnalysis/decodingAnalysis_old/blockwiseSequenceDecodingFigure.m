%% Average posterior snake plot

this_set = 'trainingSet';

m = 0; nrows = 2; ncols = length(dec_100t); 
F = default_figure([-20,0.5,20,9.9]); % default_figure([20,0.5,20,9.9]);

these_conditions = fields(rmfield(dec_100t{i}.cv.trialsForCompleteSet,'combined'));
for c=1:length(these_conditions)
    for i=1:length(dec_100t)
        
        % get data
        if strcmp(these_conditions{c},'A')
            these_trials = [find(dec_100t{i}.(this_set).trialType==1);find(dec_100t{i}.(this_set).trialType==3)];
        elseif strcmp(these_conditions{c},'X')
            these_trials = [find(dec_100t{i}.(this_set).trialType==2);find(dec_100t{i}.(this_set).trialType==4)];
        end
        this_data = nanmean(dec_100t{i}.(this_set).posterior(:,:,these_trials),3);

        % make subplot 
        m = m+1; subplot(nrows,ncols,m); hold on;
        imagesc(this_data,[0,0.1]);
        colormap('hot');
        %h = colorbar;
        taskLines(p,info,'full','decoding',1,1,0)
        xlim([0,length(dec_100t{i}.classes.t)/dec_100t{i}.supClasses.numSupClasses]+0.5)
        
        n=0;
        these_ticks = nan(1,dec_100t{i}.supClasses.numSupClasses);
        for j=1:dec_100t{i}.supClasses.numSupClasses
            temp = length(dec_100t{i}.supClasses.bins{j});
            these_ticks(j) = n+temp/2+0.5;
            n=n+temp;
            yline(n+0.5,'LineStyle','-','Color',p.col.gray,'LineWidth',3);
        end
        yticks(these_ticks);
        yticklabels(dec_100t{i}.supClasses.labels_extd);
        ylim([0,dec_100t{i}.classes.numClasses+0.5])

        xlabel('Time (s)')
        ylabel('Decoded time (s) and trial type')
        title(['Block ',num2str(i),', ',these_conditions{c},' trials'])

    end
end





