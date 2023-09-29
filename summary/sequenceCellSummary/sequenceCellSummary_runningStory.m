%% Sequence Cell Summary - Running Story

%% Extract data

% % % fractionOfCells.A       = extractVariable(d,'tng_all.passed_stats.AW.A_fractionOfCells','array','single');
% % % fractionOfCells.X       = extractVariable(d,'tng_all.passed_stats.AW.X_fractionOfCells','array','single');
% % % fractionOfCells.Aonly   = extractVariable(d,'tng_all.passed_stats.AW.Aonly_fractionOfCells','array','single');
% % % fractionOfCells.Xonly   = extractVariable(d,'tng_all.passed_stats.AW.Xonly_fractionOfCells','array','single');
% % % fractionOfCells.AandX   = extractVariable(d,'tng_all.passed_stats.AW.AandX_fractionOfCells','array','single');
% % % fractionOfCells.AorX    = extractVariable(d,'tng_all.passed_stats.AW.AorX_fractionOfCells','array','single');
% % % fractionOfCells.notAorX = 1-extractVariable(d,'tng_all.passed_stats.AW.AorX_fractionOfCells','array','single');
% % % fractionOfCells = extractDay1(fractionOfCells);
% % % 
% % % % amplitude.A       = extractVariable(d,'tng_all.firingField.A_AW.meanAmplitude_blSub_ipsi','array','mean');
% % % % amplitude.X       = extractVariable(d,'tng_all.firingField.X_AW.meanAmplitude_blSub_ipsi','array','single');
% % % % amplitude = extractDay1(amplitude);
% % % 
% % % % running
% % % running = nan(d_info.numAnimals,1);
% % % for i=1:d_info.numAnimals
% % %     if isfield(d{i,1},'paq_beh') && (~isnan(nanmean(d{i,1}.paq_beh.speed)))
% % %         running(i) = nanmean(d{i,1}.paq_beh.speed) >= 10;
% % %     end
% % % end
% % % running(16) = true; % Arwen
% % % running(27) = false; % Stanage
% % % running(36:42) = 0 % Python, correct after data is imported

% speed
speed = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if isfield(d{i,1},'paq_beh') && (~isnan(nanmean(d{i,1}.paq_beh.speed)))
        speed(i) = nanmean(d{i,1}.paq_beh.speed);
    end
end


%% ---


nrows = 2;
ncols = 4;

F = default_figure([20,0.5,20,9.9]);



subplot(nrows,ncols,1)
these_data = [fractionOfCells.Aonly,fractionOfCells.Xonly,fractionOfCells.AandX,fractionOfCells.notAorX]*100;

hold on
h=bar(1:4,diag(nanmean(these_data,1)),'stacked');
h(1).FaceColor=p.col.A; h(2).FaceColor=p.col.X; h(3).FaceColor=p.col.darkGray; h(4).FaceColor=p.col.gray;
for i=1:size(these_data,1)
    if running(i)==1
        plot([1:4],these_data(i,:),'o-','Color',p.col.runner)
    elseif running(i)==0
        plot([1:4],these_data(i,:),'o-','Color',p.col.nonrunner)
    else
        plot([1:4],these_data(i,:),'o-','Color',p.col.gray)
    end
end
yline(50,'Color',p.col.darkGray,'LineStyle',':');
hold off
xlim([0,5])
xticks([1,2,3,4])
xticklabels({'A only','X only','A and X','none'})
ylim([0,100])
yticks([0,25,50,75,100])
ytickformat('percentage')
ylabel('Fraction of all cells')
title('Participation in sequences')



subplot(nrows,ncols,5)
these_data = [fractionOfCells.A,fractionOfCells.X,fractionOfCells.AorX,fractionOfCells.notAorX]*100;

hold on
h=bar(1:4,diag(nanmean(these_data,1)),'stacked');
h(1).FaceColor=p.col.A; h(2).FaceColor=p.col.X; h(3).FaceColor=p.col.darkGray; h(4).FaceColor=p.col.gray;
for i=1:size(these_data,1)
    if running(i)==1
        plot([1:4],these_data(i,:),'o-','Color',p.col.runner)
    elseif running(i)==0
        plot([1:4],these_data(i,:),'o-','Color',p.col.nonrunner)
    else
        plot([1:4],these_data(i,:),'o-','Color',p.col.gray)
    end
end
yline(50,'Color',p.col.darkGray,'LineStyle',':');
hold off
xlim([0,5])
xticks([1,2,3,4])
xticklabels({'A','X','A or X','none'})
ylim([0,100])
yticks([0,25,50,75,100])
ytickformat('percentage')
ylabel('Fraction of all cells')
title('Participation in sequences')


%%

subplot(nrows,ncols,2)
these_data = [fractionOfCells.Aonly_runner,fractionOfCells.Aonly_nonrunner,fractionOfCells.Xonly_runner,fractionOfCells.Xonly_nonrunner]*100;
these_cols = {p.col.A,p.col.A,p.col.X,p.col.X};

hold on
h=bar(1:4,diag(nanmean(these_data,1)),'stacked');
h(1).FaceColor=these_cols{1}; h(2).FaceColor=these_cols{2}; h(3).FaceColor=these_cols{3}; h(4).FaceColor=these_cols{4};
for i=1:size(these_data,1)
    if running(i)==1
        plot([1:4],these_data(i,:),'o-','Color',p.col.runner)
    elseif running(i)==0
        plot([1:4],these_data(i,:),'o-','Color',p.col.nonrunner)
    else
        plot([1:4],these_data(i,:),'o-','Color',p.col.gray)
    end
end
yline(50,'Color',p.col.darkGray,'LineStyle',':');
hold off
xlim([0,5])
xticks([1,2,3,4])
xticklabels({'A only','X only','A and X','none'})
ylim([0,100])
yticks([0,25,50,75,100])
ytickformat('percentage')
ylabel('Fraction of all cells')
title('Participation in sequences')
