%function sequenceCellSummary(d_info,d,ops,p,path)

if ~exist([path.root_summary,'plots\sequenceCells'],'dir')
    mkdir([path.root_summary,'plots\sequenceCells']);
end


%% Extract data

% fractionOfCells.A       = extractVariable(d,'tng_all.passed_stats.AW.A_fractionOfCells','array','single');
% fractionOfCells.X       = extractVariable(d,'tng_all.passed_stats.AW.X_fractionOfCells','array','single');
% fractionOfCells.Aonly   = extractVariable(d,'tng_all.passed_stats.AW.Aonly_fractionOfCells','array','single');
% fractionOfCells.Xonly   = extractVariable(d,'tng_all.passed_stats.AW.Xonly_fractionOfCells','array','single');
% fractionOfCells.AandX   = extractVariable(d,'tng_all.passed_stats.AW.AandX_fractionOfCells','array','single');
% fractionOfCells.AorX    = extractVariable(d,'tng_all.passed_stats.AW.AorX_fractionOfCells','array','single');
% fractionOfCells.notAorX = 1-extractVariable(d,'tng_all.passed_stats.AW.AorX_fractionOfCells','array','single');
% fractionOfCells = extractDay1(fractionOfCells);

fractionOfCells.A       = extractVariable(d,'tng_all_cmpr.passed_stats.A_fractionOfCells','array','single');
fractionOfCells.X       = extractVariable(d,'tng_all_cmpr.passed_stats.X_fractionOfCells','array','single');
fractionOfCells.Aonly   = extractVariable(d,'tng_all_cmpr.passed_stats.Aonly_fractionOfCells','array','single');
fractionOfCells.Xonly   = extractVariable(d,'tng_all_cmpr.passed_stats.Xonly_fractionOfCells','array','single');
fractionOfCells.AandX   = extractVariable(d,'tng_all_cmpr.passed_stats.AandX_fractionOfCells','array','single');
fractionOfCells.AorX    = extractVariable(d,'tng_all_cmpr.passed_stats.AorX_fractionOfCells','array','single');
fractionOfCells.notAorX = 1-extractVariable(d,'tng_all_cmpr.passed_stats.AorX_fractionOfCells','array','single');
fractionOfCells = extractDay1(fractionOfCells);

% numberOfCells.A         = extractVariable(d,'tng_all.passed_stats.AW.A_num','array','single');
% numberOfCells.X         = extractVariable(d,'tng_all.passed_stats.AW.X_num','array','single');
% numberOfCells.Aonly     = extractVariable(d,'tng_all.passed_stats.AW.Aonly_num','array','single');
% numberOfCells.Xonly     = extractVariable(d,'tng_all.passed_stats.AW.Xonly_num','array','single');
% numberOfCells.AandX     = extractVariable(d,'tng_all.passed_stats.AW.AandX_num','array','single');
% numberOfCells.notAorX   = 1-extractVariable(d,'tng_all.passed_stats.AW.AorX_num','array','single');
% numberOfCells = extractDay1(numberOfCells);


%% Animal-wise pie charts: Participation of cells in sequences

refdata = fractionOfCells.A;
ncols = 8;
nrows = ceil(length(rmmissing(refdata))/ncols);

F = default_figure([20,0.5,20,9.9]);

n=0;
for i=1:d_info.numAnimals
    if ~isnan(refdata(i))
        n=n+1;
        subplot(nrows,ncols,n)
        this_data = [fractionOfCells.Aonly(i),fractionOfCells.Xonly(i),fractionOfCells.AandX(i),fractionOfCells.notAorX(i)];
        if nansum(this_data)<0.99 || nansum(this_data)>1.01
            warning('Fractions of cells dont add up to 1.')
        end
        h=pie([fractionOfCells.Aonly(i),fractionOfCells.Xonly(i),fractionOfCells.AandX(i),fractionOfCells.notAorX(i)]);
        delete(findobj(h,'Type','text')); temp = findobj(h,'Type','Patch'); temp(1).FaceColor = p.col.A; temp(2).FaceColor = p.col.X; temp(3).FaceColor = p.col.darkGray; temp(4).FaceColor = p.col.gray;
        title([d_info.animals{i},'-',d_info.dates_day1(i,:)])
    end
end

suptitle('Participation of cells in sequences')


%%

% animals split
% - stim vs no stim
% - runners vs no runners

% within-animal split
% - A vs X
% - runnning vs not running


%%

nrows = 2;
ncols = 4;

F = default_figure([20,0.5,20,9.9]);



subplot(nrows,ncols,1)
these_data = [fractionOfCells.Aonly,fractionOfCells.Xonly,fractionOfCells.AandX,fractionOfCells.notAorX]*100;

hold on
h=bar(1:4,diag(nanmean(these_data,1)),'stacked');
h(1).FaceColor=p.col.A; h(2).FaceColor=p.col.X; h(3).FaceColor=p.col.darkGray; h(4).FaceColor=p.col.gray;
for i=1:size(these_data,1)
    plot([1:4],these_data(i,:),'o-k')
end
for i=1:d_info.numAnimals
    if running(i)==1
        plot([1:4],these_data(i,:),'-','Color',p.col.runner);
    elseif running(i)==0
        plot([1:4],these_data(i,:),'-','Color',p.col.nonrunner);
    else
        plot([1:4],these_data(i,:),'k-');
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
    plot([1:4],these_data(i,:),'o-k')
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



%% --- running correlation ---

running = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if isfield(d{i,1},'paq_beh') && (~isnan(nanmean(d{i,1}.paq_beh.speed)))
        running(i) = nanmean(d{i,1}.paq_beh.speed) > 20;
    end
end
running(16) = true; % Arwen
running(27) = false; % Stanage
running(36:42) = 0 % Python, correct after data is imported

speed = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if isfield(d{i,1},'paq_beh') && (~isnan(nanmean(d{i,1}.paq_beh.speed)))
        speed(i) = nanmean(d{i,1}.paq_beh.speed);
    end
end


%% Animal-wise pie charts: Participation of cells in sequences - sorted by running speed

refdata = fractionOfCells.A;
ncols = 8;
nrows = ceil(length(rmmissing(refdata))/ncols);

F = default_figure([20,0.5,20,9.9]);

[temp,this_order] = sort(speed);
this_order = this_order(1:min(find(isnan(temp)))-1);
n=0;
for i=1:length(this_order)
    this_idx = this_order(i);
    if ~isnan(refdata(this_idx))
        n=n+1;
        subplot(nrows,ncols,n)
        this_data = [fractionOfCells.Aonly(this_idx),fractionOfCells.Xonly(this_idx),fractionOfCells.AandX(this_idx),fractionOfCells.notAorX(this_idx)];
        if nansum(this_data)<0.99 || nansum(this_data)>1.01
            warning('Fractions of cells dont add up to 1.')
        end
        h=pie([fractionOfCells.Aonly(this_idx),fractionOfCells.Xonly(this_idx),fractionOfCells.AandX(this_idx),fractionOfCells.notAorX(this_idx)]);
        delete(findobj(h,'Type','text')); temp = findobj(h,'Type','Patch'); temp(1).FaceColor = p.col.A; temp(2).FaceColor = p.col.X; temp(3).FaceColor = p.col.darkGray; temp(4).FaceColor = p.col.gray;
        title([d_info.animals{this_idx},'-',d_info.dates_day1(this_idx,:),newline,num2str(speed(this_idx),2),' cm/s'])
    end
end

suptitle('Participation of cells in sequences')

