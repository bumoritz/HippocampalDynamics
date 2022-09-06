%% Sequence Cell Summary (corrected using GLM) - Running Story

%% Extract running data

% running
running = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if isfield(d{i,1},'paq_beh') && (~isnan(nanmean(d{i,1}.paq_beh.speed)))
        running(i) = nanmean(d{i,1}.paq_beh.speed) >= 10;
    end
end
running(16) = true; % Arwen
running(27) = false; % Stanage
running(36:42) = 0 % Python, correct after data is imported

% speed
speed = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if isfield(d{i,1},'paq_beh') && (~isnan(nanmean(d{i,1}.paq_beh.speed)))
        speed(i) = nanmean(d{i,1}.paq_beh.speed);
    end
end


%% Extract tng data and clean with nem

temp = extractVariable(d,'tng_all.prop.iscell','cell','all');
iscells = nan(d_info.numAnimals,1);
for i=1:length(temp)
    if ~isempty(temp{i})
        iscells(i) = nansum(temp{i});
    end
end

passed.A = extractVariable(d,'tng_all.passed.AW.A','cell','all');
passed.X = extractVariable(d,'tng_all.passed.AW.X','cell','all');
passed.Aonly = extractVariable(d,'tng_all.passed.AW.Aonly','cell','all');
passed.Xonly = extractVariable(d,'tng_all.passed.AW.Xonly','cell','all');
passed.AandX = extractVariable(d,'tng_all.passed.AW.AandX','cell','all');
passed.AorX = extractVariable(d,'tng_all.passed.AW.AorX','cell','all');
passed.notAorX = passed.AorX;
for i=1:length(passed.notAorX)
    if ~isempty(passed.notAorX{i})
        passed.notAorX{i} = 1-passed.AorX{i};
    end
end
passed = extractDay1(passed);

temp = extractVariable(d,'nem_all_cmpr.sigM_sigTG_neg','cell','all');
negencoding = {};
negencoding.first = cell(length(temp),1);
negencoding.A = cell(length(temp),1);
negencoding.X = cell(length(temp),1);
negencoding.sens = cell(length(temp),1);
negencoding.Asens = cell(length(temp),1);
negencoding.Xsens = cell(length(temp),1);
negencoding.AXsens = cell(length(temp),1);
for i=1:length(temp)
    if ~isempty(temp{i})
        negencoding.first{i} = temp{i}{1};
        negencoding.A{i} = temp{i}{2};
        negencoding.X{i} = temp{i}{3};
        negencoding.sens{i} = temp{i}{4};
        negencoding.Asens{i} = temp{i}{6};
        negencoding.Xsens{i} = temp{i}{8};
        negencoding.AXsens{i} = ceil((temp{i}{6}+temp{i}{8})/2);
    end
end
negencoding = extractDay1(negencoding);

passed_cleaned = {};
passed_fields = fields(passed);
negencoding_fields = fields(negencoding);
for m=1:length(passed_fields)
    for n=1:length(negencoding_fields)
        passed_cleaned.([passed_fields{m},'_nneg',negencoding_fields{n}]) = cell(d_info.numAnimals,1);
        for i=1:d_info.numAnimals
            if (~isempty(passed.(passed_fields{m}){i})) && (~isempty(negencoding.(negencoding_fields{n}){i}))
                passed_cleaned.([passed_fields{m},'_nneg',negencoding_fields{n}]){i} = ...
                    setdiff(find(passed.(passed_fields{m}){i}==1),find(negencoding.(negencoding_fields{n}){i}==1));
            end
        end
    end
end

fractionOfCells = {};
passed_cleaned_fields = fields(passed_cleaned);
for m=1:length(passed_cleaned_fields)
    fractionOfCells.(passed_cleaned_fields{m}) = nan(d_info.numAnimals,1);
    for i=1:d_info.numAnimals
        if ~isempty(passed_cleaned.(passed_cleaned_fields{m}){i})
            fractionOfCells.(passed_cleaned_fields{m})(i) = length(passed_cleaned.(passed_cleaned_fields{m}){i}) / iscells(i);
        end
    end
end


%% Extract tng data and clean with threshold from tng

temp = extractVariable(d,'tng_all.prop.iscell','cell','all');
iscells = nan(d_info.numAnimals,1);
for i=1:length(temp)
    if ~isempty(temp{i})
        iscells(i) = nansum(temp{i});
    end
end

passed.A = extractVariable(d,'tng_all.passed.AW.A','cell','all');
passed.X = extractVariable(d,'tng_all.passed.AW.X','cell','all');
passed.Aonly = extractVariable(d,'tng_all.passed.AW.Aonly','cell','all');
passed.Xonly = extractVariable(d,'tng_all.passed.AW.Xonly','cell','all');
passed.AandX = extractVariable(d,'tng_all.passed.AW.AandX','cell','all');
passed.AorX = extractVariable(d,'tng_all.passed.AW.AorX','cell','all');
passed.notAorX = passed.AorX;
for i=1:length(passed.notAorX)
    if ~isempty(passed.notAorX{i})
        passed.notAorX{i} = 1-passed.AorX{i};
    end
end
passed = extractDay1(passed);

negthreshold.AX = cell(length(temp),1);
negthreshold.A = extractVariable(d,'tng_all.firingField.A_odour1.meanActivityInWindow_blSub','cell','all');
for i=1:length(negthreshold.A)
    if ~isempty(negthreshold.A{i})
        negthreshold.A{i} = negthreshold.A{i}<0;
    end
end
negthreshold.X = extractVariable(d,'tng_all.firingField.X_odour1.meanActivityInWindow_blSub','cell','all');
for i=1:length(negthreshold.X)
    if ~isempty(negthreshold.X{i})
        negthreshold.X{i} = negthreshold.X{i}<0;
        negthreshold.AX{i} = negthreshold.A{i} | negthreshold.X{i};
    end
end
negthreshold = extractDay1(negthreshold);

passed_cleaned = {};
passed_fields = fields(passed);
negthreshold_fields = fields(negthreshold);
for m=1:length(passed_fields)
    for n=1:length(negthreshold_fields)
        passed_cleaned.([passed_fields{m},'_nneg',negthreshold_fields{n}]) = cell(d_info.numAnimals,1);
        for i=1:d_info.numAnimals
            if (~isempty(passed.(passed_fields{m}){i})) && (~isempty(negthreshold.(negthreshold_fields{n}){i}))
                passed_cleaned.([passed_fields{m},'_nneg',negthreshold_fields{n}]){i} = ...
                    setdiff(find(passed.(passed_fields{m}){i}==1),find(negthreshold.(negthreshold_fields{n}){i}));
            end
        end
    end
end

fractionOfCells = {};
passed_cleaned_fields = fields(passed_cleaned);
for m=1:length(passed_cleaned_fields)
    fractionOfCells.(passed_cleaned_fields{m}) = nan(d_info.numAnimals,1);
    for i=1:d_info.numAnimals
        if ~isempty(passed_cleaned.(passed_cleaned_fields{m}){i})
            fractionOfCells.(passed_cleaned_fields{m})(i) = length(passed_cleaned.(passed_cleaned_fields{m}){i}) / iscells(i);
        end
    end
end


%% Animal-wise pie charts: Participation of cells in sequences - sorted by running speed

refdata = fractionOfCells.A_nnegA;
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
        this_data = [fractionOfCells.Aonly_nnegfirst(this_idx),fractionOfCells.Xonly_nnegfirst(this_idx),fractionOfCells.AandX_nnegfirst(this_idx),...
            1-sum([fractionOfCells.Aonly_nnegfirst(this_idx),fractionOfCells.Xonly_nnegfirst(this_idx),fractionOfCells.AandX_nnegfirst(this_idx)])];
        if nansum(this_data)<0.99 || nansum(this_data)>1.01
            warning('Fractions of cells dont add up to 1.')
        end
        h=pie(this_data);
        delete(findobj(h,'Type','text')); temp = findobj(h,'Type','Patch'); temp(1).FaceColor = p.col.A; temp(2).FaceColor = p.col.X; temp(3).FaceColor = p.col.darkGray; temp(4).FaceColor = p.col.gray;
        title([d_info.animals{this_idx},'-',d_info.dates_day1(this_idx,:),newline,num2str(speed(this_idx),2),' cm/s'])
    end
end

suptitle('Participation of cells in sequences - nnegfirst')






%% Animal-wise pie charts: Participation of cells in sequences - sorted by running speed

refdata = fractionOfCells.A_nnegA;
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
        this_data = [fractionOfCells.Aonly_nnegAsens(this_idx),fractionOfCells.Xonly_nnegXsens(this_idx),fractionOfCells.AandX_nnegAXsens(this_idx),...
            1-sum([fractionOfCells.Aonly_nnegAsens(this_idx),fractionOfCells.Xonly_nnegXsens(this_idx),fractionOfCells.AandX_nnegAXsens(this_idx)])];
        if nansum(this_data)<0.99 || nansum(this_data)>1.01
            warning('Fractions of cells dont add up to 1.')
        end
        h=pie(this_data);
        delete(findobj(h,'Type','text')); temp = findobj(h,'Type','Patch'); temp(1).FaceColor = p.col.A; temp(2).FaceColor = p.col.X; temp(3).FaceColor = p.col.darkGray; temp(4).FaceColor = p.col.gray;
        title([d_info.animals{this_idx},'-',d_info.dates_day1(this_idx,:),newline,num2str(speed(this_idx),2),' cm/s'])
    end
end

suptitle('Participation of cells in sequences - nneg ipsi sens')



%%

F = default_figure([20,0.5,20,9.9]);

these_data = [fractionOfCells.Aonly_nnegAsens,fractionOfCells.Xonly_nnegXsens,fractionOfCells.AandX_nnegAXsens,...
            1-nansum([fractionOfCells.Aonly_nnegAsens,fractionOfCells.Xonly_nnegXsens,fractionOfCells.AandX_nnegAXsens],2)]*100;

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



%%

%% Animal-wise pie charts: Participation of cells in sequences - sorted by running speed

refdata = fractionOfCells.A_nnegA;
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
        this_data = [fractionOfCells.Aonly_nnegA(this_idx),fractionOfCells.Xonly_nnegX(this_idx),fractionOfCells.AandX_nnegAX(this_idx),...
            1-sum([fractionOfCells.Aonly_nnegA(this_idx),fractionOfCells.Xonly_nnegX(this_idx),fractionOfCells.AandX_nnegAX(this_idx)])];
        if nansum(this_data)<0.99 || nansum(this_data)>1.01
            warning('Fractions of cells dont add up to 1.')
        end
        h=pie(this_data);
        delete(findobj(h,'Type','text')); temp = findobj(h,'Type','Patch'); temp(1).FaceColor = p.col.A; temp(2).FaceColor = p.col.X; temp(3).FaceColor = p.col.darkGray; temp(4).FaceColor = p.col.gray;
        title([d_info.animals{this_idx},'-',d_info.dates_day1(this_idx,:),newline,num2str(speed(this_idx),2),' cm/s'])
    end
end

suptitle('Participation of cells in sequences - nneg threshold')


