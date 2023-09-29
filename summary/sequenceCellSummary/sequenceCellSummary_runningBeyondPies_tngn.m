%% Extract running data

% running
running = d_info.running(:,1);
% running = nan(d_info.numAnimals,1);
% for i=1:d_info.numAnimals
%     if isfield(d{i,1},'paq_beh') && (~isnan(nanmean(d{i,1}.paq_beh.speed)))
%         running(i) = nanmean(d{i,1}.paq_beh.speed) >= 10;
%     end
% end
% running(16) = true; % Arwen

% speed
speed = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if isfield(d{i,1},'paq_beh') && (~isnan(nanmean(d{i,1}.paq_beh.speed)))
        speed(i) = nanmean(d{i,1}.paq_beh.speed);
    end
end


%% Extract idcs

% idcs.iscells
temp = extractVariable(d,'tngn_all.prop.iscell','cell','all');
idcs.iscells = cell(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    try
        idcs.iscells{i} = find(temp{i}==1);
    catch
    end
end

% idcs.passed
these_conditions = {'A','X','Aonly','Xonly','AandX','AorX'};
for k=1:length(these_conditions)
    temp = extractVariable(d,['tngn_all.passed.AW.',these_conditions{k}],'cell','all');
    idcs.passed.(these_conditions{k}) = cell(d_info.numAnimals,1);
    for i=1:d_info.numAnimals
        try
            idcs.passed.(these_conditions{k}){i} = find(temp{i}==1);
        catch
        end
    end
end
idcs.passed.notAorX = idcs.iscells;
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        idcs.passed.notAorX{i} = setdiff(idcs.iscells{i},idcs.passed.AorX{i});
    end
end
idcs.passed = extractDay1(idcs.passed);

% idcs.negthreshold
% idcs.negthreshold.A = extractVariable(d,'tngn_all.firingField.A_odour1.meanActivityInWindow_blSub','cell','all');
% for i=1:d_info.numAnimals
%     if ~isempty(idcs.negthreshold.A{i})
%         idcs.negthreshold.A{i} = find((idcs.negthreshold.A{i}<0)==1);
%     end
% end
% idcs.negthreshold.X = extractVariable(d,'tngn_all.firingField.X_odour1.meanActivityInWindow_blSub','cell','all');
% for i=1:d_info.numAnimals
%     if ~isempty(idcs.negthreshold.X{i})
%         idcs.negthreshold.X{i} = find((idcs.negthreshold.X{i}<0)==1);
%     end
% end
% idcs.negthreshold.AX = cell(d_info.numAnimals,1);
% for i=1:d_info.numAnimals
%     if ~isempty(idcs.iscells{i})
%         idcs.negthreshold.AX{i} = unique([idcs.negthreshold.A{i};idcs.negthreshold.X{i}]);
%     end
% end

% idcs.passed_cleaned
idcs.passed_cleaned = idcs.passed;
% idcs.passed_cleaned = {};
% passed_fields = fields(idcs.passed);
% negthreshold_fields = fields(idcs.negthreshold);
% for m=1:length(passed_fields)
%     for n=1:length(negthreshold_fields)
%         idcs.passed_cleaned.([passed_fields{m},'_nneg',negthreshold_fields{n}]) = cell(d_info.numAnimals,1);
%         for i=1:d_info.numAnimals
%             if (~isempty(idcs.passed.(passed_fields{m}){i})) && (~isempty(idcs.negthreshold.(negthreshold_fields{n}){i}))
%                 idcs.passed_cleaned.([passed_fields{m},'_nneg',negthreshold_fields{n}]){i} = ...
%                     setdiff(idcs.passed.(passed_fields{m}){i},idcs.negthreshold.(negthreshold_fields{n}){i});
%             end
%         end
%     end
% end


%% Extract and bin peak time

temp_A = extractVariable(d,'tngn_all.firingField.A_AW.peakLocation_s','cell','all');
temp_X = extractVariable(d,'tngn_all.firingField.X_AW.peakLocation_s','cell','all');
these_conditions_cells = {'A','X','Aonly','Xonly','AandX','AorX'};
these_conditions_neg = {'A','X','AX'};
for m=1:length(these_conditions_cells)
    for n=1:length(these_conditions_neg)
        peakTime.(['Atrials_',these_conditions_cells{m}]) = cell(d_info.numAnimals,1);
        peakTime.(['Xtrials_',these_conditions_cells{m}]) = cell(d_info.numAnimals,1);
        for i=1:d_info.numAnimals
            if ~isempty(idcs.passed_cleaned.(these_conditions_cells{m}){i})
                peakTime.(['Atrials_',these_conditions_cells{m}]){i} = ...
                    temp_A{i}(idcs.passed_cleaned.(these_conditions_cells{m}){i});
                peakTime.(['Xtrials_',these_conditions_cells{m}]){i} = ...
                    temp_X{i}(idcs.passed_cleaned.(these_conditions_cells{m}){i});
            end
        end
    end
end

binEdges = [0:0.1:5]; %[0:4,5.5]%[0:0.5:5.5];
these_fields = fields(peakTime);
for k=1:length(these_fields)
    binnedPeakTime.(these_fields{k}) = nan(d_info.numAnimals,length(binEdges)-1);
    for i=1:d_info.numAnimals
        if ~isempty(peakTime.(these_fields{k}){i})
            temp = discretize(peakTime.(these_fields{k}){i},binEdges);
            for n=1:length(binEdges)-1
                binnedPeakTime.(these_fields{k})(i,n) = nansum(temp==n) / length(idcs.iscells{i});
            end
        end
    end
end


%% Figure

nrows = 2; ncols = 3;
F = default_figure();

% Atrials_A
this_data_runners = binnedPeakTime.Atrials_A(find(running==1),:)*100;
this_data_nonrunners = binnedPeakTime.Atrials_A(find(running==0),:)*100;

subplot(nrows,ncols,1)
hold on
shadedErrorBar((1:length(binEdges)-1)+0.5,nanmean(this_data_runners,1),nansem(this_data_runners,1),'lineProps',p.col.runner); % dark brown: runners
shadedErrorBar((1:length(binEdges)-1)+0.5,nanmean(this_data_nonrunners,1),nansem(this_data_nonrunners,1),'lineProps',p.col.nonrunner);  % light brown: non-runners
xlim([1,length(binEdges)])
xticks(1:length(binEdges));
xticklabels(strtrim(cellstr(num2str(binEdges'))'));
xlabel('Time (s)')
ytickformat('percentage')
ylabel('Proportion of neurons')
title('A sequence cells (w/o -) in A trials')

% Atrials_Aonly
this_data_runners = binnedPeakTime.Atrials_Aonly(find(running==1),:)*100;
this_data_nonrunners = binnedPeakTime.Atrials_Aonly(find(running==0),:)*100;

subplot(nrows,ncols,2)
hold on
shadedErrorBar((1:length(binEdges)-1)+0.5,nanmean(this_data_runners,1),nansem(this_data_runners,1),'lineProps',p.col.runner); % dark brown: runners
shadedErrorBar((1:length(binEdges)-1)+0.5,nanmean(this_data_nonrunners,1),nansem(this_data_nonrunners,1),'lineProps',p.col.nonrunner);  % light brown: non-runners
xlim([1,length(binEdges)])
xticks(1:length(binEdges));
xticklabels(strtrim(cellstr(num2str(binEdges'))'));
xlabel('Time (s)')
ytickformat('percentage')
ylabel('Proportion of neurons')
title('Aonly sequence cells (w/o -) in A trials')

% Atrials_AandX
this_data_runners = binnedPeakTime.Atrials_AandX(find(running==1),:)*100;
this_data_nonrunners = binnedPeakTime.Atrials_AandX(find(running==0),:)*100;

subplot(nrows,ncols,3)
hold on
shadedErrorBar((1:length(binEdges)-1)+0.5,nanmean(this_data_runners,1),nansem(this_data_runners,1),'lineProps',p.col.runner); % dark brown: runners
shadedErrorBar((1:length(binEdges)-1)+0.5,nanmean(this_data_nonrunners,1),nansem(this_data_nonrunners,1),'lineProps',p.col.nonrunner);  % light brown: non-runners
xlim([1,length(binEdges)])
xticks(1:length(binEdges));
xticklabels(strtrim(cellstr(num2str(binEdges'))'));
xlabel('Time (s)')
ytickformat('percentage')
ylabel('Proportion of neurons')
title('AandX sequence cells (w/o -) in A trials')

% Xtrials_X
this_data_runners = binnedPeakTime.Xtrials_X(find(running==1),:)*100;
this_data_nonrunners = binnedPeakTime.Xtrials_X(find(running==0),:)*100;

subplot(nrows,ncols,ncols+1)
hold on
shadedErrorBar((1:length(binEdges)-1)+0.5,nanmean(this_data_runners,1),nansem(this_data_runners,1),'lineProps',p.col.runner); % dark brown: runners
shadedErrorBar((1:length(binEdges)-1)+0.5,nanmean(this_data_nonrunners,1),nansem(this_data_nonrunners,1),'lineProps',p.col.nonrunner);  % light brown: non-runners
xlim([1,length(binEdges)])
xticks(1:length(binEdges));
xticklabels(strtrim(cellstr(num2str(binEdges'))'));
xlabel('Time (s)')
ytickformat('percentage')
ylabel('Proportion of neurons')
title('X sequence cells (w/o -) in X trials')

% Xtrials_Xonly
this_data_runners = binnedPeakTime.Xtrials_Xonly(find(running==1),:)*100;
this_data_nonrunners = binnedPeakTime.Xtrials_Xonly(find(running==0),:)*100;

subplot(nrows,ncols,ncols+2)
hold on
shadedErrorBar((1:length(binEdges)-1)+0.5,nanmean(this_data_runners,1),nansem(this_data_runners,1),'lineProps',p.col.runner); % dark brown: runners
shadedErrorBar((1:length(binEdges)-1)+0.5,nanmean(this_data_nonrunners,1),nansem(this_data_nonrunners,1),'lineProps',p.col.nonrunner);  % light brown: non-runners
xlim([1,length(binEdges)])
xticks(1:length(binEdges));
xticklabels(strtrim(cellstr(num2str(binEdges'))'));
xlabel('Time (s)')
ytickformat('percentage')
ylabel('Proportion of neurons')
title('Xonly sequence cells (w/o -) in X trials')

% Xtrials_AandX
this_data_runners = binnedPeakTime.Xtrials_AandX(find(running==1),:)*100;
this_data_nonrunners = binnedPeakTime.Xtrials_AandX(find(running==0),:)*100;

subplot(nrows,ncols,ncols+3)
hold on
shadedErrorBar((1:length(binEdges)-1)+0.5,nanmean(this_data_runners,1),nansem(this_data_runners,1),'lineProps',p.col.runner); % dark brown: runners
shadedErrorBar((1:length(binEdges)-1)+0.5,nanmean(this_data_nonrunners,1),nansem(this_data_nonrunners,1),'lineProps',p.col.nonrunner);  % light brown: non-runners
xlim([1,length(binEdges)])
xticks(1:length(binEdges));
xticklabels(strtrim(cellstr(num2str(binEdges'))'));
xlabel('Time (s)')
ytickformat('percentage')
ylabel('Proportion of neurons')
title('AandX sequence cells (w/o -) in X trials')


%% Most important plot

F = default_figure();

this_data_runners = nanmean(cat(3,binnedPeakTime.Atrials_Aonly(running==1,:)*100,binnedPeakTime.Xtrials_Xonly(find(running==1),:)*100),3);
this_data_nonrunners = nanmean(cat(3,binnedPeakTime.Atrials_Aonly(find(running==0),:)*100,binnedPeakTime.Xtrials_Xonly(find(running==0),:)*100),3);

hold on
shadedErrorBar((1:length(binEdges)-1)+0.5,nanmean(this_data_runners,1),nansem(this_data_runners,1),'lineProps',p.col.runner); % dark brown: runners
shadedErrorBar((1:length(binEdges)-1)+0.5,nanmean(this_data_nonrunners,1),nansem(this_data_nonrunners,1),'lineProps',p.col.nonrunner);  % light brown: non-runners
xlim([1,length(binEdges)])
xticks(1:length(binEdges));
xticklabels(strtrim(cellstr(num2str(binEdges'))'));
legend('','Runners','','Non-runners')
xlabel('Time (s)')
ytickformat('percentage')
ylabel('Proportion of neurons')
title('Sequence cells')



%% Figure - individual traces

nrows = 2; ncols = 3;
F = default_figure();

% Atrials_A
this_data_runners = binnedPeakTime.Atrials_A(find(running==1),:)*100;
this_data_nonrunners = binnedPeakTime.Atrials_A(find(running==0),:)*100;

subplot(nrows,ncols,1)
hold on
plot((1:length(binEdges)-1)+0.5,this_data_runners','Color',p.col.runner); % dark brown: runners
plot((1:length(binEdges)-1)+0.5,this_data_nonrunners','Color',p.col.nonrunner); % light brown: non-runners
xlim([1,length(binEdges)])
xticks(1:length(binEdges));
xticklabels(strtrim(cellstr(num2str(binEdges'))'));
xlabel('Time (s)')
ytickformat('percentage')
ylabel('Proportion of neurons')
title('A sequence cells (w/o -) in A trials')

% Atrials_Aonly
this_data_runners = binnedPeakTime.Atrials_Aonly(find(running==1),:)*100;
this_data_nonrunners = binnedPeakTime.Atrials_Aonly(find(running==0),:)*100;

subplot(nrows,ncols,2)
hold on
plot((1:length(binEdges)-1)+0.5,this_data_runners','Color',p.col.runner); % dark brown: runners
plot((1:length(binEdges)-1)+0.5,this_data_nonrunners','Color',p.col.nonrunner); % light brown: non-runners
xlim([1,length(binEdges)])
xticks(1:length(binEdges));
xticklabels(strtrim(cellstr(num2str(binEdges'))'));
xlabel('Time (s)')
ytickformat('percentage')
ylabel('Proportion of neurons')
title('Aonly sequence cells (w/o -) in A trials')

% Atrials_AandX
this_data_runners = binnedPeakTime.Atrials_AandX(find(running==1),:)*100;
this_data_nonrunners = binnedPeakTime.Atrials_AandX(find(running==0),:)*100;

subplot(nrows,ncols,3)
hold on
plot((1:length(binEdges)-1)+0.5,this_data_runners','Color',p.col.runner); % dark brown: runners
plot((1:length(binEdges)-1)+0.5,this_data_nonrunners','Color',p.col.nonrunner); % light brown: non-runners
xlim([1,length(binEdges)])
xticks(1:length(binEdges));
xticklabels(strtrim(cellstr(num2str(binEdges'))'));
xlabel('Time (s)')
ytickformat('percentage')
ylabel('Proportion of neurons')
title('AandX sequence cells (w/o -) in A trials')

% Xtrials_X
this_data_runners = binnedPeakTime.Xtrials_X(find(running==1),:)*100;
this_data_nonrunners = binnedPeakTime.Xtrials_X(find(running==0),:)*100;

subplot(nrows,ncols,ncols+1)
hold on
plot((1:length(binEdges)-1)+0.5,this_data_runners','Color',p.col.runner); % dark brown: runners
plot((1:length(binEdges)-1)+0.5,this_data_nonrunners','Color',p.col.nonrunner); % light brown: non-runners
xlim([1,length(binEdges)])
xticks(1:length(binEdges));
xticklabels(strtrim(cellstr(num2str(binEdges'))'));
xlabel('Time (s)')
ytickformat('percentage')
ylabel('Proportion of neurons')
title('X sequence cells (w/o -) in X trials')

% Xtrials_Xonly
this_data_runners = binnedPeakTime.Xtrials_Xonly(find(running==1),:)*100;
this_data_nonrunners = binnedPeakTime.Xtrials_Xonly(find(running==0),:)*100;

subplot(nrows,ncols,ncols+2)
hold on
plot((1:length(binEdges)-1)+0.5,this_data_runners','Color',p.col.runner); % dark brown: runners
plot((1:length(binEdges)-1)+0.5,this_data_nonrunners','Color',p.col.nonrunner); % light brown: non-runners
xlim([1,length(binEdges)])
xticks(1:length(binEdges));
xticklabels(strtrim(cellstr(num2str(binEdges'))'));
xlabel('Time (s)')
ytickformat('percentage')
ylabel('Proportion of neurons')
title('Xonly sequence cells (w/o -) in X trials')

% Xtrials_AandX
this_data_runners = binnedPeakTime.Xtrials_AandX(find(running==1),:)*100;
this_data_nonrunners = binnedPeakTime.Xtrials_AandX(find(running==0),:)*100;

subplot(nrows,ncols,ncols+3)
hold on
plot((1:length(binEdges)-1)+0.5,this_data_runners','Color',p.col.runner); % dark brown: runners
plot((1:length(binEdges)-1)+0.5,this_data_nonrunners','Color',p.col.nonrunner); % light brown: non-runners
xlim([1,length(binEdges)])
xticks(1:length(binEdges));
xticklabels(strtrim(cellstr(num2str(binEdges'))'));
xlabel('Time (s)')
ytickformat('percentage')
ylabel('Proportion of neurons')
title('AandX sequence cells (w/o -) in X trials')

