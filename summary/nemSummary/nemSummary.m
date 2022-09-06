%function nemSummary(d_info,d,ops,p,path)

if ~exist([path.root_summary,'plots\nem'],'dir')
    mkdir([path.root_summary,'plots\nem']);
end


%% Preparations

% this_d{1,1} = d{8,1}; this_d{1,2} = d{8,2}; this_d{1,3} = d{8,3}; this_d{1,4} = d{8,4}; this_d{1,5} = d{8,5};
% trck = load('D:\SniffinHippo\Repo\Arasaka\Arasaka_combined\Arasaka_combined_trck_15_rigid.mat'); trck = trck.trck;

nem_ref = d{2,1}.nem_all_cmpr;
numTestGroups = length(nem_ref.testGroup.label);


%% Neuronal encoding summary figure - activation or suppression

these_testGroups = 1:numTestGroups;

nrows = 2;
ncols = 1;
default_figure([20,0.5,20,9.9])

% a)
subplot(nrows,ncols,1)

this_data = nan(d_info.numAnimals,length(these_testGroups));
this_data2 = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    try
        this_numCells = sum(~isnan(d{i,1}.nem_all_cmpr.sigM));
        this_numSigModelCells = nansum(d{i,1}.nem_all_cmpr.sigM);
        for k=1:length(these_testGroups)
            this_numEncodingCells = nansum(d{i,1}.nem_all_cmpr.sigM_sigTG{these_testGroups(k)});
            this_data(i,k) = this_numEncodingCells / this_numCells;
        end
        this_data2(i) = this_numSigModelCells / this_numCells;
    catch
    end
end

h=bar(diag(nanmean([this_data2,this_data],1)*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for k=1:length(these_testGroups)
    set(h(k+1),'FaceColor',mean([nem_ref.testGroup.cols{k};p.col.white]))
end
set(h(1),'FaceColor',mean([p.col.black;p.col.white]))
hold on
for i=1:d_info.numAnimals
    plot([this_data2(i),this_data(i,:)]*100,'k-');
end
for k=1:length(these_testGroups)
    scatter((k+1)*ones(d_info.numAnimals,1),this_data(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
end
scatter(1*ones(d_info.numAnimals,1),this_data2*100,'MarkerEdgeColor','k','MarkerFaceColor','none');

ylim([0,100])
ytickformat('percentage')
xticks(1:length(these_testGroups)+1)
xticklabels([{'*'},nem_ref.testGroup.label(these_testGroups)])
ylabel('Proportion of all neurons')
set(gca,'box','off')
title('Significant encoding (activation or suppression)')

% b)
subplot(nrows,ncols,2)

this_data = nan(d_info.numAnimals,length(these_testGroups));
for i=1:d_info.numAnimals
    try
        this_numCells = sum(~isnan(d{i,1}.nem_all_cmpr.sigM));
        this_numSigModelCells = nansum(d{i,1}.nem_all_cmpr.sigM);
        for k=1:length(these_testGroups)
            this_numEncodingCells = nansum(d{i,1}.nem_all_cmpr.sigM_sigTG{these_testGroups(k)});
            this_data(i,k) = this_numEncodingCells / this_numSigModelCells;
        end
    catch
    end
end

h=bar(diag(nanmean(this_data,1)*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for k=1:length(these_testGroups)
    set(h(k),'FaceColor',mean([nem_ref.testGroup.cols{k};p.col.white]))
end
hold on
for i=1:d_info.numAnimals
    plot(this_data(i,:)*100,'k-');
end
for k=1:length(these_testGroups)
    scatter(k*ones(d_info.numAnimals,1),this_data(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
end

ylim([0,100])
ytickformat('percentage')
xticks(1:length(these_testGroups))
xticklabels(nem_ref.testGroup.label(these_testGroups))
ylabel('Proportion of neurons with significant fit')
set(gca,'box','off')
title('Significant encoding (activation or suppression)')


%% Neuronal encoding summary figure - activation

these_testGroups = 1:numTestGroups;

nrows = 2;
ncols = 1;
default_figure([20,0.5,20,9.9])

% a)
subplot(nrows,ncols,1)

this_data = nan(d_info.numAnimals,length(these_testGroups));
this_data2 = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    try
        this_numCells = sum(~isnan(d{i,1}.nem_all_cmpr.sigM));
        this_numSigModelCells = nansum(d{i,1}.nem_all_cmpr.sigM);
        for k=1:length(these_testGroups)
            this_numEncodingCells = nansum(d{i,1}.nem_all_cmpr.sigM_sigTG_pos{these_testGroups(k)});
            this_data(i,k) = this_numEncodingCells / this_numCells;
        end
        this_data2(i) = this_numSigModelCells / this_numCells;
    catch
    end
end

h=bar(diag(nanmean([this_data2,this_data],1)*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for k=1:length(these_testGroups)
    set(h(k+1),'FaceColor',mean([nem_ref.testGroup.cols{k};p.col.white]))
end
set(h(1),'FaceColor',mean([p.col.black;p.col.white]))
hold on
for i=1:d_info.numAnimals
    plot([this_data2(i),this_data(i,:)]*100,'k-');
end
for k=1:length(these_testGroups)
    scatter((k+1)*ones(d_info.numAnimals,1),this_data(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
end
scatter(1*ones(d_info.numAnimals,1),this_data2*100,'MarkerEdgeColor','k','MarkerFaceColor','none');

ylim([0,100])
ytickformat('percentage')
xticks(1:length(these_testGroups)+1)
xticklabels([{'*'},nem_ref.testGroup.label(these_testGroups)])
ylabel('Proportion of all neurons')
set(gca,'box','off')
title('Significant encoding (activation)')

% b)
subplot(nrows,ncols,2)

this_data = nan(d_info.numAnimals,length(these_testGroups));
for i=1:d_info.numAnimals
    try
        this_numCells = sum(~isnan(d{i,1}.nem_all_cmpr.sigM));
        this_numSigModelCells = nansum(d{i,1}.nem_all_cmpr.sigM);
        for k=1:length(these_testGroups)
            this_numEncodingCells = nansum(d{i,1}.nem_all_cmpr.sigM_sigTG_pos{these_testGroups(k)});
            this_data(i,k) = this_numEncodingCells / this_numSigModelCells;
        end
    catch
    end
end

h=bar(diag(nanmean(this_data,1)*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for k=1:length(these_testGroups)
    set(h(k),'FaceColor',mean([nem_ref.testGroup.cols{k};p.col.white]))
end
hold on
for i=1:d_info.numAnimals
    plot(this_data(i,:)*100,'k-');
end
for k=1:length(these_testGroups)
    scatter(k*ones(d_info.numAnimals,1),this_data(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
end

ylim([0,100])
ytickformat('percentage')
xticks(1:length(these_testGroups))
xticklabels(nem_ref.testGroup.label(these_testGroups))
ylabel('Proportion of neurons with significant fit')
set(gca,'box','off')
title('Significant encoding (activation)')


%% Neuronal encoding summary figure - suppression

these_testGroups = 1:numTestGroups;

nrows = 2;
ncols = 1;
default_figure([20,0.5,20,9.9])

% a)
subplot(nrows,ncols,1)

this_data = nan(d_info.numAnimals,length(these_testGroups));
this_data2 = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    try
        this_numCells = sum(~isnan(d{i,1}.nem_all_cmpr.sigM));
        this_numSigModelCells = nansum(d{i,1}.nem_all_cmpr.sigM);
        for k=1:length(these_testGroups)
            this_numEncodingCells = nansum(d{i,1}.nem_all_cmpr.sigM_sigTG_neg{these_testGroups(k)});
            this_data(i,k) = this_numEncodingCells / this_numCells;
        end
        this_data2(i) = this_numSigModelCells / this_numCells;
    catch
    end
end

h=bar(diag(nanmean([this_data2,this_data],1)*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for k=1:length(these_testGroups)
    set(h(k+1),'FaceColor',mean([nem_ref.testGroup.cols{k};p.col.white]))
end
set(h(1),'FaceColor',mean([p.col.black;p.col.white]))
hold on
for i=1:d_info.numAnimals
    plot([this_data2(i),this_data(i,:)]*100,'k-');
end
for k=1:length(these_testGroups)
    scatter((k+1)*ones(d_info.numAnimals,1),this_data(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
end
scatter(1*ones(d_info.numAnimals,1),this_data2*100,'MarkerEdgeColor','k','MarkerFaceColor','none');

ylim([0,100])
ytickformat('percentage')
xticks(1:length(these_testGroups)+1)
xticklabels([{'*'},nem_ref.testGroup.label(these_testGroups)])
ylabel('Proportion of all neurons')
set(gca,'box','off')
title('Significant encoding (suppression)')

% b)
subplot(nrows,ncols,2)

this_data = nan(d_info.numAnimals,length(these_testGroups));
for i=1:d_info.numAnimals
    try
        this_numCells = sum(~isnan(d{i,1}.nem_all_cmpr.sigM));
        this_numSigModelCells = nansum(d{i,1}.nem_all_cmpr.sigM);
        for k=1:length(these_testGroups)
            this_numEncodingCells = nansum(d{i,1}.nem_all_cmpr.sigM_sigTG_neg{these_testGroups(k)});
            this_data(i,k) = this_numEncodingCells / this_numSigModelCells;
        end
    catch
    end
end

h=bar(diag(nanmean(this_data,1)*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for k=1:length(these_testGroups)
    set(h(k),'FaceColor',mean([nem_ref.testGroup.cols{k};p.col.white]))
end
hold on
for i=1:d_info.numAnimals
    plot(this_data(i,:)*100,'k-');
end
for k=1:length(these_testGroups)
    scatter(k*ones(d_info.numAnimals,1),this_data(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
end

ylim([0,100])
ytickformat('percentage')
xticks(1:length(these_testGroups))
xticklabels(nem_ref.testGroup.label(these_testGroups))
ylabel('Proportion of neurons with significant fit')
set(gca,'box','off')
title('Significant encoding (suppression)')



