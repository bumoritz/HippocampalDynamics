%% Preparations

nem_ref = d{2,1}.nem_all_cmpr;

%%


nem_A = {};
nem_Asens = {};
nem_Atemp = {};
nem_X = {};
nem_Xsens = {};
nem_Xtemp = {};
tng_A = {};
tng_Asens = {};
tng_Atemp = {};
tng_X = {};
tng_Xsens = {};
tng_Xtemp = {};
numCells = {};
for i=1:d_info.numAnimals
    try
        
        % extract nem and tng
        nem_A{i,1} = d{i,1}.nem_all_cmpr.sigM_sigTG_pos{2};
        nem_Asens{i,1} = d{i,1}.nem_all_cmpr.sigM_sigTG_pos{6};
        nem_Atemp{i,1} = d{i,1}.nem_all_cmpr.sigM_sigTG_pos{7};
        nem_X{i,1} = d{i,1}.nem_all_cmpr.sigM_sigTG_pos{3};
        nem_Xsens{i,1} = d{i,1}.nem_all_cmpr.sigM_sigTG_pos{8};
        nem_Xtemp{i,1} = d{i,1}.nem_all_cmpr.sigM_sigTG_pos{9};
        tng_A{i,1} = d{i,1}.tng_all.passed.AW.A;
        tng_Asens{i,1} = d{i,1}.tng_all.passed.AW.A_early;
        tng_Atemp{i,1} = d{i,1}.tng_all.passed.AW.A_late;
        tng_X{i,1} = d{i,1}.tng_all.passed.AW.X;
        tng_Xsens{i,1} = d{i,1}.tng_all.passed.AW.X_early;
        tng_Xtemp{i,1} = d{i,1}.tng_all.passed.AW.X_late;
        numCells{i,1} = d{i,1}.tng_all.prop.numCells;
       
        % 
        
    catch
    end
end


%% Figure

nrows = 2;
ncols = 3;
default_figure([20,0.5,20,9.9])

% a) A
subplot(nrows,ncols,1)

this_data = nan(d_info.numAnimals,3);
for i=1:d_info.numAnimals
    try
        temp = intersect(find(nem_A{i}==1),find(tng_A{i}==1));
        this_data(i,1) = length(temp) / length(find(nem_A{i}==1));
        this_data(i,2) = length(temp) / length(find(tng_A{i}==1));
        this_data(i,3) = length(temp) / numCells{i};
    catch
    end
end

h=bar(diag(nanmean(this_data,1)*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for k=1:3
    set(h(k),'FaceColor',p.col.A) % mean([p.col.A;p.col.white])
end
hold on
for i=1:d_info.numAnimals
    plot(this_data(i,:)*100,'k-');
end
for k=1:3
    scatter(k*ones(d_info.numAnimals,1),this_data(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
end

ylim([0,100])
ytickformat('percentage')
xticks(1:3)
xticklabels({'GLM','seq analysis','all cells'})
ylabel('Proportion of neurons that pass both GLM and seq analysis')
set(gca,'box','off')
title(nem_ref.testGroup.label{2})

% b) A_sens
subplot(nrows,ncols,2)

this_data = nan(d_info.numAnimals,3);
for i=1:d_info.numAnimals
    try
        temp = intersect(find(nem_Asens{i}==1),find(tng_Asens{i}==1));
        this_data(i,1) = length(temp) / length(find(nem_Asens{i}==1));
        this_data(i,2) = length(temp) / length(find(tng_Asens{i}==1));
        this_data(i,3) = length(temp) / numCells{i};
    catch
    end
end

h=bar(diag(nanmean(this_data,1)*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for k=1:3
    set(h(k),'FaceColor',p.col.A) % mean([p.col.A;p.col.white])
end
hold on
for i=1:d_info.numAnimals
    plot(this_data(i,:)*100,'k-');
end
for k=1:3
    scatter(k*ones(d_info.numAnimals,1),this_data(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
end

ylim([0,100])
ytickformat('percentage')
xticks(1:3)
xticklabels({'GLM','seq analysis','all cells'})
ylabel('Proportion of neurons that pass both GLM and seq analysis')
set(gca,'box','off')
title(nem_ref.testGroup.label{6})

% c) A_temp
subplot(nrows,ncols,3)

this_data = nan(d_info.numAnimals,3);
for i=1:d_info.numAnimals
    try
        temp = intersect(find(nem_Atemp{i}==1),find(tng_Atemp{i}==1));
        this_data(i,1) = length(temp) / length(find(nem_Atemp{i}==1));
        this_data(i,2) = length(temp) / length(find(tng_Atemp{i}==1));
        this_data(i,3) = length(temp) / numCells{i};
    catch
    end
end

h=bar(diag(nanmean(this_data,1)*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for k=1:3
    set(h(k),'FaceColor',p.col.A) % mean([p.col.A;p.col.white])
end
hold on
for i=1:d_info.numAnimals
    plot(this_data(i,:)*100,'k-');
end
for k=1:3
    scatter(k*ones(d_info.numAnimals,1),this_data(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
end

ylim([0,100])
ytickformat('percentage')
xticks(1:3)
xticklabels({'GLM','seq analysis','all cells'})
ylabel('Proportion of neurons that pass both GLM and seq analysis')
set(gca,'box','off')
title(nem_ref.testGroup.label{7})


% d) X
subplot(nrows,ncols,4)

this_data = nan(d_info.numAnimals,3);
for i=1:d_info.numAnimals
    try
        temp = intersect(find(nem_X{i}==1),find(tng_X{i}==1));
        this_data(i,1) = length(temp) / length(find(nem_X{i}==1));
        this_data(i,2) = length(temp) / length(find(tng_X{i}==1));
        this_data(i,3) = length(temp) / numCells{i};
    catch
    end
end

h=bar(diag(nanmean(this_data,1)*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for k=1:3
    set(h(k),'FaceColor',p.col.X) % mean([p.col.A;p.col.white])
end
hold on
for i=1:d_info.numAnimals
    plot(this_data(i,:)*100,'k-');
end
for k=1:3
    scatter(k*ones(d_info.numAnimals,1),this_data(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
end

ylim([0,100])
ytickformat('percentage')
xticks(1:3)
xticklabels({'GLM','seq analysis','all cells'})
ylabel('Proportion of neurons that pass both GLM and seq analysis')
set(gca,'box','off')
title(nem_ref.testGroup.label{3})

% e) X_sens
subplot(nrows,ncols,5)

this_data = nan(d_info.numAnimals,3);
for i=1:d_info.numAnimals
    try
        temp = intersect(find(nem_Xsens{i}==1),find(tng_Xsens{i}==1));
        this_data(i,1) = length(temp) / length(find(nem_Xsens{i}==1));
        this_data(i,2) = length(temp) / length(find(tng_Xsens{i}==1));
        this_data(i,3) = length(temp) / numCells{i};
    catch
    end
end

h=bar(diag(nanmean(this_data,1)*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for k=1:3
    set(h(k),'FaceColor',p.col.X) % mean([p.col.A;p.col.white])
end
hold on
for i=1:d_info.numAnimals
    plot(this_data(i,:)*100,'k-');
end
for k=1:3
    scatter(k*ones(d_info.numAnimals,1),this_data(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
end

ylim([0,100])
ytickformat('percentage')
xticks(1:3)
xticklabels({'GLM','seq analysis','all cells'})
ylabel('Proportion of neurons that pass both GLM and seq analysis')
set(gca,'box','off')
title(nem_ref.testGroup.label{8})

% f) X_temp
subplot(nrows,ncols,6)

this_data = nan(d_info.numAnimals,3);
for i=1:d_info.numAnimals
    try
        temp = intersect(find(nem_Xtemp{i}==1),find(tng_Xtemp{i}==1));
        this_data(i,1) = length(temp) / length(find(nem_Xtemp{i}==1));
        this_data(i,2) = length(temp) / length(find(tng_Xtemp{i}==1));
        this_data(i,3) = length(temp) / numCells{i};
    catch
    end
end

h=bar(diag(nanmean(this_data,1)*100),'stacked','BaseValue',0,'FaceColor','none','EdgeColor','none');
for k=1:3
    set(h(k),'FaceColor',p.col.X) % mean([p.col.A;p.col.white])
end
hold on
for i=1:d_info.numAnimals
    plot(this_data(i,:)*100,'k-');
end
for k=1:3
    scatter(k*ones(d_info.numAnimals,1),this_data(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
end

ylim([0,100])
ytickformat('percentage')
xticks(1:3)
xticklabels({'GLM','seq analysis','all cells'})
ylabel('Proportion of neurons that pass both GLM and seq analysis')
set(gca,'box','off')
title(nem_ref.testGroup.label{9})








