%% Data import

this_d = {};
trck = {};

% i=1; this_animal = 3; % Biontech
% this_d{i,1} = d{this_animal,1}; this_d{i,2} = d{this_animal,2}; this_d{i,3} = d{this_animal,3}; this_d{i,4} = d{this_animal,4}; this_d{i,5} = d{this_animal,5};
% temp = load('D:\SniffinHippo\Repo\Biontech\Biontech_combined\Biontech_combined_trck_15_rigid.mat'); trck{i,1} = temp.trck;

i=1; this_animal = 8; % Arasaka
this_d{i,1} = d{this_animal,1}; this_d{i,2} = d{this_animal,2}; this_d{i,3} = d{this_animal,3}; this_d{i,4} = d{this_animal,4}; this_d{i,5} = d{this_animal,5};
temp = load('D:\SniffinHippo\Repo\Arasaka\Arasaka_combined\Arasaka_combined_trck_15_rigid.mat'); trck{i,1} = temp.trck;

i=2; this_animal = 31; % Maria
this_d{i,1} = d{this_animal,1}; this_d{i,2} = d{this_animal,2}; this_d{i,3} = d{this_animal,3}; this_d{i,4} = d{this_animal,4}; this_d{i,5} = d{this_animal,5};
temp = load('D:\SniffinHippo\Repo\Maria\Maria_combined\Maria_combined_trck_15_nonrigid.mat'); trck{i,1} = temp.trck;

i=3; this_animal = 33; % William
this_d{i,1} = d{this_animal,1}; this_d{i,2} = d{this_animal,2}; this_d{i,3} = d{this_animal,3}; this_d{i,4} = d{this_animal,4}; this_d{i,5} = d{this_animal,5};
temp = load('D:\SniffinHippo\Repo\William\William_combined\William_combined_trck_15_nonrigid.mat'); trck{i,1} = temp.trck;


%% Preparations

nem_ref = this_d{1,1}.nem_all_cmpr;
numTestGroups = length(nem_ref.testGroup.label);

idcs = {};
numCells = nan(size(this_d,1),1);
for i=1:size(this_d,1)
    idcs{i,1} = find(sum(trck{i}.cell_to_index_map>0,2)==size(trck{i}.cell_to_index_map,2));
    numCells(i) = length(idcs{i,1});
end


%% Calculate probability of same encoding

for n=1:size(this_d,1)
    for k=1:numTestGroups
        Psame_all{n,k} = nan(d_info.numDays,d_info.numDays);
        Psame_pos{n,k} = nan(d_info.numDays,d_info.numDays);
        Psame_neg{n,k} = nan(d_info.numDays,d_info.numDays);
        for i=1:d_info.numDays
            for j=1:d_info.numDays

                % get basic structs
                this_nem_i = this_d{n,i}.nem_all_cmpr;
                this_nem_j = this_d{n,j}.nem_all_cmpr;
                % MB20220713: changed that from: temp = intersect(find(this_nem_i.sigM(trck.cell_to_index_map(idcs,i))==1),find(this_nem_j.sigM(trck.cell_to_index_map(idcs,i))==1));
                temp = intersect(find(this_nem_i.sigM(trck{n}.cell_to_index_map(idcs{n},i))==1),find(this_nem_j.sigM(trck{n}.cell_to_index_map(idcs{n},j))==1));
                % this is basically just the intersection of cells that pass model across two days
                
                % proportion of fully tracked cells encoding this feature in both block i and block j
                temp_i = this_nem_i.sigM_sigTG{k}(trck{n}.cell_to_index_map(idcs{n},i));
                temp_j = this_nem_j.sigM_sigTG{k}(trck{n}.cell_to_index_map(idcs{n},j));
                Psame_all{n,k}(i,j) = nansum(temp_j(find(temp_i))) / length(temp); %numCells; %length(find(temp_i));

                % proportion of fully tracked cells encoding this feature in both block i and block j
                temp_i = this_nem_i.sigM_sigTG_pos{k}(trck{n}.cell_to_index_map(idcs{n},i));
                temp_j = this_nem_j.sigM_sigTG_pos{k}(trck{n}.cell_to_index_map(idcs{n},j));
                Psame_pos{n,k}(i,j) = nansum(temp_j(find(temp_i))) / length(temp); %numCells; %length(find(temp_i));

                % proportion of fully tracked cells encoding this feature in both block i and block j
                temp_i = this_nem_i.sigM_sigTG_neg{k}(trck{n}.cell_to_index_map(idcs{n},i));
                temp_j = this_nem_j.sigM_sigTG_neg{k}(trck{n}.cell_to_index_map(idcs{n},j));
                Psame_neg{n,k}(i,j) = nansum(temp_j(find(temp_i))) / length(temp); %numCells; %length(find(temp_i));
            end
        end
    end
end


%% Compare Psame_pos between within and across switch days

Psame_pos_withinSwitch = nan(size(this_d,1),numTestGroups);
Psame_pos_acrossSwitches = nan(size(this_d,1),numTestGroups);
for n=1:size(this_d,1)
    for k=1:numTestGroups   
        Psame_pos_withinSwitch(n,k) = nanmean([Psame_pos{n,k}(2,3),Psame_pos{n,k}(4,5)]);
        Psame_pos_acrossSwitches(n,k) = nanmean([Psame_pos{n,k}(2,1),Psame_pos{n,k}(4,3)]);
    end
end


%%

nrows = 2; ncols = 4;
default_figure([20,0.5,20,9.9]);


k = 1;
subplot(nrows,ncols,1)

h=bar(diag(nanmean([Psame_pos_withinSwitch(:,k),Psame_pos_acrossSwitches(:,k)],1)*100),'stacked','BaseValue',0,'FaceColor',p.col.gray,'EdgeColor','none');
hold on
for i=1:size(this_d,1)
    plot([Psame_pos_withinSwitch(i,k),Psame_pos_acrossSwitches(i,k)]*100,'k-');
end
scatter(1*ones(size(this_d,1),1),Psame_pos_withinSwitch(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
scatter(2*ones(size(this_d,1),1),Psame_pos_acrossSwitches(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');

% ylim([0,100])
ytickformat('percentage')
xticks(1:2)
xticklabels({'within','across'})
ylabel('Proportion of fully tracked neurons')
set(gca,'box','off')
title([nem_ref.testGroup.label(k)])


k = 10;
subplot(nrows,ncols,5)

h=bar(diag(nanmean([Psame_pos_withinSwitch(:,k),Psame_pos_acrossSwitches(:,k)],1)*100),'stacked','BaseValue',0,'FaceColor',p.col.gray,'EdgeColor','none');
hold on
for i=1:size(this_d,1)
    plot([Psame_pos_withinSwitch(i,k),Psame_pos_acrossSwitches(i,k)]*100,'k-');
end
scatter(1*ones(size(this_d,1),1),Psame_pos_withinSwitch(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
scatter(2*ones(size(this_d,1),1),Psame_pos_acrossSwitches(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');

% ylim([0,100])
ytickformat('percentage')
xticks(1:2)
xticklabels({'within','across'})
ylabel('Proportion of fully tracked neurons')
set(gca,'box','off')
title([nem_ref.testGroup.label(k)])


k = 2;
subplot(nrows,ncols,2)

h=bar(diag(nanmean([Psame_pos_withinSwitch(:,k),Psame_pos_acrossSwitches(:,k)],1)*100),'stacked','BaseValue',0,'FaceColor',p.col.gray,'EdgeColor','none');
hold on
for i=1:size(this_d,1)
    plot([Psame_pos_withinSwitch(i,k),Psame_pos_acrossSwitches(i,k)]*100,'k-');
end
scatter(1*ones(size(this_d,1),1),Psame_pos_withinSwitch(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
scatter(2*ones(size(this_d,1),1),Psame_pos_acrossSwitches(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');

% ylim([0,100])
ytickformat('percentage')
xticks(1:2)
xticklabels({'within','across'})
ylabel('Proportion of fully tracked neurons')
set(gca,'box','off')
title([nem_ref.testGroup.label(k)])


k = 11;
subplot(nrows,ncols,6)

h=bar(diag(nanmean([Psame_pos_withinSwitch(:,k),Psame_pos_acrossSwitches(:,k)],1)*100),'stacked','BaseValue',0,'FaceColor',p.col.gray,'EdgeColor','none');
hold on
for i=1:size(this_d,1)
    plot([Psame_pos_withinSwitch(i,k),Psame_pos_acrossSwitches(i,k)]*100,'k-');
end
scatter(1*ones(size(this_d,1),1),Psame_pos_withinSwitch(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
scatter(2*ones(size(this_d,1),1),Psame_pos_acrossSwitches(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');

% ylim([0,100])
ytickformat('percentage')
xticks(1:2)
xticklabels({'within','across'})
ylabel('Proportion of fully tracked neurons')
set(gca,'box','off')
title([nem_ref.testGroup.label(k)])


k = 3;
subplot(nrows,ncols,3)

h=bar(diag(nanmean([Psame_pos_withinSwitch(:,k),Psame_pos_acrossSwitches(:,k)],1)*100),'stacked','BaseValue',0,'FaceColor',p.col.gray,'EdgeColor','none');
hold on
for i=1:size(this_d,1)
    plot([Psame_pos_withinSwitch(i,k),Psame_pos_acrossSwitches(i,k)]*100,'k-');
end
scatter(1*ones(size(this_d,1),1),Psame_pos_withinSwitch(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
scatter(2*ones(size(this_d,1),1),Psame_pos_acrossSwitches(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');

% ylim([0,100])
ytickformat('percentage')
xticks(1:2)
xticklabels({'within','across'})
ylabel('Proportion of fully tracked neurons')
set(gca,'box','off')
title([nem_ref.testGroup.label(k)])


k = 12;
subplot(nrows,ncols,7)

h=bar(diag(nanmean([Psame_pos_withinSwitch(:,k),Psame_pos_acrossSwitches(:,k)],1)*100),'stacked','BaseValue',0,'FaceColor',p.col.gray,'EdgeColor','none');
hold on
for i=1:size(this_d,1)
    plot([Psame_pos_withinSwitch(i,k),Psame_pos_acrossSwitches(i,k)]*100,'k-');
end
scatter(1*ones(size(this_d,1),1),Psame_pos_withinSwitch(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
scatter(2*ones(size(this_d,1),1),Psame_pos_acrossSwitches(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');

% ylim([0,100])
ytickformat('percentage')
xticks(1:2)
xticklabels({'within','across'})
ylabel('Proportion of fully tracked neurons')
set(gca,'box','off')
title([nem_ref.testGroup.label(k)])


k = 17;
subplot(nrows,ncols,4)

h=bar(diag(nanmean([Psame_pos_withinSwitch(:,k),Psame_pos_acrossSwitches(:,k)],1)*100),'stacked','BaseValue',0,'FaceColor',p.col.gray,'EdgeColor','none');
hold on
for i=1:size(this_d,1)
    plot([Psame_pos_withinSwitch(i,k),Psame_pos_acrossSwitches(i,k)]*100,'k-');
end
scatter(1*ones(size(this_d,1),1),Psame_pos_withinSwitch(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
scatter(2*ones(size(this_d,1),1),Psame_pos_acrossSwitches(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');

% ylim([0,100])
ytickformat('percentage')
xticks(1:2)
xticklabels({'within','across'})
ylabel('Proportion of fully tracked neurons')
set(gca,'box','off')
title([nem_ref.testGroup.label(k)])


k = 16;
subplot(nrows,ncols,8)

h=bar(diag(nanmean([Psame_pos_withinSwitch(:,k),Psame_pos_acrossSwitches(:,k)],1)*100),'stacked','BaseValue',0,'FaceColor',p.col.gray,'EdgeColor','none');
hold on
for i=1:size(this_d,1)
    plot([Psame_pos_withinSwitch(i,k),Psame_pos_acrossSwitches(i,k)]*100,'k-');
end
scatter(1*ones(size(this_d,1),1),Psame_pos_withinSwitch(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');
scatter(2*ones(size(this_d,1),1),Psame_pos_acrossSwitches(:,k)*100,'MarkerEdgeColor','k','MarkerFaceColor','none');

% ylim([0,100])
ytickformat('percentage')
xticks(1:2)
xticklabels({'within','across'})
ylabel('Proportion of fully tracked neurons')
set(gca,'box','off')
title([nem_ref.testGroup.label(k)])


%%