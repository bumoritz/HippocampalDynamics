%% Data import
this_d = {};

% this_animal = 3; % Biontech
% this_d{1,1} = d{this_animal,1}; this_d{1,2} = d{this_animal,2}; this_d{1,3} = d{this_animal,3}; this_d{1,4} = d{this_animal,4}; this_d{1,5} = d{this_animal,5};
% trck = load('D:\SniffinHippo\Repo\Biontech\Biontech_combined\Biontech_combined_trck_15_rigid.mat'); trck = trck.trck;

% this_animal = 8; % Arasaka
% this_d{1,1} = d{this_animal,1}; this_d{1,2} = d{this_animal,2}; this_d{1,3} = d{this_animal,3}; this_d{1,4} = d{this_animal,4}; this_d{1,5} = d{this_animal,5};
% trck = load('D:\SniffinHippo\Repo\Arasaka\Arasaka_combined\Arasaka_combined_trck_15_rigid.mat'); trck = trck.trck;

% this_animal = 31; % Maria
% this_d{1,1} = d{this_animal,1}; this_d{1,2} = d{this_animal,2}; this_d{1,3} = d{this_animal,3}; this_d{1,4} = d{this_animal,4}; this_d{1,5} = d{this_animal,5};
% trck = load('D:\SniffinHippo\Repo\Maria\Maria_combined\Maria_combined_trck_15_nonrigid.mat'); trck = trck.trck;

this_animal = 33; % William
this_d{1,1} = d{this_animal,1}; this_d{1,2} = d{this_animal,2}; this_d{1,3} = d{this_animal,3}; this_d{1,4} = d{this_animal,4}; this_d{1,5} = d{this_animal,5};
trck = load('D:\SniffinHippo\Repo\William\William_combined\William_combined_trck_15_nonrigid.mat'); trck = trck.trck;


%% Preparations

idcs = find(sum(trck.cell_to_index_map>0,2)==size(trck.cell_to_index_map,2));
numCells = length(idcs);
nem_ref = this_d{1}.nem_all_cmpr;
numTestGroups = length(nem_ref.testGroup.label);


%% Calculate probability of same encoding

for k=1:numTestGroups
    Psame_all{k} = nan(d_info.numDays,d_info.numDays);
    Psame_pos{k} = nan(d_info.numDays,d_info.numDays);
    Psame_neg{k} = nan(d_info.numDays,d_info.numDays);
    for i=1:d_info.numDays
        for j=1:d_info.numDays

            % get basic structs
            this_nem_i = this_d{i}.nem_all_cmpr;
            this_nem_j = this_d{j}.nem_all_cmpr;
            % MB20220713: changed that from: temp = intersect(find(this_nem_i.sigM(trck.cell_to_index_map(idcs,i))==1),find(this_nem_j.sigM(trck.cell_to_index_map(idcs,i))==1));
            temp = intersect(find(this_nem_i.sigM(trck.cell_to_index_map(idcs,i))==1),find(this_nem_j.sigM(trck.cell_to_index_map(idcs,j))==1));

            % proportion of fully tracked cells encoding this feature in both block i and block j
            temp_i = this_nem_i.sigM_sigTG{k}(trck.cell_to_index_map(idcs,i));
            temp_j = this_nem_j.sigM_sigTG{k}(trck.cell_to_index_map(idcs,j));
            Psame_all{k}(i,j) = nansum(temp_j(find(temp_i))) / length(temp); %numCells; %length(find(temp_i));
            
            % proportion of fully tracked cells encoding this feature in both block i and block j
            temp_i = this_nem_i.sigM_sigTG_pos{k}(trck.cell_to_index_map(idcs,i));
            temp_j = this_nem_j.sigM_sigTG_pos{k}(trck.cell_to_index_map(idcs,j));
            Psame_pos{k}(i,j) = nansum(temp_j(find(temp_i))) / length(temp); %numCells; %length(find(temp_i));
            
            % proportion of fully tracked cells encoding this feature in both block i and block j
            temp_i = this_nem_i.sigM_sigTG_neg{k}(trck.cell_to_index_map(idcs,i));
            temp_j = this_nem_j.sigM_sigTG_neg{k}(trck.cell_to_index_map(idcs,j));
            Psame_neg{k}(i,j) = nansum(temp_j(find(temp_i))) / length(temp); %numCells; %length(find(temp_i));
        end
    end
end

% problem if divided by number of fully tracked neurons: more significant model fits on some days than others
% maybe try dividing by number of fully tracked neurons, that have significant model fits on both days

%% Calculate number of cells with same encoding

for k=1:numTestGroups
    Nsame_all{k} = nan(d_info.numDays,d_info.numDays);
    Nsame_pos{k} = nan(d_info.numDays,d_info.numDays);
    Nsame_neg{k} = nan(d_info.numDays,d_info.numDays);
    for i=1:d_info.numDays
        for j=1:d_info.numDays

            % get basic structs
            this_nem_i = this_d{i}.nem_all_cmpr;
            this_nem_j = this_d{j}.nem_all_cmpr;

            % proportion of fully tracked cells encoding this feature in both block i and block j
            temp_i = this_nem_i.sigM_sigTG{k}(trck.cell_to_index_map(idcs,i));
            temp_j = this_nem_j.sigM_sigTG{k}(trck.cell_to_index_map(idcs,j));
            Nsame_all{k}(i,j) = nansum(temp_j(find(temp_i)));
            
            % proportion of fully tracked cells encoding this feature in both block i and block j
            temp_i = this_nem_i.sigM_sigTG_pos{k}(trck.cell_to_index_map(idcs,i));
            temp_j = this_nem_j.sigM_sigTG_pos{k}(trck.cell_to_index_map(idcs,j));
            Nsame_pos{k}(i,j) = nansum(temp_j(find(temp_i)));
            
            % proportion of fully tracked cells encoding this feature in both block i and block j
            temp_i = this_nem_i.sigM_sigTG_neg{k}(trck.cell_to_index_map(idcs,i));
            temp_j = this_nem_j.sigM_sigTG_neg{k}(trck.cell_to_index_map(idcs,j));
            Nsame_neg{k}(i,j) = nansum(temp_j(find(temp_i)));      
        end
    end
end


%% Figure

nrows = 3;
ncols = 7;
default_figure([20,0.5,20,9.9]);

for k=1:numTestGroups
    subplot(nrows,ncols,k)

    imagesc(Psame_pos{k})
    %imagesc(Nsame_pos{k})
    hold on
    xline(1+0.5,'LineStyle','-');
    xline(3+0.5,'LineStyle','-');
    yline(1+0.5,'LineStyle','-');
    yline(3+0.5,'LineStyle','-');

    set(gca,'CLim',[0,inf]); %set(gca,'CLim',[0,inf]); % 0.35
    colormap(gca,flipud(bone));
    daspect([1,1,1])
    colorbar;
    title([nem_ref.testGroup.label{k},' (',num2str(Nsame_pos{k}(1,1)),',',num2str(Nsame_pos{k}(2,2)),',',num2str(Nsame_pos{k}(3,3)),',',num2str(Nsame_pos{k}(4,4)),',',num2str(Nsame_pos{k}(5,5)),')'])
end


%% Figure

nrows = 2;
ncols = 4;
default_figure([20,0.5,20,9.9]);

for k=1:numTestGroups
    temp = [1,2,3,17,10,11,12,16];
    if ismember(k,temp)
        [~,temp2] = find(temp==k);
        subplot(nrows,ncols,temp2)

        imagesc(Psame_pos{k})
        %imagesc(Nsame_pos{k})
        hold on
        xline(1+0.5,'LineStyle','-');
        xline(3+0.5,'LineStyle','-');
        yline(1+0.5,'LineStyle','-');
        yline(3+0.5,'LineStyle','-');

        set(gca,'CLim',[0,inf]); %set(gca,'CLim',[0,inf]); % 0.35
        colormap(gca,flipud(bone));
        daspect([1,1,1])
        h=colorbar;
        ylabel(h,'Proportion of fully tracked neurons','FontSize',10,'Rotation',270);
        h.Label.Position(1) = 3;
        xlabel('Experiment day')
        ylabel('Experiment day')
        %title([nem_ref.testGroup.label{k},' (',num2str(Nsame_pos{k}(1,1)),',',num2str(Nsame_pos{k}(2,2)),',',num2str(Nsame_pos{k}(3,3)),',',num2str(Nsame_pos{k}(4,4)),',',num2str(Nsame_pos{k}(5,5)),')'])
        title(nem_ref.testGroup.label{k})
    end
end

%%

% likelihood of encoding 1st odour if:
% - if encoding 1st odour before
% - if endoing 2nd odour before

% for both switches: day before switch -> day of switch











