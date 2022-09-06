%% Data import
this_d = {};

this_animal = 8; % Arasaka
this_d{1,1} = d{this_animal,1}; this_d{1,2} = d{this_animal,2}; this_d{1,3} = d{this_animal,3}; this_d{1,4} = d{this_animal,4}; this_d{1,5} = d{this_animal,5};
trck = load('D:\SniffinHippo\Repo\Arasaka\Arasaka_combined\Arasaka_combined_trck_15_rigid.mat'); trck = trck.trck;
 
% this_animal = 31; % Maria
% this_d{1,1} = d{this_animal,1}; this_d{1,2} = d{this_animal,2}; this_d{1,3} = d{this_animal,3}; this_d{1,4} = d{this_animal,4}; this_d{1,5} = d{this_animal,5};
% trck = load('D:\SniffinHippo\Repo\Maria\Maria_combined\Maria_combined_trck_15_nonrigid.mat'); trck = trck.trck;

% this_animal = 33; % William
% this_d{1,1} = d{this_animal,1}; this_d{1,2} = d{this_animal,2}; this_d{1,3} = d{this_animal,3}; this_d{1,4} = d{this_animal,4}; this_d{1,5} = d{this_animal,5};
% trck = load('D:\SniffinHippo\Repo\William\William_combined\William_combined_trck_15_nonrigid.mat'); trck = trck.trck;
% 

%% Preparations

idcs = find(sum(trck.cell_to_index_map>0,2)==size(trck.cell_to_index_map,2));
numCells = length(idcs);
nem_ref = this_d{1}.nem_100t_cmpr{1};
numTestGroups = length(nem_ref.testGroup.label);


%% Calculate probability of same encoding

blocks = [5,8,5,8,5];
numBlocks = sum(blocks);
for k=1:numTestGroups
    Psame_all{k} = nan(numBlocks,numBlocks);
    Psame_pos{k} = nan(numBlocks,numBlocks);
    Psame_neg{k} = nan(numBlocks,numBlocks);
    for i=1:numBlocks
        for j=1:numBlocks
            temp = cumsum(blocks);
            m_i = min(find(temp-i>=0)); % day
            n_i = blocks(m_i)-(temp(m_i)-i); % block
            m_j = min(find(temp-j>=0)); % day
            n_j = blocks(m_j)-(temp(m_j)-j); % block
            
            % get basic structs
            this_nem_i = this_d{m_i}.nem_100t_cmpr{n_i};
            this_nem_j = this_d{m_j}.nem_100t_cmpr{n_j};

            % proportion of fully tracked cells encoding this feature in both block i and block j
            temp_i = this_nem_i.sigM_sigTG{k}(trck.cell_to_index_map(idcs,m_i));
            temp_j = this_nem_j.sigM_sigTG{k}(trck.cell_to_index_map(idcs,m_j));
            Psame_all{k}(i,j) = nansum(temp_j(find(temp_i))) / numCells; %length(find(temp_i));
            
            % proportion of fully tracked cells encoding this feature in both block i and block j
            temp_i = this_nem_i.sigM_sigTG_pos{k}(trck.cell_to_index_map(idcs,m_i));
            temp_j = this_nem_j.sigM_sigTG_pos{k}(trck.cell_to_index_map(idcs,m_j));
            Psame_pos{k}(i,j) = nansum(temp_j(find(temp_i))) / numCells; %length(find(temp_i));
            
            % proportion of fully tracked cells encoding this feature in both block i and block j
            temp_i = this_nem_i.sigM_sigTG_neg{k}(trck.cell_to_index_map(idcs,m_i));
            temp_j = this_nem_j.sigM_sigTG_neg{k}(trck.cell_to_index_map(idcs,m_j));
            Psame_neg{k}(i,j) = nansum(temp_j(find(temp_i))) / numCells; %length(find(temp_i));
        end
    end
end


%% Calculate number of cells with same encoding

blocks = [5,8,5,8,5];
numBlocks = sum(blocks);
for k=1:numTestGroups
    Nsame_all{k} = nan(numBlocks,numBlocks);
    Nsame_pos{k} = nan(numBlocks,numBlocks);
    Nsame_neg{k} = nan(numBlocks,numBlocks);
    for i=1:numBlocks
        for j=1:numBlocks
            temp = cumsum(blocks);
            m_i = min(find(temp-i>=0)); % day
            n_i = blocks(m_i)-(temp(m_i)-i); % block
            m_j = min(find(temp-j>=0)); % day
            n_j = blocks(m_j)-(temp(m_j)-j); % block
            
            % get basic structs
            this_nem_i = this_d{m_i}.nem_100t_cmpr{n_i};
            this_nem_j = this_d{m_j}.nem_100t_cmpr{n_j};

            % proportion of fully tracked cells encoding this feature in both block i and block j
            temp_i = this_nem_i.sigM_sigTG{k}(trck.cell_to_index_map(idcs,m_i));
            temp_j = this_nem_j.sigM_sigTG{k}(trck.cell_to_index_map(idcs,m_j));
            Nsame_all{k}(i,j) = nansum(temp_j(find(temp_i)));
            
            % proportion of fully tracked cells encoding this feature in both block i and block j
            temp_i = this_nem_i.sigM_sigTG_pos{k}(trck.cell_to_index_map(idcs,m_i));
            temp_j = this_nem_j.sigM_sigTG_pos{k}(trck.cell_to_index_map(idcs,m_j));
            Nsame_pos{k}(i,j) = nansum(temp_j(find(temp_i)));
            
            % proportion of fully tracked cells encoding this feature in both block i and block j
            temp_i = this_nem_i.sigM_sigTG_neg{k}(trck.cell_to_index_map(idcs,m_i));
            temp_j = this_nem_j.sigM_sigTG_neg{k}(trck.cell_to_index_map(idcs,m_j));
            Nsame_neg{k}(i,j) = nansum(temp_j(find(temp_i)));
        end
    end
end


%% Figure

nrows = 3;
ncols = 7;
% default_figure([20,0.5,20,9.9]);

for k=1:numTestGroups

    if ismember(k,[1,2,3,10,11,12,16,17])
        figure; %subplot(nrows,ncols,k)

        imagesc(Psame_pos{k})
        %imagesc(Nsame_pos{k})
        hold on
        temp = cumsum(blocks);
        xline(temp(1)+0.5,'LineStyle','-');
        xline(temp(2)+0.5,'LineStyle',':');
        xline(temp(3)+0.5,'LineStyle','-');
        xline(temp(4)+0.5,'LineStyle',':');
        yline(temp(1)+0.5,'LineStyle','-');
        yline(temp(2)+0.5,'LineStyle',':');
        yline(temp(3)+0.5,'LineStyle','-');
        yline(temp(4)+0.5,'LineStyle',':');

        set(gca,'CLim',[0,inf]); %set(gca,'CLim',[0,inf]);%set(gca,'CLim',[0,1]); 0.25
        colormap(gca,flipud(bone));
        daspect([1,1,1])
        h=colorbar;
        ylabel(h,'Proportion of fully tracked neurons','FontSize',10,'Rotation',270);
        h.Label.Position(1) = 4;
        xlabel('Trial block (100 trials each)')
        ylabel('Trial block (100 trials each)')
        title([nem_ref.testGroup.label{k}])
        %title([nem_ref.testGroup.label{k},' (',num2str(Nsame_pos{k}(1,1)),',',num2str(Nsame_pos{k}(temp(1)+1,temp(1)+1)),',',num2str(Nsame_pos{k}(temp(2)+1,temp(2)+1)),',',num2str(Nsame_pos{k}(temp(3)+1,temp(3)+1)),',',num2str(Nsame_pos{k}(temp(4)+1,temp(4)+1)),')'])
    end
end





