%% SfN

% d_info = selectDataset(d_info,'-g78-d2345-e01-r01-p01-l01-i00',sheet,path,ops);


%% Preparations

blocks = [5,8,5,8,5];
numBlocks = sum(blocks);

binEdges = [0,1.5,3,5.3];
numBins = length(binEdges)-1;
binLabels = {'Early peak (0 s - 1.5 s)','Middle peak (1.5 s - 3 s)','Late peak (3 s - 5.3 s)'};

showIndividuals = false;
showAverages = true;
showCorrelation = true;

tng_struct = 'tng_100t_stimVersion';
warp_struct = 'warp_100t_stimVersion';

firingField_metric = 'peak'; %'com'; 'peak'


%% Get performance metrics by block and animal

correct = nan(d_info.numAnimals,numBlocks); % [animal, block]
correct_catch = nan(d_info.numAnimals,numBlocks); % [animal, block]
hitRate = nan(d_info.numAnimals,numBlocks); % [animal, block]
crRate = nan(d_info.numAnimals,numBlocks); % [animal, block]
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if j==1
            temp0=0;
        else
            temp0 = cumsum(blocks);
            temp0 = temp0(j-1);
        end
        if ~isempty(d{i,j})
            for k=1:blocks(j)
                this_block = temp0+k;
                these_blocks_20t = (k-1)*5+1:k*5;
                
                % get performance metrics
                correct(i,this_block) = nanmean(d{i,j}.perf.blocks_general.correct(these_blocks_20t));
                correct_catch(i,this_block) = nanmean(d{i,j}.perf.blocks_general.correct_catch(these_blocks_20t));
                hitRate(i,this_block) = nanmean(d{i,j}.perf.blocks_general.H(these_blocks_20t));
                crRate(i,this_block) = nanmean(d{i,j}.perf.blocks_general.CR(these_blocks_20t));
            end
        end
    end
end


%% Get data

% for each 100t block across all 5 days
seqCellsByTime_num = nan(numBins,numBlocks,d_info.numAnimals); % [peak time bin, block, animal]
seqCellsByTime_prop = nan(numBins,numBlocks,d_info.numAnimals); % [peak time bin, block, animal]
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if j==1
            temp0=0;
        else
            temp0 = cumsum(blocks);
            temp0 = temp0(j-1);
        end
        if isfield(d{i,j},tng_struct) %~isempty(d{i,j})
            for k=1:length(d{i,j}.(tng_struct))
                this_block = temp0+k;

                % get idcs
                this_numCells = length(find(d{i,j}.(tng_struct){k}.prop.iscell==1));
                these_idcs_A = find(d{i,j}.(tng_struct){k}.passed_catch.AW.Acatchonly==1);
                these_idcs_X = find(d{i,j}.(tng_struct){k}.passed_catch.AW.Xcatchonly==1);

                if strcmp(firingField_metric,'peak')
                    % get peak locations
                    these_peakTimes_A = d{i,j}.(tng_struct){k}.firingField.Acatch_AW.peakLocation_s(these_idcs_A);
                    these_peakTimes_X = d{i,j}.(tng_struct){k}.firingField.Xcatch_AW.peakLocation_s(these_idcs_X);
                    these_peakTimes = [these_peakTimes_A; these_peakTimes_X];
                elseif strcmp(firingField_metric,'com')
                    % get com locations
                    these_peakTimes_A = d{i,j}.(tng_struct){k}.firingField.Acatch_AW.comLocation_s(these_idcs_A);
                    these_peakTimes_X = d{i,j}.(tng_struct){k}.firingField.Xcatch_AW.comLocation_s(these_idcs_X);
                    these_peakTimes = [these_peakTimes_A; these_peakTimes_X];
                end

                % bin peak locations
%                 temp2 = sort(these_peakTimes);
%                 these_binEdges = [0,temp2(round(length(temp2)/3)),temp2(2*round(length(temp2)/3)),5.3];
%                 temp = discretize(these_peakTimes,these_binEdges);
                temp = discretize(these_peakTimes,binEdges);
                for n=1:numBins
                    seqCellsByTime_num(n,this_block,i) = nansum(temp==n);
                    seqCellsByTime_prop(n,this_block,i) = nansum(temp==n) / this_numCells;
                end
            end
        end
    end
end
sequenceCells_num_avg = nanmean(seqCellsByTime_num,3);
sequenceCells_prop_avg = nanmean(seqCellsByTime_prop,3);

% for each 100t block within switch
blocks_switch = [8,5];
numBlocks_switch = sum(blocks_switch);
seqCellsByTime_switch_num = nan(numBins,numBlocks_switch,d_info.numAnimals); % [peak time bin, block, animal]
seqCellsByTime_switch_prop = nan(numBins,numBlocks_switch,d_info.numAnimals); % [peak time bin, block, animal]
for i=1:d_info.numAnimals
    for k=1:numBlocks_switch
        seqCellsByTime_switch_num(:,k,i) = nanmean([seqCellsByTime_num(:,blocks(1)+k,i),seqCellsByTime_num(:,sum(blocks(1:3))+k,i)],2);
        seqCellsByTime_switch_prop(:,k,i) = nanmean([seqCellsByTime_prop(:,blocks(1)+k,i),seqCellsByTime_prop(:,sum(blocks(1:3))+k,i)],2);
    end
end
seqCellsByTime_switch_num_avg = nanmean(seqCellsByTime_switch_num,3);
seqCellsByTime_switch_prop_avg = nanmean(seqCellsByTime_switch_prop,3);


%% Get sequence warp by block and animal

% for each 100t block across all 5 days
warp_exp1_a = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_exp1_b = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_power2_a = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_power2_b = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_power2_c = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_exp1const_a = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_exp1const_b = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_exp1const_c = nan(d_info.numAnimals,numBlocks); % [animal, block]
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if j==1
            temp0=0;
        else
            temp0 = cumsum(blocks);
            temp0 = temp0(j-1);
        end
        if isfield(d{i,j},warp_struct) %~isempty(d{i,j})
            for k=1:length(d{i,j}.(warp_struct))
                this_block = temp0+k;
                
                try
                    warp_exp1_a(i,this_block) = d{i,j}.(warp_struct){k}.(firingField_metric).exp1.a;
                    warp_exp1_b(i,this_block) = d{i,j}.(warp_struct){k}.(firingField_metric).exp1.b;
                    warp_power2_a(i,this_block) = d{i,j}.(warp_struct){k}.(firingField_metric).power2.a;
                    warp_power2_b(i,this_block) = d{i,j}.(warp_struct){k}.(firingField_metric).power2.b;
                    warp_power2_c(i,this_block) = d{i,j}.(warp_struct){k}.(firingField_metric).power2.c;
                catch
                    warp_exp1_a(i,this_block) = d{i,j}.(warp_struct){k}.exp1.a;
                    warp_exp1_b(i,this_block) = d{i,j}.(warp_struct){k}.exp1.b;
                    warp_power2_a(i,this_block) = d{i,j}.(warp_struct){k}.power2.a;
                    warp_power2_b(i,this_block) = d{i,j}.(warp_struct){k}.power2.b;
                    warp_power2_c(i,this_block) = d{i,j}.(warp_struct){k}.power2.c;
                end
%                 warp_exp1const_a(i,this_block) = d{i,j}.(warp_struct){k}.exp1const.a;
%                 warp_exp1const_b(i,this_block) = d{i,j}.(warp_struct){k}.exp1const.b;
%                 warp_exp1const_c(i,this_block) = d{i,j}.(warp_struct){k}.exp1const.c;
            end
        end
    end
end



%% Simplify correct_catch and warp_exp1 data

correct_catch_switch1 = correct_catch(:,6:18);
correct_catch_switch2 = correct_catch(:,19:31);
correct_catch_all = [correct_catch_switch1;correct_catch_switch1];
correct_catch_switch1_seq = correct_catch_switch1(find(d_info.group==7),:);
correct_catch_switch1_ctrl = correct_catch_switch1(find(d_info.group==8),:);
correct_catch_switch2_seq = correct_catch_switch2(find(d_info.group==8),:);
correct_catch_switch2_ctrl = correct_catch_switch2(find(d_info.group==7),:);
correct_catch_seq = [correct_catch_switch1_seq;correct_catch_switch2_seq];
correct_catch_ctrl = [correct_catch_switch1_ctrl;correct_catch_switch2_ctrl];

warp_exp1_a_switch1 = warp_exp1_a(:,6:18);
warp_exp1_a_switch2 = warp_exp1_a(:,19:31);
warp_exp1_a_all = [warp_exp1_a_switch1;warp_exp1_a_switch1];
warp_exp1_a_switch1_seq = warp_exp1_a_switch1(find(d_info.group==7),:);
warp_exp1_a_switch1_ctrl = warp_exp1_a_switch1(find(d_info.group==8),:);
warp_exp1_a_switch2_seq = warp_exp1_a_switch2(find(d_info.group==8),:);
warp_exp1_a_switch2_ctrl = warp_exp1_a_switch2(find(d_info.group==7),:);
warp_exp1_a_seq = [warp_exp1_a_switch1_seq;warp_exp1_a_switch2_seq];
warp_exp1_a_ctrl = [warp_exp1_a_switch1_ctrl;warp_exp1_a_switch2_ctrl];

warp_exp1_b_switch1 = warp_exp1_b(:,6:18);
warp_exp1_b_switch2 = warp_exp1_b(:,19:31);
warp_exp1_b_all = [warp_exp1_b_switch1;warp_exp1_b_switch1];
warp_exp1_b_switch1_seq = warp_exp1_b_switch1(find(d_info.group==7),:);
warp_exp1_b_switch1_ctrl = warp_exp1_b_switch1(find(d_info.group==8),:);
warp_exp1_b_switch2_seq = warp_exp1_b_switch2(find(d_info.group==8),:);
warp_exp1_b_switch2_ctrl = warp_exp1_b_switch2(find(d_info.group==7),:);
warp_exp1_b_seq = [warp_exp1_b_switch1_seq;warp_exp1_b_switch2_seq];
warp_exp1_b_ctrl = [warp_exp1_b_switch1_ctrl;warp_exp1_b_switch2_ctrl];


%% SUMMARY FIGURE - discarded

% select sessions
these_animals = 'stim'; 
these_sessions = ~d_info.excl & d_info.presponsive==1; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1; & d_info.running==1
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

% select data
this_data = warp_exp1_a; %seqCellsByTime_prop;

nrows = 1; ncols = 5; m=0;
F = default_figure([0,3.5,18,2.5]);
m = m+1; subplot(nrows,ncols,m); hold on;

temp_data = [];
for idx=1:d_info.numAnimals
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
        temp = this_data(idx,6:18);
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
        temp = this_data(idx,19:31);
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
end
temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.ctrl);
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r_ctrl,this_corr_p_ctrl] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete');  %fitLine(temp(:),temp_data(:),'k');
this_corr_n_ctrl = size(temp_data,1);

temp_data = [];
for idx=1:d_info.numAnimals
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
        temp = this_data(idx,6:18);
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
        temp = this_data(idx,19:31);
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
end
temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.seq);
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r_seq,this_corr_p_seq] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete'); %fitLine(temp(:),temp_data(:),'k');
this_corr_n_seq = size(temp_data,1);

xline(blocks_switch(1)+0.5,'k:');
xlim([0,numBlocks_switch])
ylim([0,0.25])
xlabel('Trial block (in blocks of 100 trials)')
ylim([0,0.12])
yticks([0,0.12])
xticks([1,8,13])
ylabel(['Coefficient a',newline,'(initial proportion)'])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['seq: n = ',num2str(this_corr_n_seq),', \rho = ',num2str(this_corr_r_seq,2),', p = ',num2str(this_corr_p_seq,2),newline,...
    'ctrl: n = ',num2str(this_corr_n_ctrl),', \rho = ',num2str(this_corr_r_ctrl,2),', p = ',num2str(this_corr_p_ctrl,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);

% select data
this_data = warp_exp1_b; %seqCellsByTime_prop;
m = m+1; subplot(nrows,ncols,m); hold on;

temp_data = [];
for idx=1:d_info.numAnimals
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
        temp = this_data(idx,6:18);
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
        temp = this_data(idx,19:31);
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
end
temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.ctrl);
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r_ctrl,this_corr_p_ctrl] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete');  %fitLine(temp(:),temp_data(:),'k');
this_corr_n_ctrl = size(temp_data,1);

temp_data = [];
for idx=1:d_info.numAnimals
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
        temp = this_data(idx,6:18);
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
        temp = this_data(idx,19:31);
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
end
temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.seq);
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r_seq,this_corr_p_seq] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete'); %fitLine(temp(:),temp_data(:),'k');
this_corr_n_seq = size(temp_data,1);
xline(blocks_switch(1)+0.5,'k:');
xlim([0,numBlocks_switch])
ylim([-0.7,0])
yticks([-0.7,0])
xticks([1,8,13])
ylabel(['Coefficient b',newline,'(decay rate, s-1)'])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['seq: n = ',num2str(this_corr_n_seq),', \rho = ',num2str(this_corr_r_seq,2),', p = ',num2str(this_corr_p_seq,2),newline,...
    'ctrl: n = ',num2str(this_corr_n_ctrl),', \rho = ',num2str(this_corr_r_ctrl,2),', p = ',num2str(this_corr_p_ctrl,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);

% select data
this_data = seqCellsByTime_prop; %seqCellsByTime_prop;

% select sessions
these_animals = 'stim';
these_sessions = ~d_info.excl & d_info.presponsive==1; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1;

% preparations
%these_rgbs = discretisedColourMap('winter',true,numBins);
this_n_lower = min(sum(these_sessions(:,2:5),1));
this_n_upper = max(sum(these_sessions(:,2:5),1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

% nrows = 1; ncols = numBins; m=0;
% F = default_figure([-20,0.5,15,5]);

early_middle_late_seq = {};
early_middle_late_seq_corr = {};
early_middle_late_seq_corr{n} = {};
for n=1:numBins
    %m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
            temp = this_data(n,6:18,idx);
%             plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
            temp = this_data(n,19:31,idx);
%             plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
%     shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.seq);
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:)*100,'k');
%     xline(blocks_switch(1)+0.5,'k:');
%     xlim([0,numBlocks_switch])
%     xlabel('Trial block')
%     ylim([0,ceil(nanmax(this_data(:))*100)])
%     ytickformat('percentage')
%     ylabel('Proportion of cells')   
%     SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
%     text(SE(1),SE(2),...
%         ['n = ',num2str(size(temp_data,1)),', \rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
%         'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
%     title([binLabels{n}])
    early_middle_late_seq{n} = temp_data;
    early_middle_late_seq_corr{n} = [size(temp_data,1),this_corr_r,this_corr_p];
end

early_middle_late_ctrl = {};
early_middle_late_ctrl_corr = {};

for n=1:numBins
    %m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
            temp = this_data(n,6:18,idx);
%             plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
            temp = this_data(n,19:31,idx);
%             plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
%     shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.ctrl);
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:)*100,'k');
%     xline(blocks_switch(1)+0.5,'k:');
%     xlim([0,numBlocks_switch])
%     xlabel('Trial block')
%     ylim([0,ceil(nanmax(this_data(:))*100)])
%     ytickformat('percentage')
%     ylabel('Proportion of cells')
%     SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; % NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
%     text(SE(1),SE(2),...
%         ['n = ',num2str(size(temp_data,1)),', \rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
%         'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
%     title([binLabels{n}])
    early_middle_late_ctrl{n} = temp_data;
    early_middle_late_ctrl_corr{n} = [size(temp_data,1),this_corr_r,this_corr_p];
end

for n=1:numBins
    m = m+1; subplot(nrows,ncols,m); hold on;
    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_ctrl{n},1)*100,nansem(early_middle_late_ctrl{n},1)*100,'lineProps',p.col.ctrl);
    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_seq{n},1)*100,nansem(early_middle_late_seq{n},1)*100,'lineProps',p.col.seq);
    xline(blocks_switch(1)+0.5,'k:');
    xlim([0,numBlocks_switch])
    if n==2
        xlabel('Trial block (in blocks of 100 trials)')
    end
    if n==1
        ylim([0,7])
        yticks([0,7])
    else
        ylim([0,5])
        yticks([0,5])
    end
    xticks([1,8,13])
    ylabel('Proportion of cells')
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),...
        ['seq: n = ',num2str(early_middle_late_seq_corr{n}(1)),', \rho = ',num2str(early_middle_late_seq_corr{n}(2),2),', p = ',num2str(early_middle_late_seq_corr{n}(3),2),newline,...
        'ctrl: n = ',num2str(early_middle_late_ctrl_corr{n}(1)),', \rho = ',num2str(early_middle_late_ctrl_corr{n}(2),2),', p = ',num2str(early_middle_late_ctrl_corr{n}(3),2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title([binLabels{n}])
end

%suptitle(['Sequence development (catch trials) across days (',these_animals,')'])
drawnow;
savefig(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Stimulation\new\learning.fig']);
saveas(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Stimulation\new\learning.png']);


%% SUMMARY FIGURE - NEWEST

% select sessions
these_animals = 'stim'; 
these_sessions = ~d_info.excl & d_info.presponsive==1; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1; & d_info.running==1
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

% select data
this_data = warp_exp1_a; %seqCellsByTime_prop;

nrows = 1; ncols = 3; m=0;
F = default_figure([0,3.5,16,3.5]);
m = m+1; subplot(nrows,ncols,m); hold on;

temp_data = [];
for idx=1:d_info.numAnimals
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
        temp = this_data(idx,6:18);
        if showIndividuals
            plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.ctrl]));
        end
        temp_data = [temp_data; temp];
    end
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
        temp = this_data(idx,19:31);
        if showIndividuals
            plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.ctrl]));
        end
        temp_data = [temp_data; temp];
    end
end
temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
if showAverages
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.ctrl);
end
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r_ctrl,this_corr_p_ctrl] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete');  %fitLine(temp(:),temp_data(:),'k');
this_corr_n_ctrl = size(temp_data,1);

temp_data = [];
for idx=1:d_info.numAnimals
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
        temp = this_data(idx,6:18);
        if showIndividuals
            plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.seq]));
        end
        temp_data = [temp_data; temp];
    end
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
        temp = this_data(idx,19:31);
        if showIndividuals
            plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.seq]));
        end
        temp_data = [temp_data; temp];
    end
end
temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
if showAverages
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.seq);
end
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r_seq,this_corr_p_seq] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete'); %fitLine(temp(:),temp_data(:),'k');
this_corr_n_seq = size(temp_data,1);

xline(blocks_switch(1)+0.5,'k:');
xlim([0,numBlocks_switch])
xlabel('Trial block (in blocks of 100 trials)')
ylim([0,0.1])
yticks([0,0.1])
% ylim([0,0.25])
% yticks([0,0.25])
% ylim([0,0.15])
% yticks([0,0.15])
xticks([1,8,13])
ylabel(['Coefficient a',newline,'(initial proportion)'])
if showCorrelation
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),...
        ['seq: n = ',num2str(this_corr_n_seq),', \rho = ',num2str(this_corr_r_seq,2),', p = ',num2str(this_corr_p_seq,2),newline,...
        'ctrl: n = ',num2str(this_corr_n_ctrl),', \rho = ',num2str(this_corr_r_ctrl,2),', p = ',num2str(this_corr_p_ctrl,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
end

% select data
this_data = warp_exp1_b; %seqCellsByTime_prop;
m = m+1; subplot(nrows,ncols,m); hold on;

temp_data = [];
for idx=1:d_info.numAnimals
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
        temp = this_data(idx,6:18);
        if showIndividuals
            plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.ctrl]));
        end
        temp_data = [temp_data; temp];
    end
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
        temp = this_data(idx,19:31);
        if showIndividuals
            plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.ctrl]));
        end
        temp_data = [temp_data; temp];
    end
end
temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
if showAverages
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.ctrl);
end
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r_ctrl,this_corr_p_ctrl] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete');  %fitLine(temp(:),temp_data(:),'k');
this_corr_n_ctrl = size(temp_data,1);

temp_data = [];
yline(0,':');
for idx=1:d_info.numAnimals
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
        temp = this_data(idx,6:18);
        if showIndividuals
            plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.seq]));
        end
        temp_data = [temp_data; temp];
    end
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
        temp = this_data(idx,19:31);
        if showIndividuals
            plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.seq]));
        end
        temp_data = [temp_data; temp];
    end
end
temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
if showAverages
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.seq);
end
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r_seq,this_corr_p_seq] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete'); %fitLine(temp(:),temp_data(:),'k');
this_corr_n_seq = size(temp_data,1);
xline(blocks_switch(1)+0.5,'k:');
xlim([0,numBlocks_switch])
ylim([-0.8,0])
yticks([-0.8,0])
% ylim([-1.8,0.2])
% yticks([-1.8,0.2])
% ylim([-1,0])
% yticks([-1,0])
xticks([1,8,13])
ylabel(['Coefficient b',newline,'(decay rate, s-1)'])
if showCorrelation
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),...
        ['seq: n = ',num2str(this_corr_n_seq),', \rho = ',num2str(this_corr_r_seq,2),', p = ',num2str(this_corr_p_seq,2),newline,...
        'ctrl: n = ',num2str(this_corr_n_ctrl),', \rho = ',num2str(this_corr_r_ctrl,2),', p = ',num2str(this_corr_p_ctrl,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
end

% select data
this_data = seqCellsByTime_prop; %seqCellsByTime_prop;

% select sessions
these_animals = 'stim';
% these_sessions = ~d_info.excl & d_info.presponsive==1; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1;

% preparations
%these_rgbs = discretisedColourMap('winter',true,numBins);
this_n_lower = min(sum(these_sessions(:,2:5),1));
this_n_upper = max(sum(these_sessions(:,2:5),1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

% nrows = 1; ncols = numBins; m=0;
% F = default_figure([-20,0.5,15,5]);

early_middle_late_seq = {};
early_middle_late_seq_corr = {}
early_middle_late_seq_corr{n} = {};
for n=1:numBins
    %m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
            temp = this_data(n,6:18,idx);
%             plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
            temp = this_data(n,19:31,idx);
%             plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
%     shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.seq);
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:)*100,'k');
    early_middle_late_seq{n} = temp_data;
    early_middle_late_seq_corr{n} = [size(temp_data,1),this_corr_r,this_corr_p];
end

early_middle_late_ctrl = {};
early_middle_late_ctrl_corr = {};

for n=1:numBins
    %m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
            temp = this_data(n,6:18,idx);
%             plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
            temp = this_data(n,19:31,idx);
%             plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
%     shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.ctrl);
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:)*100,'k');
    early_middle_late_ctrl{n} = temp_data;
    early_middle_late_ctrl_corr{n} = [size(temp_data,1),this_corr_r,this_corr_p];
end

m = m+1; subplot(nrows,ncols,m); hold on;
xline(blocks_switch(1)+0.5,'k:');
yline(0,'k:');
if showIndividuals
    plot(1:numBlocks_switch,(early_middle_late_ctrl{1}-early_middle_late_ctrl{3})*100,'Color',mean([p.col.gray;p.col.ctrl]));
    plot(1:numBlocks_switch,(early_middle_late_seq{1}-early_middle_late_seq{3})*100,'Color',mean([p.col.gray;p.col.seq]));
end
if showAverages
    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_ctrl{1}-early_middle_late_ctrl{3},1)*100,nansem(early_middle_late_ctrl{1}-early_middle_late_ctrl{3},1)*100,'lineProps',p.col.ctrl);
    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_seq{1}-early_middle_late_seq{3},1)*100,nansem(early_middle_late_seq{1}-early_middle_late_seq{3},1)*100,'lineProps',p.col.seq);
end
ylim([-2,4])
yticks([-2,0,4])
% ylim([-3,8])
% yticks([-3,0,8])
% ylim([-3,5])
% yticks([-3,0,5])
xlim([0,numBlocks_switch])
xlabel('Trial block (in blocks of 100 trials)')
title([binLabels{1},' - ',binLabels{3}])
ylabel('Delta Proportion of cells')

temp_data = early_middle_late_ctrl{1} - early_middle_late_ctrl{3};
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r,this_corr_p] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete');
early_middle_late_ctrl_corr = [size(temp_data,1),this_corr_r,this_corr_p];

temp_data = early_middle_late_seq{1} - early_middle_late_seq{3};
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r,this_corr_p] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete');
early_middle_late_seq_corr = [size(temp_data,1),this_corr_r,this_corr_p];

if showCorrelation
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),...
        ['seq: n = ',num2str(early_middle_late_seq_corr(1)),', \rho = ',num2str(early_middle_late_seq_corr(2),2),', p = ',num2str(early_middle_late_seq_corr(3),2),newline,...
        'ctrl: n = ',num2str(early_middle_late_ctrl_corr(1)),', \rho = ',num2str(early_middle_late_ctrl_corr(2),2),', p = ',num2str(early_middle_late_ctrl_corr(3),2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
end
        
%suptitle(['Sequence development (catch trials) across days (',these_animals,')'])
drawnow;
savefig(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Stimulation\new\learning.fig']);
saveas(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Stimulation\new\learning.png']);


%% LinIdx

% select sessions
these_animals = 'stim'; 
these_sessions_raw = (~d_info.excl & d_info.presponsive==1); % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1; & d_info.running==1

these_sessions = these_sessions_raw;
%these_sessions(~all(these_sessions_raw(:,2:5),2),:) = false;

these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

% select data
this_data = warp_other; %warp_nonlinIdx; %seqCellsByTime_prop; warp_poly1_rsquare; warp_other
this_title = 'other';

nrows = 4; ncols = 3;
F = default_figure([0,0,16,9.9]);

showIndividuals = false;

for m=1:12
    subplot(nrows,ncols,m); hold on;

    temp_data = [];
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
            temp = this_data(idx,6:18,m);
            if showIndividuals
                plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.ctrl]));
            end
            temp_data = [temp_data; temp];
        end
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
            temp = this_data(idx,19:31,m);
            if showIndividuals
                plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.ctrl]));
            end
            temp_data = [temp_data; temp];
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
    if showAverages
        shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.ctrl);
    end
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r_ctrl,this_corr_p_ctrl] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete');  %fitLine(temp(:),temp_data(:),'k');
    this_corr_n_ctrl = size(temp_data,1);

    temp_data = [];
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
            temp = this_data(idx,6:18,m);
            if showIndividuals
                plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.seq]));
            end
            temp_data = [temp_data; temp];
        end
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
            temp = this_data(idx,19:31,m);
            if showIndividuals
                plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.seq]));
            end
            temp_data = [temp_data; temp];
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
    if showAverages
        shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.seq);
    end
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r_seq,this_corr_p_seq] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete'); %fitLine(temp(:),temp_data(:),'k');
    this_corr_n_seq = size(temp_data,1);

    xline(blocks_switch(1)+0.5,'k:');
    xlim([0,numBlocks_switch])
    xlabel('Trial block (in blocks of 100 trials)')

    xticks([1,8,13])
    ylabel([this_title,' ',num2str(m)])
    if showCorrelation
        SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
        text(SE(1),SE(2),...
            ['seq: n = ',num2str(this_corr_n_seq),', \rho = ',num2str(this_corr_r_seq,2),', p = ',num2str(this_corr_p_seq,2),newline,...
            'ctrl: n = ',num2str(this_corr_n_ctrl),', \rho = ',num2str(this_corr_r_ctrl,2),', p = ',num2str(this_corr_p_ctrl,2)],...
            'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    end
end




%% SUMMARY FIGURE - NEWEST

% select sessions
these_animals = 'stim'; 
these_sessions = ~d_info.excl & d_info.presponsive==1;
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

% select data
this_data = seqCellsByTime_prop;  %numSeqCells; %seqCellsByTime_prop;

nrows = 1; ncols = 1; m=0;
F = default_figure([0,3.5,5,3.5]);
m = m+1; subplot(nrows,ncols,m); hold on;

temp_data = [];
for idx=1:d_info.numAnimals
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
        temp = this_data(idx,6:18);
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
        temp = this_data(idx,19:31);
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
end
temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.ctrl);
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r_ctrl,this_corr_p_ctrl] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete');  %fitLine(temp(:),temp_data(:),'k');
this_corr_n_ctrl = size(temp_data,1);

temp_data = [];
for idx=1:d_info.numAnimals
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
        temp = this_data(idx,6:18);
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
        temp = this_data(idx,19:31);
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
end
temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.seq);
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r_seq,this_corr_p_seq] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete'); %fitLine(temp(:),temp_data(:),'k');
this_corr_n_seq = size(temp_data,1);

xline(blocks_switch(1)+0.5,'k:');
xlim([0,numBlocks_switch])
xlabel('Trial block (in blocks of 100 trials)')
ylim([0,200])
yticks([0:100:200])
xticks([1,8,13])
ylabel(['Number of sequence cells'])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['seq: n = ',num2str(this_corr_n_seq),', \rho = ',num2str(this_corr_r_seq,2),', p = ',num2str(this_corr_p_seq,2),newline,...
    'ctrl: n = ',num2str(this_corr_n_ctrl),', \rho = ',num2str(this_corr_r_ctrl,2),', p = ',num2str(this_corr_p_ctrl,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
savefig(F,'C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Supp\numSeqCells_stim.fig');
saveas(F,'C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Supp\numSeqCells_stim.png');


%% !!!!!!!!!!!!!!! Sequence fit across blocks - seq stim vs ctrl stim

nrows = 1; ncols = numBlocks_switch; m=0;
F = default_figure([0,3.5,18,2.5]);

these_sessions = ((~d_info.excl & d_info.presponsive==1))% & d_info.running==0);

x = 0:0.01:5.3;
for m=1:numBlocks_switch
    subplot(nrows,ncols,m); hold on;

    temp_data = [];
    temp_a_ctrl = [];
    temp_b_ctrl = [];
    for idx=1:d_info.numAnimals % [1:14,16:d_info.numAnimals]%
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
            these_blocks = 6:18;
            this_a = warp_exp1_a(idx,these_blocks(m));
            this_b = warp_exp1_b(idx,these_blocks(m));
            this_ypred = this_a*exp(this_b*x);
%             this_a = warp_power2_a(idx,these_blocks(m));
%             this_b = warp_power2_b(idx,these_blocks(m));
%             this_c = warp_power2_c(idx,these_blocks(m));
%             this_ypred = this_a*x.^this_b+this_c;
            temp_data = [temp_data; this_ypred];
            temp_a_ctrl = [temp_a_ctrl; this_a];
            temp_b_ctrl = [temp_b_ctrl; this_b];
            plot(x,this_ypred,'Color',mean([mean([p.col.ctrl;p.col.white]);p.col.white]));
        end
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
            these_blocks = 19:31;
            this_a = warp_exp1_a(idx,these_blocks(m));
            this_b = warp_exp1_b(idx,these_blocks(m));
            this_ypred = this_a*exp(this_b*x);
%             this_a = warp_power2_a(idx,these_blocks(m));
%             this_b = warp_power2_b(idx,these_blocks(m));
%             this_c = warp_power2_c(idx,these_blocks(m));
%             this_ypred = this_a*x.^this_b+this_c;
            temp_data = [temp_data; this_ypred];
            temp_a_ctrl = [temp_a_ctrl; this_a];
            temp_b_ctrl = [temp_b_ctrl; this_b];
            plot(x,this_ypred,'Color',mean([mean([p.col.ctrl;p.col.white]);p.col.white]));
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
    shadedErrorBar(x,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.ctrl);
    disp(['ctrl, this block n=',num2str(size(temp_data,1))])
    
    temp_data = [];
    temp_a_seq = [];
    temp_b_seq = [];
    for idx=1:d_info.numAnimals % [1:14,16:d_info.numAnimals]
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
            these_blocks = 6:18;
            this_a = warp_exp1_a(idx,these_blocks(m));
            this_b = warp_exp1_b(idx,these_blocks(m));
            this_ypred = this_a*exp(this_b*x);
%             this_a = warp_power2_a(idx,these_blocks(m));
%             this_b = warp_power2_b(idx,these_blocks(m));
%             this_c = warp_power2_c(idx,these_blocks(m));
%             this_ypred = this_a*x.^this_b+this_c;
            temp_data = [temp_data; this_ypred];
            temp_a_seq = [temp_a_seq; this_a];
            temp_b_seq = [temp_b_seq; this_b];
            plot(x,this_ypred,'Color',mean([mean([p.col.seq;p.col.white]);p.col.white]));
        end
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
            these_blocks = 19:31;
            this_a = warp_exp1_a(idx,these_blocks(m));
            this_b = warp_exp1_b(idx,these_blocks(m));
            this_ypred = this_a*exp(this_b*x);
%             this_a = warp_power2_a(idx,these_blocks(m));
%             this_b = warp_power2_b(idx,these_blocks(m));
%             this_c = warp_power2_c(idx,these_blocks(m));
%             this_ypred = this_a*x.^this_b+this_c;
            temp_data = [temp_data; this_ypred];
            temp_a_seq = [temp_a_seq; this_a];
            temp_b_seq = [temp_b_seq; this_b];
            plot(x,this_ypred,'Color',mean([mean([p.col.seq;p.col.white]);p.col.white]));
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
    shadedErrorBar(x,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.seq);
    disp(['seq, this block n=',num2str(size(temp_data,1))])

    xlim([0,5.3])
    xticks([0,5.3])
    xticklabels({'0 s','5.3 s'})
    ylim([0,0.1])
    yticks([0,0.1])
    
%     xlim([0,5.3])
%     ylim([0,0.1])
    
    if m==1 || m==numBlocks_switch
        NE = [max(xlim) max(ylim)-max(ylim)/2]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
        text(NE(1),NE(2),...
            ['a = ',num2str(nanmean(temp_a_seq),2),newline,'b = ',num2str(nanmean(temp_b_seq),2),newline,...
            'a = ',num2str(nanmean(temp_a_ctrl),2),newline,'b = ',num2str(nanmean(temp_b_ctrl),2)],...
            'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    end


   % xlabel('Time (s)')
    if m==1
        ylabel('Proportion of seq. cells')
    end
    if m~=1
        set(gca,'yticklabel',{[]})
    end
    title(['Block ',num2str(m)])
end
drawnow;
savefig(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Stimulation\new\development.fig']);
saveas(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Stimulation\new\development.png']);


%% Warp -> performance correlation ---

%F = default_figure([-20,0.5,10,3.5]); hold on;
F = default_figure([0,3.5,8,3.5]); hold on;

this_data_x = warp_exp1_a; %this_data_x = warp_exp1_b;
this_data_y = correct_catch;

% % NEW
% these_sessions = ~d_info.excl & d_info.presponsive==1;
% these_sessions_blocks = [repmat(these_sessions(:,1),1,5),repmat(these_sessions(:,2),1,8),repmat(these_sessions(:,3),1,5),repmat(these_sessions(:,4),1,8),repmat(these_sessions(:,5),1,5)];
% this_data_x(~these_sessions_blocks)=NaN;
% this_data_y(~these_sessions_blocks)=NaN;

this_data_x_seq = [warp_exp1_a(d_info.group==7,6:18); warp_exp1_a(d_info.group==8,19:31)];
this_data_y_seq = [correct_catch(d_info.group==7,6:18); correct_catch(d_info.group==8,19:31)];
this_data_x_ctrl = [warp_exp1_a(d_info.group==8,6:18); warp_exp1_a(d_info.group==7,19:31)];
this_data_y_ctrl = [correct_catch(d_info.group==8,6:18); correct_catch(d_info.group==7,19:31)];

subplot(1,2,1); hold on;
yline(50,'k:');
% plot(this_data_x(:),this_data_y(:)*100,'ko')
plot(this_data_x_ctrl(:),this_data_y_ctrl(:)*100,'.','MarkerSize',10,'Color',p.col.ctrl)
plot(this_data_x_seq(:),this_data_y_seq(:)*100,'.','MarkerSize',10,'Color',p.col.seq)
[this_corr_r_ctrl,this_corr_p_ctrl] = fitLine(this_data_x_ctrl(:),this_data_y_ctrl(:)*100,p.col.ctrl);
[this_corr_r_seq,this_corr_p_seq] = fitLine(this_data_x_seq(:),this_data_y_seq(:)*100,p.col.seq);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k');
xlim([0,0.3])
xticks([0,0.3])
ylim([0,100])
yticks([0,50,100])
xlabel(['Coefficient a',newline,'(initial proportion)'])
ylabel('Performance')
ytickformat('percentage')
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['overall: rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2),newline,...
    'seq: rho = ',num2str(this_corr_r_seq,2),', p = ',num2str(this_corr_p_seq,2),newline,...
    'ctrl: rho = ',num2str(this_corr_r_ctrl,2),', p = ',num2str(this_corr_p_ctrl,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);


% %xlim([0,0.25])
% ylim([0,100])
% xlabel('Exponential coefficient a (initial proportion)')
% ylabel('Performance')
% ytickformat('percentage')
% SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
% text(SE(1),SE(2),...
%     ['\rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
%     'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);

this_data_x = warp_exp1_b; %this_data_x = warp_exp1_b;
this_data_y = correct_catch;

% NEW
% these_sessions = ~d_info.excl & d_info.presponsive==1;
% these_sessions_blocks = [repmat(these_sessions(:,1),1,5),repmat(these_sessions(:,2),1,8),repmat(these_sessions(:,3),1,5),repmat(these_sessions(:,4),1,8),repmat(these_sessions(:,5),1,5)];
% this_data_x(~these_sessions_blocks)=NaN;
% this_data_y(~these_sessions_blocks)=NaN;

this_data_x_seq = [warp_exp1_b(d_info.group==7,6:18); warp_exp1_b(d_info.group==8,19:31)];
this_data_y_seq = [correct_catch(d_info.group==7,6:18); correct_catch(d_info.group==8,19:31)];
this_data_x_ctrl = [warp_exp1_b(d_info.group==8,6:18); warp_exp1_b(d_info.group==7,19:31)];
this_data_y_ctrl = [correct_catch(d_info.group==8,6:18); correct_catch(d_info.group==7,19:31)];

subplot(1,2,2); hold on;
yline(50,'k:');
% plot(this_data_x(:),this_data_y(:)*100,'ko')
plot(this_data_x_ctrl(:),this_data_y_ctrl(:)*100,'.','MarkerSize',10,'Color',p.col.ctrl)
plot(this_data_x_seq(:),this_data_y_seq(:)*100,'.','MarkerSize',10,'Color',p.col.seq)
% [this_corr_r,this_corr_p] = fitLine(this_data_x_ctrl(:),this_data_y_ctrl(:)*100,p.col.ctrl);
% [this_corr_r,this_corr_p] = fitLine(this_data_x_seq(:),this_data_y_seq(:)*100,p.col.seq);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k');
%xlim([0,0.25])
% ylim([0,100])
% xlabel('Exponential coefficient b (decay rate)')
% ylabel('Performance')
% ytickformat('percentage')
% SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
% text(SE(1),SE(2),...
%     ['\rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
%     'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
[this_corr_r_ctrl,this_corr_p_ctrl] = fitLine(this_data_x_ctrl(:),this_data_y_ctrl(:)*100,p.col.ctrl);
[this_corr_r_seq,this_corr_p_seq] = fitLine(this_data_x_seq(:),this_data_y_seq(:)*100,p.col.seq);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k');
xlim([-2,0.5])
xticks([-2,0,0.5])
ylim([0,100])
yticks([0,50,100])
xlabel(['Coefficient b',newline,'(decay rate, s-1)'])
ylabel('Performance')
ytickformat('percentage')
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['overall: rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2),newline,...
    'seq: rho = ',num2str(this_corr_r_seq,2),', p = ',num2str(this_corr_p_seq,2),newline,...
    'ctrl: rho = ',num2str(this_corr_r_ctrl,2),', p = ',num2str(this_corr_p_ctrl,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);

drawnow;
savefig(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Stimulation\new\corr.fig']);
saveas(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Stimulation\new\corr.png']);

disp('Looks like non-responders not yet removed for this plot')




