%% SfN

% d_info = selectDataset(d_info,'-g2-d12345-e01-r01-p01-l01-i00',sheet,path,ops);


%% Preparations

blocks = [5,8,5,8,5];
numBlocks = sum(blocks);

binEdges = 0:5.3/3:5.3
%binEdges = [0,1.5,3,5.3];
numBins = length(binEdges)-1;
binLabels = {'Early peak (0 s - 1.5 s)','Middle peak (1.5 s - 3 s)','Late peak (3 s - 5.3 s)'};

showIndividuals = true;
showAverages = false;
showCorrelation = false;


%% Get performance metrics by block and animal

correct = nan(d_info.numAnimals,numBlocks); % [animal, block]
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
                hitRate(i,this_block) = nanmean(d{i,j}.perf.blocks_general.H(these_blocks_20t));
                crRate(i,this_block) = nanmean(d{i,j}.perf.blocks_general.CR(these_blocks_20t));
            end
        end
    end
end
correct_avg = nanmean(correct,1);
hitRate_avg = nanmean(hitRate,1);
crRate_avg = nanmean(crRate,1);

nonrunners = ~d_info.excl & d_info.running==0;
runners = ~d_info.excl & d_info.running==1;
nonrunners_switch_d1 = false(d_info.numAnimals,1);
runners_switch_d1 = false(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~d_info.excl(i,2) && ~d_info.excl(i,4)
        runners_switch_d1(i) = runners(i,2) && runners(i,4);
        nonrunners_switch_d1(i) = nonrunners(i,2) && nonrunners(i,4);
    elseif ~d_info.excl(i,2)
        runners_switch_d1(i) = runners(i,2);
        nonrunners_switch_d1(i) = nonrunners(i,2);
    elseif ~d_info.excl(i,4)
        runners_switch_d1(i) = runners(i,4);
        nonrunners_switch_d1(i) = nonrunners(i,4);
    end
end
nonrunners_switch_d2 = false(d_info.numAnimals,1);
runners_switch_d2 = false(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~d_info.excl(i,3) && ~d_info.excl(i,5)
        runners_switch_d2(i) = runners(i,3) && runners(i,5);
        nonrunners_switch_d2(i) = nonrunners(i,3) && nonrunners(i,5);
    elseif ~d_info.excl(i,3)
        runners_switch_d2(i) = runners(i,3);
        nonrunners_switch_d2(i) = nonrunners(i,3);
    elseif ~d_info.excl(i,5)
        runners_switch_d2(i) = runners(i,5);
        nonrunners_switch_d2(i) = nonrunners(i,5);
    end
end
nonrunners_switch_d12 = false(d_info.numAnimals,1);
runners_switch_d12 = false(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if (~d_info.excl(i,2) || ~d_info.excl(i,4)) && (~d_info.excl(i,3) || ~d_info.excl(i,5))
        runners_switch_d12(i) = runners_switch_d1(i) && runners_switch_d2(i);
        nonrunners_switch_d12(i) = nonrunners_switch_d1(i) && nonrunners_switch_d2(i);       
    elseif (~d_info.excl(i,2) || ~d_info.excl(i,4))
        runners_switch_d12(i) = runners_switch_d1(i);
        nonrunners_switch_d12(i) = nonrunners_switch_d1(i);
    elseif (~d_info.excl(i,3) || ~d_info.excl(i,5))
        runners_switch_d12(i) = runners_switch_d2(i);
        nonrunners_switch_d12(i) = nonrunners_switch_d2(i);
    end
end


%% Get sequence cells by peak time bin, block and animal

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
        try
            for k=1:blocks(j)
                this_block = temp0+k;

                % get idcs
                this_numCells = length(find(d{i,j}.tng_100t{k}.prop.iscell==1));
                these_idcs_A = find(d{i,j}.tng_100t{k}.passed.AW.Aonly==1);
                these_idcs_X = find(d{i,j}.tng_100t{k}.passed.AW.Xonly==1);

                % get peak locations
                these_peakTimes_A = d{i,j}.tng_100t{k}.firingField.A_AW.peakLocation_s(these_idcs_A);
                these_peakTimes_X = d{i,j}.tng_100t{k}.firingField.X_AW.peakLocation_s(these_idcs_X);
                these_peakTimes = [these_peakTimes_A; these_peakTimes_X];

                % bin peak locations
                temp = discretize(these_peakTimes,binEdges);
                for n=1:numBins
                    seqCellsByTime_num(n,this_block,i) = nansum(temp==n);
                    seqCellsByTime_prop(n,this_block,i) = nansum(temp==n) / this_numCells;
                end
            end
        catch
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
% warp_exp1const_a = nan(d_info.numAnimals,numBlocks); % [animal, block]
% warp_exp1const_b = nan(d_info.numAnimals,numBlocks); % [animal, block]
% warp_exp1const_c = nan(d_info.numAnimals,numBlocks); % [animal, block]
numSeqCells = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if j==1
            temp0=0;
        else
            temp0 = cumsum(blocks);
            temp0 = temp0(j-1);
        end
        if isfield(d{i,j},'warp_100t') %~isempty(d{i,j})
            for k=1:length(d{i,j}.warp_100t)
                this_block = temp0+k;
                
                warp_exp1_a(i,this_block) = d{i,j}.warp_100t{k}.peak.exp1.a;
                warp_exp1_b(i,this_block) = d{i,j}.warp_100t{k}.peak.exp1.b;
                warp_power2_a(i,this_block) = d{i,j}.warp_100t{k}.peak.power2.a;
                warp_power2_b(i,this_block) = d{i,j}.warp_100t{k}.peak.power2.b;
                warp_power2_c(i,this_block) = d{i,j}.warp_100t{k}.peak.power2.c;
%                 warp_exp1const_a(i,this_block) = d{i,j}.warp_100t{k}.peak.exp1const.a;
%                 warp_exp1const_b(i,this_block) = d{i,j}.warp_100t{k}.peak.exp1const.b;
%                 warp_exp1const_c(i,this_block) = d{i,j}.warp_100t{k}.peak.exp1const.c;
                try
                    numSeqCells(i,this_block) = d{i,j}.warp_100t{k}.peak.input.numSequenceCells;
                catch
                end
            end
        end
    end
end

% for each 100t block within switch
blocks_switch = [8,5];
numBlocks_switch = sum(blocks_switch);
warp_exp1_a_switch = nan(d_info.numAnimals,numBlocks_switch); % [animal, block]
warp_exp1_b_switch = nan(d_info.numAnimals,numBlocks_switch); % [animal, block]
warp_power2_a_switch = nan(d_info.numAnimals,numBlocks_switch); % [animal, block]
warp_power2_b_switch = nan(d_info.numAnimals,numBlocks_switch); % [animal, block]
warp_power2_c_switch = nan(d_info.numAnimals,numBlocks_switch); % [animal, block]
% warp_exp1const_a_switch = nan(d_info.numAnimals,numBlocks_switch); % [animal, block]
% warp_exp1const_b_switch = nan(d_info.numAnimals,numBlocks_switch); % [animal, block]
% warp_exp1const_c_switch = nan(d_info.numAnimals,numBlocks_switch); % [animal, block]
numSeqCells_switch = nan(d_info.numAnimals,numBlocks_switch); % [animal, block]
for i=1:d_info.numAnimals
    for k=1:numBlocks_switch
        warp_exp1_a_switch(i,k) = nanmean([warp_exp1_a(i,blocks(1)+k),warp_exp1_a(i,sum(blocks(1:3))+k)],2);
        warp_exp1_b_switch(i,k) = nanmean([warp_exp1_b(i,blocks(1)+k),warp_exp1_b(i,sum(blocks(1:3))+k)],2);
        warp_power2_a_switch(i,k) = nanmean([warp_power2_a(i,blocks(1)+k),warp_power2_a(i,sum(blocks(1:3))+k)],2);
        warp_power2_b_switch(i,k) = nanmean([warp_power2_b(i,blocks(1)+k),warp_power2_b(i,sum(blocks(1:3))+k)],2);
        warp_power2_c_switch(i,k) = nanmean([warp_power2_c(i,blocks(1)+k),warp_power2_c(i,sum(blocks(1:3))+k)],2);
%         warp_exp1const_a_switch(i,k) = nanmean([warp_exp1const_a(i,blocks(1)+k),warp_exp1const_a(i,sum(blocks(1:3))+k)],2);
%         warp_exp1const_b_switch(i,k) = nanmean([warp_exp1const_b(i,blocks(1)+k),warp_exp1const_b(i,sum(blocks(1:3))+k)],2);
%         warp_exp1const_c_switch(i,k) = nanmean([warp_exp1const_c(i,blocks(1)+k),warp_exp1const_c(i,sum(blocks(1:3))+k)],2);
        numSeqCells_switch(i,k) = nanmean([numSeqCells(i,blocks(1)+k),numSeqCells(i,sum(blocks(1:3))+k)],2);
    end
end



%% 5 panels - old version

% select data
this_data_switch = warp_exp1_a_switch;
these_sessions = ~d_info.excl; %& d_info.running==0;
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

nrows = 1; ncols = 5; m=0;
F = default_figure([0,3.5,18,2.5]);

m = m+1; subplot(nrows,ncols,m); hold on;
temp_data = [];
for idx=1:d_info.numAnimals
    if nonrunners_switch_d12(idx)
        temp = nan(1,size(these_blocks_switch,2));
        temp(1,find(these_blocks_switch(idx,:))) = this_data_switch(idx,find(these_blocks_switch(idx,:)));
        if showIndividuals
            plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.nonrunner])); 
        end
        temp_data = [temp_data; temp];
    end
end
if showAverages
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.nonrunner);
end
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r,this_corr_p] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete'); %fitLine(temp(:),temp_data(:),'k');
exponential_nonrunner_corr = [size(temp_data,1),this_corr_r,this_corr_p];
these_sessions_running = ~d_info.excl & d_info.running==1; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1; & d_info.running==1;
temp_data = [];
for idx=1:d_info.numAnimals
    if runners_switch_d12(idx)
        temp = nan(1,size(these_blocks_switch,2));
        temp(1,find(these_blocks_switch(idx,:))) = this_data_switch(idx,find(these_blocks_switch(idx,:)));
        if showIndividuals
            plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.runner]));
        end
        temp_data = [temp_data; temp];
    end
end
if showAverages
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.runner);
end
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r,this_corr_p] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete'); %fitLine(temp(:),temp_data(:),'k');
exponential_runner_corr = [size(temp_data,1),this_corr_r,this_corr_p];
xline(blocks_switch(1)+0.5,'k:');
xlim([0,numBlocks_switch])
ylim([0,0.25])
xlabel('Trial block (in blocks of 100 trials)')
ylim([0,0.25])
yticks([0,0.25])
xticks([1,8,13])
if showCorrelation
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),...
        ['n = ',num2str(exponential_nonrunner_corr(1)),', rho = ',num2str(exponential_nonrunner_corr(2),2),', p = ',num2str(exponential_nonrunner_corr(3),2),newline,...
        'n = ',num2str(exponential_runner_corr(1)),', rho = ',num2str(exponential_runner_corr(2),2),', p = ',num2str(exponential_runner_corr(3),2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
end
ylabel(['Coefficient a',newline,'(initial proportion)'])

this_data_switch = warp_exp1_b_switch;
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

m = m+1; subplot(nrows,ncols,m); hold on;
yline(0,'k:')
these_sessions_nonrunning = ~d_info.excl & d_info.running==0; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1; & d_info.running==1;
temp_data = [];
for idx=1:d_info.numAnimals
    if nonrunners_switch_d12(idx)
        temp = nan(1,size(these_blocks_switch,2));
        temp(1,find(these_blocks_switch(idx,:))) = this_data_switch(idx,find(these_blocks_switch(idx,:)));
        if showIndividuals
            plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.nonrunner]));
        end
        temp_data = [temp_data; temp];
    end
end
if showAverages
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.nonrunner);
end
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r,this_corr_p] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete'); %fitLine(temp(:),temp_data(:),'k');
exponential_nonrunner_corr = [size(temp_data,1),this_corr_r,this_corr_p];
these_sessions_running = ~d_info.excl & d_info.running==1; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1; & d_info.running==1;
temp_data = [];
for idx=1:d_info.numAnimals
    if runners_switch_d12(idx)
        temp = nan(1,size(these_blocks_switch,2));
        temp(1,find(these_blocks_switch(idx,:))) = this_data_switch(idx,find(these_blocks_switch(idx,:)));
        if showIndividuals
            plot(1:numBlocks_switch,temp,'Color',mean([p.col.gray;p.col.runner]));
        end
        temp_data = [temp_data; temp];
    end
end
shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.runner);
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r,this_corr_p] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete'); %fitLine(temp(:),temp_data(:),'k');
exponential_runner_corr = [size(temp_data,1),this_corr_r,this_corr_p];
xline(blocks_switch(1)+0.5,'k:');
xlim([0,numBlocks_switch])
ylim([-2,0.5]) %ylim([-1.5,0.5])
yticks([-2,0.5]) % yticks([-1.5,0.5])
xticks([1,8,13])
if showCorrelation
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),...
        ['n = ',num2str(exponential_nonrunner_corr(1)),', rho = ',num2str(exponential_nonrunner_corr(2),2),', p = ',num2str(exponential_nonrunner_corr(3),2),newline,...
        'n = ',num2str(exponential_runner_corr(1)),', rho = ',num2str(exponential_runner_corr(2),2),', p = ',num2str(exponential_runner_corr(3),2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
end
ylabel(['Coefficient b',newline,'(decay rate, s-1)'])

% select data
%this_data = seqCellsByTime_prop; %seqCellsByTime_prop;
this_data_switch = seqCellsByTime_switch_prop; %seqCellsByTime_switch_prop;

% select sessions
these_sessions = ~d_info.excl;% & d_info.running==1;
%these_animals = 'non-runners'; these_sessions = ~d_info.excl & d_info.running==0;

% preparations
%these_rgbs = discretisedColourMap('winter',true,numBins);
this_n_lower = min(sum(these_sessions,1));
this_n_upper = max(sum(these_sessions,1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

early_middle_late_runner = {};
early_middle_late_runner_corr = {};
for n=1:numBins
%    m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if runners_switch_d12(idx)
            temp = nan(1,numBlocks_switch);
            temp(find(these_blocks_switch(idx,:))) = this_data_switch(n,find(these_blocks_switch(idx,:)),idx);
%             if showIndividuals
%                 plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
%             end
            temp_data = [temp_data; temp];
        end
    end
    %shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.runner);
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:)*100,'k');
%     xline(blocks_switch(1)+0.5,'k:');
%     xlim([0,numBlocks_switch])
%     xlabel('Trial block')
%     ylim([0,15])
%     ytickformat('percentage')
%     ylabel('Proportion of cells')
%     SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
%     text(SE(1),SE(2),...
%         ['n = ',num2str(size(temp_data,1)),', rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
%         'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
%     title(binLabels{n})
    early_middle_late_runner_corr{n} = [size(temp_data,1),this_corr_r,this_corr_p];
    early_middle_late_runner{n} = temp_data;
end


% select sessions
%these_animals = 'runners'; these_sessions = ~d_info.excl & d_info.running==1;
%these_animals = 'non-runners'; these_sessions = ~d_info.excl & d_info.running==0;

% preparations
%these_rgbs = discretisedColourMap('winter',true,numBins);
this_n_lower = min(sum(these_sessions,1));
this_n_upper = max(sum(these_sessions,1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

early_middle_late_nonrunner = {};
early_middle_late_nonrunner_corr = {};
for n=1:numBins
%    m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if nonrunners_switch_d12(idx)
            temp = nan(1,numBlocks_switch);
            temp(find(these_blocks_switch(idx,:))) = this_data_switch(n,find(these_blocks_switch(idx,:)),idx);
%             if showIndividuals
%                 plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
%             end
            temp_data = [temp_data; temp];
        end
    end
    %shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.nonrunner);
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:)*100,'k');
%     xline(blocks_switch(1)+0.5,'k:');
%     xlim([0,numBlocks_switch])
%     xlabel('Trial block')
%     ylim([0,15])
%     ytickformat('percentage')
%     ylabel('Proportion of cells')
%     SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
%     text(SE(1),SE(2),...
%         ['n = ',num2str(size(temp_data,1)),', rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
%         'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
%     title(binLabels{n})
    early_middle_late_nonrunner_corr{n} = [size(temp_data,1),this_corr_r,this_corr_p];
    early_middle_late_nonrunner{n} = temp_data;
end

for n=1:numBins
    m = m+1; subplot(nrows,ncols,m); hold on;
    if showIndividuals
        plot(1:numBlocks_switch,early_middle_late_nonrunner{n}*100,'Color',mean([p.col.gray;p.col.nonrunner]));
        plot(1:numBlocks_switch,early_middle_late_runner{n}*100,'Color',mean([p.col.gray;p.col.runner]));
    end
    if showAverages
        shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_nonrunner{n},1)*100,nansem(early_middle_late_nonrunner{n},1)*100,'lineProps',p.col.nonrunner);
        shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_runner{n},1)*100,nansem(early_middle_late_runner{n},1)*100,'lineProps',p.col.runner);
    end
    xline(blocks_switch(1)+0.5,'k:');
    xlim([0,numBlocks_switch])
    if n==2
        xlabel('Trial block (in blocks of 100 trials)')
    end
    if n==1
        ylim([0,25]) %ylim([0,15])
        yticks([0,25]) %yticks([0,15])
    else
        ylim([0,15]) %ylim([0,10])
        yticks([0,15]) %yticks([0,10])
    end
    xticks([1,8,13])
    ytickformat('percentage')
    ylabel('Proportion of cells')
    if showCorrelation
        SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
        text(SE(1),SE(2),...
            ['n = ',num2str(early_middle_late_nonrunner_corr{n}(1)),', rho = ',num2str(early_middle_late_nonrunner_corr{n}(2),2),', p = ',num2str(early_middle_late_nonrunner_corr{n}(3),2),newline,...
            'n = ',num2str(early_middle_late_runner_corr{n}(1)),', rho = ',num2str(early_middle_late_runner_corr{n}(2),2),', p = ',num2str(early_middle_late_runner_corr{n}(3),2)],...
            'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    end
    title(binLabels{n})
end

%suptitle(['Sequence development across days'])
drawnow;
savefig(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Imaging\new\learning.fig']);
saveas(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Imaging\new\learning.png']);




%% !!!!!!!!!!!!!!! Sequence fit across blocks - runners vs. non-runners

nrows = 1; ncols = numBlocks_switch; m=0;
F = default_figure([0,3.5,18,2.5]);

x = 0:0.01:5.3;
for m=1:numBlocks_switch
    subplot(nrows,ncols,m); hold on;

    temp_data = [];
    temp_a_nonrunners = [];
    temp_b_nonrunners = [];
    for idx=1:d_info.numAnimals % [1:14,16:d_info.numAnimals]%
        if nonrunners_switch_d12(idx)
            these_blocks_1 = 6:18; these_blocks_2 = 19:31;
            this_a = nanmean([warp_exp1_a(idx,these_blocks_1(m)),warp_exp1_a(idx,these_blocks_2(m))]);
            this_b = nanmean([warp_exp1_b(idx,these_blocks_1(m)),warp_exp1_b(idx,these_blocks_2(m))]);
            this_ypred = this_a*exp(this_b*x);
            temp_data = [temp_data; this_ypred];
            temp_a_nonrunners = [temp_a_nonrunners; this_a];
            temp_b_nonrunners = [temp_b_nonrunners; this_b];
        end
    end
    %temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
    shadedErrorBar(x,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.nonrunner);
    
    temp_data = [];
    temp_a_runners = [];
    temp_b_runners = [];
    for idx=1:d_info.numAnimals % [1:14,16:d_info.numAnimals]
        if runners_switch_d12(idx)
            these_blocks_1 = 6:18; these_blocks_2 = 19:31;
            this_a = nanmean([warp_exp1_a(idx,these_blocks_1(m)),warp_exp1_a(idx,these_blocks_2(m))]);
            this_b = nanmean([warp_exp1_b(idx,these_blocks_1(m)),warp_exp1_b(idx,these_blocks_2(m))]);
            this_ypred = this_a*exp(this_b*x);
            temp_data = [temp_data; this_ypred];
            temp_a_runners = [temp_a_runners; this_a];
            temp_b_runners = [temp_b_runners; this_b];
        end
    end
    %temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
    shadedErrorBar(x,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.runner);
    xlim([0,5.3])
    xticks([0,5.3])
    xticklabels({'0 s','5.3 s'})
    ylim([0,0.1])
    yticks([0,0.1])
    
    if m==1 || m==numBlocks_switch
        NE = [max(xlim) max(ylim)-max(ylim)/2]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
        text(NE(1),NE(2),...
            ['a = ',num2str(nanmean(temp_a_nonrunners),2),newline,'b = ',num2str(nanmean(temp_b_nonrunners),2),newline,...
            'a = ',num2str(nanmean(temp_a_runners),2),newline,'b = ',num2str(nanmean(temp_b_runners),2)],...
            'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    end

    %xlabel('Time (s)')
    if m==1
        ylabel('Proportion of seq. cells')
    end
    if m~=1
        set(gca,'yticklabel',{[]})
    end
    title(['Block ',num2str(m)])
end
drawnow;
savefig(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Imaging\new\development.fig']);
saveas(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Imaging\new\development.png']);


%% Warp -> performance correlation ---

%F = default_figure([-20,0.5,10,3.5]); hold on;
F = default_figure([0,3.5,8,3.5]); hold on;

these_sessions_runners = ~d_info.excl & d_info.running==1;
these_blocks_runners = [repmat(these_sessions_runners(:,1),1,blocks(1)),repmat(these_sessions_runners(:,2),1,blocks(2)),repmat(these_sessions_runners(:,3),1,blocks(3)),repmat(these_sessions_runners(:,4),1,blocks(4)),repmat(these_sessions_runners(:,5),1,blocks(5))];
these_sessions_nonrunners = ~d_info.excl & d_info.running==0;
these_blocks_nonrunners = [repmat(these_sessions_nonrunners(:,1),1,blocks(1)),repmat(these_sessions_nonrunners(:,2),1,blocks(2)),repmat(these_sessions_nonrunners(:,3),1,blocks(3)),repmat(these_sessions_nonrunners(:,4),1,blocks(4)),repmat(these_sessions_nonrunners(:,5),1,blocks(5))];

this_data_x = warp_exp1_a;
this_data_y = correct;
this_data_x_runners = [warp_exp1_a(these_blocks_runners); warp_exp1_a(these_blocks_runners)];
this_data_y_runners = [correct(these_blocks_runners); correct(these_blocks_runners)];
this_data_x_nonrunners = [warp_exp1_a(these_blocks_nonrunners); warp_exp1_a(these_blocks_nonrunners)];
this_data_y_nonrunners = [correct(these_blocks_nonrunners); correct(these_blocks_nonrunners)];

subplot(1,2,1); hold on;
yline(50,'k:');
% plot(this_data_x(:),this_data_y(:)*100,'ko')
plot(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'.','MarkerSize',10,'Color',p.col.nonrunner)
plot(this_data_x_runners(:),this_data_y_runners(:)*100,'.','MarkerSize',10,'Color',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k');
xlim([0,0.35])
xticks([0,0.35])
ylim([0,100])
yticks([0,50,100])
xlabel(['Coefficient a',newline,'(initial proportion)'])
ylabel('Performance')
ytickformat('percentage')
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['overall: rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2),newline,...
    'non-runners: rho = ',num2str(this_corr_r_nonrunners,2),', p = ',num2str(this_corr_p_nonrunners,2),newline,...
    'runners: rho = ',num2str(this_corr_r_runners,2),', p = ',num2str(this_corr_p_runners,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);

this_data_x = warp_exp1_b;
this_data_y = correct;
this_data_x_runners = [warp_exp1_b(these_blocks_runners); warp_exp1_b(these_blocks_runners)];
this_data_y_runners = [correct(these_blocks_runners); correct(these_blocks_runners)];
this_data_x_nonrunners = [warp_exp1_b(these_blocks_nonrunners); warp_exp1_b(these_blocks_nonrunners)];
this_data_y_nonrunners = [correct(these_blocks_nonrunners); correct(these_blocks_nonrunners)];

subplot(1,2,2); hold on;
xline(0,'k:');
yline(50,'k:');
% plot(this_data_x(:),this_data_y(:)*100,'ko')
plot(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'.','MarkerSize',10,'Color',p.col.nonrunner)
plot(this_data_x_runners(:),this_data_y_runners(:)*100,'.','MarkerSize',10,'Color',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k');
xlim([-2.5,0.5])
xticks([-2.5,0,0.5])
ylim([0,100])
yticks([0,50,100])
xlabel(['Coefficient b',newline,'(decay rate, s-1)'])
ylabel('Performance')
ytickformat('percentage')
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['overall: rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2),newline,...
    'non-runners: rho = ',num2str(this_corr_r_nonrunners,2),', p = ',num2str(this_corr_p_nonrunners,2),newline,...
    'runners: rho = ',num2str(this_corr_r_runners,2),', p = ',num2str(this_corr_p_runners,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);

drawnow;
savefig(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Imaging\new\corr.fig']);
saveas(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Imaging\new\corr.png']);




%% 5 (-> 3!) panels - new version

% select data
this_data_switch = warp_exp1_a_switch;
these_sessions = ~d_info.excl; %& d_info.running==0;
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

nrows = 1; ncols = 3; m=0;
F = default_figure([0,3.5,16,3.5]);

m = m+1; subplot(nrows,ncols,m); hold on;
temp_data = [];
for idx=1:d_info.numAnimals
    if nonrunners_switch_d12(idx)
        temp = nan(1,size(these_blocks_switch,2));
        temp(1,find(these_blocks_switch(idx,:))) = this_data_switch(idx,find(these_blocks_switch(idx,:)));
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
end
shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.nonrunner);
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r,this_corr_p] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete'); %fitLine(temp(:),temp_data(:),'k');
exponential_nonrunner_corr = [size(temp_data,1),this_corr_r,this_corr_p];
these_sessions_running = ~d_info.excl & d_info.running==1; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1; & d_info.running==1;
temp_data = [];
for idx=1:d_info.numAnimals
    if runners_switch_d12(idx)
        temp = nan(1,size(these_blocks_switch,2));
        temp(1,find(these_blocks_switch(idx,:))) = this_data_switch(idx,find(these_blocks_switch(idx,:)));
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
end
shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.runner);
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r,this_corr_p] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete'); %fitLine(temp(:),temp_data(:),'k');
exponential_runner_corr = [size(temp_data,1),this_corr_r,this_corr_p];
xline(blocks_switch(1)+0.5,'k:');
xlim([0,numBlocks_switch])
ylim([0,0.25])
xlabel('Trial block (in blocks of 100 trials)')
ylim([0,0.25])
yticks([0,0.25])
xticks([1,8,13])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['n = ',num2str(exponential_nonrunner_corr(1)),', rho = ',num2str(exponential_nonrunner_corr(2),2),', p = ',num2str(exponential_nonrunner_corr(3),2),newline,...
    'n = ',num2str(exponential_runner_corr(1)),', rho = ',num2str(exponential_runner_corr(2),2),', p = ',num2str(exponential_runner_corr(3),2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
% text(SE(1),SE(2),...
%     ['non-runner: n = ',num2str(exponential_nonrunner_corr(1)),', rho = ',num2str(exponential_nonrunner_corr(2),2),', p = ',num2str(exponential_nonrunner_corr(3),2),newline,...
%     'runner: n = ',num2str(exponential_runner_corr(1)),', rho = ',num2str(exponential_runner_corr(2),2),', p = ',num2str(exponential_runner_corr(3),2)],...
%     'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
ylabel(['Coefficient a',newline,'(initial proportion)'])

this_data_switch = warp_exp1_b_switch;
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

m = m+1; subplot(nrows,ncols,m); hold on;
yline(0,'k:')
these_sessions_nonrunning = ~d_info.excl & d_info.running==0; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1; & d_info.running==1;
temp_data = [];
for idx=1:d_info.numAnimals
    if nonrunners_switch_d12(idx)
        temp = nan(1,size(these_blocks_switch,2));
        temp(1,find(these_blocks_switch(idx,:))) = this_data_switch(idx,find(these_blocks_switch(idx,:)));
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
end
shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.nonrunner);
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r,this_corr_p] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete'); %fitLine(temp(:),temp_data(:),'k');
exponential_nonrunner_corr = [size(temp_data,1),this_corr_r,this_corr_p];
these_sessions_running = ~d_info.excl & d_info.running==1; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1; & d_info.running==1;
temp_data = [];
for idx=1:d_info.numAnimals
    if runners_switch_d12(idx)
        temp = nan(1,size(these_blocks_switch,2));
        temp(1,find(these_blocks_switch(idx,:))) = this_data_switch(idx,find(these_blocks_switch(idx,:)));
        %plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
end
shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.runner);
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r,this_corr_p] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete'); %fitLine(temp(:),temp_data(:),'k');
exponential_runner_corr = [size(temp_data,1),this_corr_r,this_corr_p];
xline(blocks_switch(1)+0.5,'k:');
xlim([0,numBlocks_switch])
ylim([-1.5,0.5])
yticks([-1.5,0.5])
xticks([1,8,13])
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['n = ',num2str(exponential_nonrunner_corr(1)),', rho = ',num2str(exponential_nonrunner_corr(2),2),', p = ',num2str(exponential_nonrunner_corr(3),2),newline,...
    'n = ',num2str(exponential_runner_corr(1)),', rho = ',num2str(exponential_runner_corr(2),2),', p = ',num2str(exponential_runner_corr(3),2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
% text(SE(1),SE(2),...
%     ['non-runner: n = ',num2str(exponential_nonrunner_corr(1)),', rho = ',num2str(exponential_nonrunner_corr(2),2),', p = ',num2str(exponential_nonrunner_corr(3),2),newline,...
%     'runner: n = ',num2str(exponential_runner_corr(1)),', rho = ',num2str(exponential_runner_corr(2),2),', p = ',num2str(exponential_runner_corr(3),2)],...
%     'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
ylabel(['Coefficient b',newline,'(decay rate, s-1)'])

% select data
%this_data = seqCellsByTime_prop; %seqCellsByTime_prop;
this_data_switch = seqCellsByTime_switch_prop; %seqCellsByTime_switch_prop;

% select sessions
these_sessions = ~d_info.excl;% & d_info.running==1;
%these_animals = 'non-runners'; these_sessions = ~d_info.excl & d_info.running==0;

% preparations
%these_rgbs = discretisedColourMap('winter',true,numBins);
this_n_lower = min(sum(these_sessions,1));
this_n_upper = max(sum(these_sessions,1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

early_middle_late_runner = {};
early_middle_late_runner_corr = {};
for n=1:numBins
%    m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if runners_switch_d12(idx)
            temp = nan(1,numBlocks_switch);
            temp(find(these_blocks_switch(idx,:))) = this_data_switch(n,find(these_blocks_switch(idx,:)),idx);
            %plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    %shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.runner);
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:)*100,'k');
%     xline(blocks_switch(1)+0.5,'k:');
%     xlim([0,numBlocks_switch])
%     xlabel('Trial block')
%     ylim([0,15])
%     ytickformat('percentage')
%     ylabel('Proportion of cells')
%     SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
%     text(SE(1),SE(2),...
%         ['n = ',num2str(size(temp_data,1)),', rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
%         'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
%     title(binLabels{n})
    early_middle_late_runner_corr{n} = [size(temp_data,1),this_corr_r,this_corr_p];
    early_middle_late_runner{n} = temp_data;
end


% select sessions
%these_animals = 'runners'; these_sessions = ~d_info.excl & d_info.running==1;
%these_animals = 'non-runners'; these_sessions = ~d_info.excl & d_info.running==0;

% preparations
%these_rgbs = discretisedColourMap('winter',true,numBins);
this_n_lower = min(sum(these_sessions,1));
this_n_upper = max(sum(these_sessions,1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

early_middle_late_nonrunner = {};
early_middle_late_nonrunner_corr = {};
for n=1:numBins
%    m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if nonrunners_switch_d12(idx)
            temp = nan(1,numBlocks_switch);
            temp(find(these_blocks_switch(idx,:))) = this_data_switch(n,find(these_blocks_switch(idx,:)),idx);
            %plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    %shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.nonrunner);
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:)*100,'k');
%     xline(blocks_switch(1)+0.5,'k:');
%     xlim([0,numBlocks_switch])
%     xlabel('Trial block')
%     ylim([0,15])
%     ytickformat('percentage')
%     ylabel('Proportion of cells')
%     SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
%     text(SE(1),SE(2),...
%         ['n = ',num2str(size(temp_data,1)),', rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
%         'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
%     title(binLabels{n})
    early_middle_late_nonrunner_corr{n} = [size(temp_data,1),this_corr_r,this_corr_p];
    early_middle_late_nonrunner{n} = temp_data;
end

m = m+1; subplot(nrows,ncols,m); hold on;
% for n=1:numBins
%     m = m+1; subplot(nrows,ncols,m); hold on;
%     shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_nonrunner{n},1)*100,nansem(early_middle_late_nonrunner{n},1)*100,'lineProps',p.col.nonrunner);
%     shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_runner{n},1)*100,nansem(early_middle_late_runner{n},1)*100,'lineProps',p.col.runner);
%     xline(blocks_switch(1)+0.5,'k:');
%     xlim([0,numBlocks_switch])
%     if n==2
%         xlabel('Trial block (in blocks of 100 trials)')
%     end
%     if n==1
%         ylim([0,15])
%         yticks([0,15])
%     else
%         ylim([0,10])
%         yticks([0,10])
%     end
%     xticks([1,8,13])
%     ytickformat('percentage')
%     ylabel('Proportion of cells')
%     SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
%     %NE = [max(xlim) max(ylim)-max(ylim)/5]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
% %     text(NE(1),NE(2),...
% %         ['non-runners: n = ',num2str(early_middle_late_nonrunner_corr{n}(1)),', rho = ',num2str(early_middle_late_nonrunner_corr{n}(2),2),', p = ',num2str(early_middle_late_nonrunner_corr{n}(3),2),newline,...
% %         'runners: n = ',num2str(early_middle_late_runner_corr{n}(1)),', rho = ',num2str(early_middle_late_runner_corr{n}(2),2),', p = ',num2str(early_middle_late_runner_corr{n}(3),2)],...
% %         'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
%     text(SE(1),SE(2),...
%         ['n = ',num2str(early_middle_late_nonrunner_corr{n}(1)),', rho = ',num2str(early_middle_late_nonrunner_corr{n}(2),2),', p = ',num2str(early_middle_late_nonrunner_corr{n}(3),2),newline,...
%         'n = ',num2str(early_middle_late_runner_corr{n}(1)),', rho = ',num2str(early_middle_late_runner_corr{n}(2),2),', p = ',num2str(early_middle_late_runner_corr{n}(3),2)],...
%         'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
%     title(binLabels{n})
% end


xline(blocks_switch(1)+0.5,'k:');
yline(0,'k:');
shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_nonrunner{1}-early_middle_late_nonrunner{3},1)*100,nansem(early_middle_late_nonrunner{1}-early_middle_late_nonrunner{3},1)*100,'lineProps',p.col.nonrunner);
shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_runner{1}-early_middle_late_runner{3},1)*100,nansem(early_middle_late_runner{1}-early_middle_late_runner{3},1)*100,'lineProps',p.col.runner);
ylim([-6,12])
yticks([-6,0,12])
xlim([0,numBlocks_switch])
xlabel('Trial block (in blocks of 100 trials)')
title([binLabels{1},' - ',binLabels{3}])
ylabel('Delta Proportion of cells')

temp_data = early_middle_late_nonrunner{1} - early_middle_late_nonrunner{3};
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r,this_corr_p] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete');
early_middle_late_nonrunner_corr = [size(temp_data,1),this_corr_r,this_corr_p];

temp_data = early_middle_late_runner{1} - early_middle_late_runner{3};
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r,this_corr_p] = corr(temp(:),temp_data(:),'Type','Pearson','Rows','Complete');
early_middle_late_runner_corr = [size(temp_data,1),this_corr_r,this_corr_p];

SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['non-runner: n = ',num2str(early_middle_late_nonrunner_corr(1)),', \rho = ',num2str(early_middle_late_nonrunner_corr(2),2),', p = ',num2str(early_middle_late_nonrunner_corr(3),2),newline,...
    'runner: n = ',num2str(early_middle_late_runner_corr(1)),', \rho = ',num2str(early_middle_late_runner_corr(2),2),', p = ',num2str(early_middle_late_runner_corr(3),2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);

       

%suptitle(['Sequence development across days'])
drawnow;
savefig(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Imaging\new\learning.fig']);
saveas(F,['C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Warping\Switch Imaging\new\learning.png']);

