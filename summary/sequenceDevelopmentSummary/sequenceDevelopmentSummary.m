%% Preparations

blocks = [5,8,5,8,5];
numBlocks = sum(blocks);

% binEdges = [0,1,3,5.4]; %[0:5.3/3:5.3]; %[0:1:5];
% numBins = length(binEdges)-1;
% %binLabels = {'Early sequence cells (sensory cells)','Middle sequence cells (early time cells)','Late sequence cells (late time cells)'};
% %binLabels = {'Sequence cells with peak time 0 s - 1 s','Sequence cells with peak time 1 s - 2 s','Sequence cells with peak time 2 s - 3 s','Sequence cells with peak time 3 s - 4 s','Sequence cells with peak time 4 s - 5 s'};
% binLabels = {'Sequence cells with peak time 0 s - 1 s','Sequence cells with peak time 1 s - 3 s','Sequence cells with peak time 3 s - 5.3 s'};

binEdges = [0,1.5,3,5.3];
numBins = length(binEdges)-1;
binLabels = {'Early peak (0 s - 1.5 s)','Middle peak (1.5 s - 3 s)','Late peak (3 s - 5.3 s)'};


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
        try %if ~isempty(d{i,j})
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
warp_exp1const_a_switch = nan(d_info.numAnimals,numBlocks_switch); % [animal, block]
warp_exp1const_b_switch = nan(d_info.numAnimals,numBlocks_switch); % [animal, block]
warp_exp1const_c_switch = nan(d_info.numAnimals,numBlocks_switch); % [animal, block]
numSeqCells_switch = nan(d_info.numAnimals,numBlocks_switch); % [animal, block]
for i=1:d_info.numAnimals
    for k=1:numBlocks_switch
        warp_exp1_a_switch(i,k) = nanmean([warp_exp1_a(i,blocks(1)+k),warp_exp1_a(i,sum(blocks(1:3))+k)],2);
        warp_exp1_b_switch(i,k) = nanmean([warp_exp1_b(i,blocks(1)+k),warp_exp1_b(i,sum(blocks(1:3))+k)],2);
        warp_power2_a_switch(i,k) = nanmean([warp_power2_a(i,blocks(1)+k),warp_power2_a(i,sum(blocks(1:3))+k)],2);
        warp_power2_b_switch(i,k) = nanmean([warp_power2_b(i,blocks(1)+k),warp_power2_b(i,sum(blocks(1:3))+k)],2);
        warp_power2_c_switch(i,k) = nanmean([warp_power2_c(i,blocks(1)+k),warp_power2_c(i,sum(blocks(1:3))+k)],2);
        warp_exp1const_a_switch(i,k) = nanmean([warp_exp1const_a(i,blocks(1)+k),warp_exp1const_a(i,sum(blocks(1:3))+k)],2);
        warp_exp1const_b_switch(i,k) = nanmean([warp_exp1const_b(i,blocks(1)+k),warp_exp1const_b(i,sum(blocks(1:3))+k)],2);
        warp_exp1const_c_switch(i,k) = nanmean([warp_exp1const_c(i,blocks(1)+k),warp_exp1const_c(i,sum(blocks(1:3))+k)],2);
        numSeqCells_switch(i,k) = nanmean([numSeqCells(i,blocks(1)+k),numSeqCells(i,sum(blocks(1:3))+k)],2);
    end
end


%% Correlate sequence cells by peak time bin with performance metrics

% corr_aw_d12345.seqCellsByTime_num_correct = nan(d_info.numAnimals,numBins); 
corr_aw_d12345.seqCellsByTime_prop_correct = nan(d_info.numAnimals,numBins);
corr_aw_d1.seqCellsByTime_prop_correct = nan(d_info.numAnimals,numBins);
corr_aw_d23.seqCellsByTime_prop_correct = nan(d_info.numAnimals,numBins);
corr_aw_d45.seqCellsByTime_prop_correct = nan(d_info.numAnimals,numBins);
corr_aw_d23d45.seqCellsByTime_prop_correct = nan(d_info.numAnimals,numBins);
for i=1:d_info.numAnimals
    for n=1:numBins
%         corr_aw_d12345.seqCellsByTime_num_correct(i,n) = ...
%             corr(seqCellsByTime_num(n,:,i)',correct(i,:)','Type','Pearson','Rows','Complete');
        corr_aw_d12345.seqCellsByTime_prop_correct(i,n) = ...
            corr(seqCellsByTime_prop(n,:,i)',correct(i,:)','Type','Pearson','Rows','Complete');
        corr_aw_d1.seqCellsByTime_prop_correct(i,n) = ...
            corr(seqCellsByTime_prop(n,1:5,i)',correct(i,1:5)','Type','Pearson','Rows','Complete');
        corr_aw_d23.seqCellsByTime_prop_correct(i,n) = ...
            corr(seqCellsByTime_prop(n,6:18,i)',correct(i,6:18)','Type','Pearson','Rows','Complete');
        corr_aw_d45.seqCellsByTime_prop_correct(i,n) = ...
            corr(seqCellsByTime_prop(n,19:31,i)',correct(i,19:31)','Type','Pearson','Rows','Complete');
        corr_aw_d23d45.seqCellsByTime_prop_correct(i,n) = nanmean([corr_aw_d23.seqCellsByTime_prop_correct(i,n),corr_aw_d45.seqCellsByTime_prop_correct(i,n)]);
    end
end


%% --- %%% --- %%% --- %%% --- %%% --- %%% ---


%% --- SUMMARY FIGURE --- Proportion of cells as a function of trial block

% select sessions
%these_animals = 'all animals'; these_sessions = ~d_info.excl & d_info.engagement==1;
%these_animals = 'disengaged animals'; these_sessions = ~d_info.excl & d_info.engagement==0;
%these_animals = 'runners'; these_sessions = ~d_info.excl & d_info.engagement==1 & d_info.running==1;
%these_animals = 'non-runners'; these_sessions = ~d_info.excl & d_info.engagement==1 & d_info.running==0;
%these_animals = 'learners'; these_sessions = ~d_info.excl & d_info.engagement==1 & d_info.learning==1;
%these_animals = 'non-learners'; these_sessions = ~d_info.excl & d_info.engagement==1 & d_info.learning==0;

%these_animals = 'runners'; these_sessions = ~d_info.excl & d_info.running==1;
these_animals = 'non-runners'; these_sessions = ~d_info.excl & d_info.running==0;

% select data
this_data = seqCellsByTime_prop; %seqCellsByTime_prop;
this_data_switch = seqCellsByTime_switch_prop; %seqCellsByTime_switch_prop;

% preparations
these_rgbs = discretisedColourMap('winter',true,numBins);
this_n_lower = min(sum(these_sessions,1));
this_n_upper = max(sum(these_sessions,1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

nrows = 2; ncols = numBins+1; m=0;
F = default_figure();

for n=1:numBins
    m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0
            temp = nan(1,numBlocks);
            temp(find(these_blocks(idx,:))) = this_data(n,find(these_blocks(idx,:)),idx);
            plot(1:numBlocks,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    shadedErrorBar(1:numBlocks,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',these_rgbs(n,:));
    temp = cumsum(blocks(1:end-1))+0.5;
    for i=1:length(temp)
        if i==1 || i==3
            xline(temp(i),'k-');
        else
            xline(temp(i),'k:');
        end
    end
    xlim([0,numBlocks])
    xlabel('Trial block')
    ylim([0,ceil(nanmax(this_data(:))*100)])
    ytickformat('percentage')
    ylabel('Proportion of cells')
    title(binLabels{n})
end

m = m+1; subplot(nrows,ncols,m); hold on;
for n=1:numBins
    temp = nan(size(these_blocks));
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0
            temp(idx,find(these_blocks(idx,:))) = this_data(n,find(these_blocks(idx,:)),idx);
        end
    end
    shadedErrorBar(1:numBlocks,nanmean(temp,1)*100,nansem(temp,1)*100,'lineProps',these_rgbs(n,:));
end
temp = cumsum(blocks(1:end-1))+0.5;
for i=1:length(temp)
    if i==1 || i==3
        xline(temp(i),'k-');
    else
        xline(temp(i),'k:');
    end
end
xlim([0,numBlocks])
xlabel('Trial block')
% temp = nanmean(this_data(:,:,:),3);
% ylim([0,ceil(nanmax(temp(:))*100)])
ylim([0,ceil(nanmax(this_data(:))*100)])
ytickformat('percentage')
ylabel('Proportion of cells')
title('Days 1-5')

for n=1:numBins
    m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0
            temp = nan(1,numBlocks_switch);
            temp(find(these_blocks_switch(idx,:))) = this_data_switch(n,find(these_blocks_switch(idx,:)),idx);
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',these_rgbs(n,:));
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:)*100,'k');
    xline(blocks_switch(1)+0.5,'k:');
    xlim([0,numBlocks_switch])
    xlabel('Trial block')
    ylim([0,ceil(nanmax(this_data(:))*100)])
    ytickformat('percentage')
    ylabel('Proportion of cells')
    NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(NE(1),NE(2),...
        ['n = ',num2str(size(temp_data,1)),', rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title(binLabels{n})
end

m = m+1; subplot(nrows,ncols,m); hold on;
for n=1:numBins
    temp = nan(size(these_blocks_switch));
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0
            temp(idx,find(these_blocks_switch(idx,:))) = this_data_switch(n,find(these_blocks_switch(idx,:)),idx);
        end
    end
    shadedErrorBar(1:numBlocks_switch,nanmean(temp,1)*100,nansem(temp,1)*100,'lineProps',these_rgbs(n,:));
end
xline(blocks_switch(1)+0.5,'k:');
xlim([0,numBlocks_switch])
xlabel('Trial block')
% temp = nanmean(this_data(:,:,:),3);
% ylim([0,ceil(nanmax(temp(:))*100)])
ylim([0,ceil(nanmax(this_data(:))*100)])
ytickformat('percentage')
ylabel('Proportion of cells')
title('Switches averaged')

if this_n_lower==this_n_upper
    suptitle(['Sequence development across days (',these_animals,', n=',num2str(this_n_lower),')'])
else
    suptitle(['Sequence development across days (',these_animals,', n=',num2str(this_n_lower),'-',num2str(this_n_upper),')'])
end
drawnow;
savefig(F,[path.root_summary,'plots\sequenceDevelopment\sequenceDevelopment_',num2str(numBins),'_',these_animals,'.fig']);
saveas(F,[path.root_summary,'plots\sequenceDevelopment\sequenceDevelopment_',num2str(numBins),'_',these_animals,'.png']);




%% --- SUMMARY FIGURE --- Proportion of cells as a function of trial block - for SfN

% select data
%this_data = seqCellsByTime_prop; %seqCellsByTime_prop;
this_data_switch = seqCellsByTime_switch_prop; %seqCellsByTime_switch_prop;

% select sessions
these_animals = 'runners'; these_sessions = ~d_info.excl & d_info.running==1;
%these_animals = 'non-runners'; these_sessions = ~d_info.excl & d_info.running==0;

% preparations
%these_rgbs = discretisedColourMap('winter',true,numBins);
this_n_lower = min(sum(these_sessions,1));
this_n_upper = max(sum(these_sessions,1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

nrows = 3; ncols = numBins; m=0;
F = default_figure([-20,0.5,10,9.9]);

early_middle_late_runner = {};
for n=1:numBins
    m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0
            temp = nan(1,numBlocks_switch);
            temp(find(these_blocks_switch(idx,:))) = this_data_switch(n,find(these_blocks_switch(idx,:)),idx);
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.runner);
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:)*100,'k');
    xline(blocks_switch(1)+0.5,'k:');
    xlim([0,numBlocks_switch])
    xlabel('Trial block')
    ylim([0,ceil(nanmax(this_data(:))*100)])
    ytickformat('percentage')
    ylabel('Proportion of cells')
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),...
        ['n = ',num2str(size(temp_data,1)),', rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title(binLabels{n})
    early_middle_late_runner{n} = temp_data;
end


% select sessions
%these_animals = 'runners'; these_sessions = ~d_info.excl & d_info.running==1;
these_animals = 'non-runners'; these_sessions = ~d_info.excl & d_info.running==0;

% preparations
%these_rgbs = discretisedColourMap('winter',true,numBins);
this_n_lower = min(sum(these_sessions,1));
this_n_upper = max(sum(these_sessions,1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

early_middle_late_nonrunner = {};
for n=1:numBins
    m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0
            temp = nan(1,numBlocks_switch);
            temp(find(these_blocks_switch(idx,:))) = this_data_switch(n,find(these_blocks_switch(idx,:)),idx);
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.nonrunner);
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:)*100,'k');
    xline(blocks_switch(1)+0.5,'k:');
    xlim([0,numBlocks_switch])
    xlabel('Trial block')
    ylim([0,ceil(nanmax(this_data(:))*100)])
    ytickformat('percentage')
    ylabel('Proportion of cells')
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),...
        ['n = ',num2str(size(temp_data,1)),', rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title(binLabels{n})
    early_middle_late_nonrunner{n} = temp_data;
end

for n=1:numBins
    m = m+1; subplot(nrows,ncols,m); hold on;
    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_nonrunner{n},1)*100,nansem(early_middle_late_nonrunner{n},1)*100,'lineProps',p.col.nonrunner);
    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_runner{n},1)*100,nansem(early_middle_late_runner{n},1)*100,'lineProps',p.col.runner);
    xline(blocks_switch(1)+0.5,'k:');
    xlim([0,numBlocks_switch])
    xlabel('Trial block')
    % temp = nanmean(this_data(:,:,:),3);
    % ylim([0,ceil(nanmax(temp(:))*100)])
    ylim([0,ceil(nanmax(this_data(:))*100)])
    ytickformat('percentage')
    ylabel('Proportion of cells')
    title(binLabels{n})
end

%suptitle(['Sequence development across days'])
drawnow;
savefig(F,[path.root_summary,'plots\sequenceDevelopment\sequenceDevelopment_',num2str(numBins),'.fig']);
saveas(F,[path.root_summary,'plots\sequenceDevelopment\sequenceDevelopment_',num2str(numBins),'.png']);


%% Plot different types of correlation with performance, animal-by-animal
% 
% for idx=1:d_info.numAnimals
%     if ~isempty(d{idx,1})
%         
%         nrows = 2; ncols = 4; m=0;
%         default_figure();
%         
%         m = m+1; subplot(nrows,ncols,m); hold on;
%         yline(50,'k:');
%         this_data = correct;
%         plot(1:numBlocks,this_data(idx,:)*100,'Color','k','LineWidth',3)
%         temp = cumsum(blocks(1:end-1))+0.5;
%         for i=1:length(temp)
%             if i==1 || i==3
%                 xline(temp(i),'k-');
%             else
%                 xline(temp(i),'k:');
%             end
%         end
%         xlim([0,numBlocks])
%         xlabel('Trial block')
%         ylim([40,100])
%         ytickformat('percentage')
%         ylabel('Performance (%correct)')
%         title('Behaviour')
% 
%         m = m+1; subplot(nrows,ncols,m); hold on;
%         this_data = seqCellsByTime_prop;
%         for n=1:numBins
%             temp0 = flipud(eval('winter'));
%             [~,temp1] = discretize([1,numBins],size(temp0,1));
%             temp2 = discretize(1:numBins,temp1);
%             plot(1:numBlocks,this_data(n,:,idx)*100,'Color',temp0(temp2(n),:));
%         end
%         temp = cumsum(blocks(1:end-1))+0.5;
%         for i=1:length(temp)
%             if i==1 || i==3
%                 xline(temp(i),'k-');
%             else
%                 xline(temp(i),'k:');
%             end
%         end
%         xlim([0,numBlocks])
%         xlabel('Trial block')
%         ytickformat('percentage')
%         ylabel('Proportion of cells')
%         if numBins==3
%             legend('early','middle','late')
%         end
%         title('Sequence')
%         
%         m = m+1; subplot(nrows,ncols,m); hold on;
%         these_blocks = 1:5;
%         this_data_x = seqCellsByTime_prop;
%         this_data_y = correct;
%         for n=1:numBins
%             temp0 = flipud(eval('winter'));
%             [~,temp1] = discretize([1,numBins],size(temp0,1));
%             temp2 = discretize(1:numBins,temp1);
%             scatter(this_data_x(n,these_blocks,idx)*100,this_data_y(idx,these_blocks)*100,'MarkerEdgeColor',temp0(temp2(n),:));
%             fitLine(this_data_x(n,these_blocks,idx)'*100,this_data_y(idx,these_blocks)'*100,temp0(temp2(n),:));
%         end
%         xlim([0,inf])
%         xtickformat('percentage')
%         xlabel('Proportion of cells')
%         ylim([0,100])
%         ytickformat('percentage')
%         ylabel('Performance (%correct)')
%         title('Correlation (Expert: Day 1)')
%         
%                 
%         m = m+1; subplot(nrows,ncols,m); hold on;
%         these_blocks = 6:31;
%         this_data_x = seqCellsByTime_prop;
%         this_data_y = correct;
%         for n=1:numBins
%             temp0 = flipud(eval('winter'));
%             [~,temp1] = discretize([1,numBins],size(temp0,1));
%             temp2 = discretize(1:numBins,temp1);
%             scatter(this_data_x(n,these_blocks,idx)*100,this_data_y(idx,these_blocks)*100,'MarkerEdgeColor',temp0(temp2(n),:));
%             fitLine(this_data_x(n,these_blocks,idx)'*100,this_data_y(idx,these_blocks)'*100,temp0(temp2(n),:));
%         end
%         xlim([0,inf])
%         xtickformat('percentage')
%         xlabel('Proportion of cells')
%         ylim([0,100])
%         ytickformat('percentage')
%         ylabel('Performance (%correct)')
%         title('Correlation (Switches: Days 2-5)')
%         
%         m = m+1; subplot(nrows,ncols,m); hold on;
%         these_blocks = 6:18;
%         this_data_x = seqCellsByTime_prop;
%         this_data_y = correct;
%         for n=1:numBins
%             temp0 = flipud(eval('winter'));
%             [~,temp1] = discretize([1,numBins],size(temp0,1));
%             temp2 = discretize(1:numBins,temp1);
%             scatter(this_data_x(n,these_blocks,idx)*100,this_data_y(idx,these_blocks)*100,'MarkerEdgeColor',temp0(temp2(n),:));
%             fitLine(this_data_x(n,these_blocks,idx)'*100,this_data_y(idx,these_blocks)'*100,temp0(temp2(n),:));
%         end
%         xlim([0,inf])
%         xtickformat('percentage')
%         xlabel('Proportion of cells')
%         ylim([0,100])
%         ytickformat('percentage')
%         ylabel('Performance (%correct)')
%         title('Correlation (Switch 1: Days 2-3)')
%         
%         m = m+1; subplot(nrows,ncols,m); hold on;
%         these_blocks = 19:31;
%         this_data_x = seqCellsByTime_prop;
%         this_data_y = correct;
%         for n=1:numBins
%             temp0 = flipud(eval('winter'));
%             [~,temp1] = discretize([1,numBins],size(temp0,1));
%             temp2 = discretize(1:numBins,temp1);
%             scatter(this_data_x(n,these_blocks,idx)*100,this_data_y(idx,these_blocks)*100,'MarkerEdgeColor',temp0(temp2(n),:));
%             fitLine(this_data_x(n,these_blocks,idx)'*100,this_data_y(idx,these_blocks)'*100,temp0(temp2(n),:));
%         end
%         xlim([0,inf])
%         xtickformat('percentage')
%         xlabel('Proportion of cells')
%         ylim([0,100])
%         ytickformat('percentage')
%         ylabel('Performance (%correct)')
%         title('Correlation (Switch 2: Days 4-5)')
%         
%         m = m+1; subplot(nrows,ncols,m); hold on;
%         these_blocks = [6:13,19:26];
%         this_data_x = seqCellsByTime_prop;
%         this_data_y = correct;
%         for n=1:numBins
%             temp0 = flipud(eval('winter'));
%             [~,temp1] = discretize([1,numBins],size(temp0,1));
%             temp2 = discretize(1:numBins,temp1);
%             scatter(this_data_x(n,these_blocks,idx)*100,this_data_y(idx,these_blocks)*100,'MarkerEdgeColor',temp0(temp2(n),:));
%             fitLine(this_data_x(n,these_blocks,idx)'*100,this_data_y(idx,these_blocks)'*100,temp0(temp2(n),:));
%         end
%         xlim([0,inf])
%         xtickformat('percentage')
%         xlabel('Proportion of cells')
%         ylim([0,100])
%         ytickformat('percentage')
%         ylabel('Performance (%correct)')
%         title('Correlation (First switch days: Days 2 and 4)')
%         
%         m = m+1; subplot(nrows,ncols,m); hold on;
%         these_blocks = [14:18,27:31];
%         this_data_x = seqCellsByTime_prop;
%         this_data_y = correct;
%         for n=1:numBins
%             temp0 = flipud(eval('winter'));
%             [~,temp1] = discretize([1,numBins],size(temp0,1));
%             temp2 = discretize(1:numBins,temp1);
%             scatter(this_data_x(n,these_blocks,idx)*100,this_data_y(idx,these_blocks)*100,'MarkerEdgeColor',temp0(temp2(n),:));
%             fitLine(this_data_x(n,these_blocks,idx)'*100,this_data_y(idx,these_blocks)'*100,temp0(temp2(n),:));
%         end
%         xlim([0,inf])
%         xtickformat('percentage')
%         xlabel('Proportion of cells')
%         ylim([0,100])
%         ytickformat('percentage')
%         ylabel('Performance (%correct)')
%         title('Correlation (Second switch days: Days 3 and 5)')
% 
%         suptitle(['Sequences across days (',d_info.animals{idx},')'])
%     end
% end


%% Plot different types of correlation with performance, summary figure

nrows = 2; ncols = 2; m=0;
default_figure();

these_labels = {'early','middle','late'};
these_data = corr_aw_d1.seqCellsByTime_prop_correct;
temp0 = flipud(eval('winter'));
[~,temp1] = discretize([1,numBins],size(temp0,1));
temp2 = discretize(1:numBins,temp1);
these_cols = {temp0(temp2(1),:),temp0(temp2(2),:),temp0(temp2(3),:)};

m = m+1; subplot(nrows,ncols,m); hold on;
barsWithLines(these_labels,these_data,these_cols);
ylabel('Correlation (proportion of cells x %correct)')
title('Expert')

these_labels = {'early','middle','late'};
these_data = corr_aw_d23d45.seqCellsByTime_prop_correct;
temp0 = flipud(eval('winter'));
[~,temp1] = discretize([1,numBins],size(temp0,1));
temp2 = discretize(1:numBins,temp1);
these_cols = {temp0(temp2(1),:),temp0(temp2(2),:),temp0(temp2(3),:)};

m = m+1; subplot(nrows,ncols,m); hold on;
barsWithLines(these_labels,these_data,these_cols);
ylabel('Correlation (proportion of cells x %correct)')
title('Switches (within switch, then both switches averaged)')

these_labels = {'early','middle','late'};
these_data = corr_aw_d23.seqCellsByTime_prop_correct;
temp0 = flipud(eval('winter'));
[~,temp1] = discretize([1,numBins],size(temp0,1));
temp2 = discretize(1:numBins,temp1);
these_cols = {temp0(temp2(1),:),temp0(temp2(2),:),temp0(temp2(3),:)};

m = m+1; subplot(nrows,ncols,m); hold on;
barsWithLines(these_labels,these_data,these_cols);
ylabel('Correlation (proportion of cells x %correct)')
title('Switch 1')

these_labels = {'early','middle','late'};
these_data = corr_aw_d45.seqCellsByTime_prop_correct;
temp0 = flipud(eval('winter'));
[~,temp1] = discretize([1,numBins],size(temp0,1));
temp2 = discretize(1:numBins,temp1);
these_cols = {temp0(temp2(1),:),temp0(temp2(2),:),temp0(temp2(3),:)};

m = m+1; subplot(nrows,ncols,m); hold on;
barsWithLines(these_labels,these_data,these_cols);
ylabel('Correlation (proportion of cells x %correct)')
title('Switch 2')

suptitle(['Sequence -> Behaviour summary'])



%% --- SUMMARY FIGURE --- Warp as a function of trial block

% select data
this_data = warp_exp1_b; %seqCellsByTime_prop;
this_data_switch = warp_exp1_b_switch; %seqCellsByTime_switch_prop;

% select sessions
these_animals = 'nostim-r1-exp1b'; these_sessions = ~d_info.excl & d_info.running==1; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1; & d_info.running==1;

% select plotting range
%this_lower = -0.9; this_upper = 0.9; % power2a
%this_lower = -22; this_upper = 22; % power2b
%this_lower = -1; this_upper = 3; % power2c
%this_lower = 0; this_upper = 0.25; % exp1a
this_lower = -2; this_upper = 1; % exp1b

% preparations
this_n_lower = min(sum(these_sessions(:,2:5),1));
this_n_upper = max(sum(these_sessions(:,2:5),1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

nrows = 1; ncols = 2; m=0;
F = default_figure();

m = m+1; subplot(nrows,ncols,m); hold on;
temp_data = [];
for idx=1:d_info.numAnimals
    if sum(these_sessions(idx,:))>0
        temp = this_data(idx,:);
        plot(1:numBlocks,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
end
temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks);
shadedErrorBar(1:numBlocks,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.black);
temp = repmat(1:numBlocks,size(temp_data,1),1);
this_temp = temp(:,6:18);
this_temp_data = temp_data(:,6:18);
[this_corr_r_1,this_corr_p_1] = fitLine(this_temp(:),this_temp_data(:),'b');
this_temp = temp(:,19:31);
this_temp_data = temp_data(:,19:31);
[this_corr_r_2,this_corr_p_2] = fitLine(this_temp(:),this_temp_data(:),'b');
temp = cumsum(blocks(1:end-1))+0.5;
for i=1:length(temp)
    if i==1 || i==3
        xline(temp(i),'k-');
    else
        xline(temp(i),'k:');
    end
end
xlim([0,numBlocks])
xlabel('Trial block')
ylim([this_lower,this_upper])
ylabel('Coefficient')
NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(NE(1),NE(2),...
    ['switch 1: n = ',num2str(size(temp_data,1)),', rho = ',num2str(this_corr_r_1,2),', p = ',num2str(this_corr_p_1,2),newline,...
    'switch 2: n = ',num2str(size(temp_data,1)),', rho = ',num2str(this_corr_r_2,2),', p = ',num2str(this_corr_p_2,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);

m = m+1; subplot(nrows,ncols,m); hold on;
temp_data = [];
for i=1:d_info.numAnimals
    if sum(these_sessions(idx,:))>0
        temp = nan(1,size(these_blocks_switch,2));
        temp(1,find(these_blocks_switch(idx,:))) = this_data_switch(idx,find(these_blocks_switch(idx,:)));
        plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
end
shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.black);
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
[this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:),'k');
this_temp = temp(:,1:8);
this_temp_data = temp_data(:,1:8);
[this_corr_r_1,this_corr_p_1] = fitLine(this_temp(:),this_temp_data(:),'b');
this_temp = temp(:,9:13);
this_temp_data = temp_data(:,9:13);
[this_corr_r_2,this_corr_p_2] = fitLine(this_temp(:),this_temp_data(:),'b');
xline(blocks_switch(1)+0.5,'k:');
xlim([0,numBlocks_switch])
xlabel('Trial block')
ylim([this_lower,this_upper])
ylabel('Coefficient')
NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(NE(1),NE(2),...
    ['overall: n = ',num2str(size(temp_data,1)),', rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2),newline,...
    'switch 1: n = ',num2str(size(temp_data,1)),', rho = ',num2str(this_corr_r_1,2),', p = ',num2str(this_corr_p_1,2),newline,...
    'switch 2: n = ',num2str(size(temp_data,1)),', rho = ',num2str(this_corr_r_2,2),', p = ',num2str(this_corr_p_2,2),...
    ],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
ylabel('Proportion of cells')
title('Switches averaged')

suptitle(['Sequence development across days (',these_animals,')'])
drawnow;
savefig(F,[path.root_summary,'plots\warping\warping_',these_animals,'.fig']);
saveas(F,[path.root_summary,'plots\warping\warping_',these_animals,'.png']);



%% Warp vs performance

% select data
this_data_x = warp_exp1_b;
this_data_y = correct;

figure; hold on;
yline(50,'k:');
plot(this_data_x(:),this_data_y(:)*100,'ko')
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k');

xlim([-2,1])
ylim([30,100])
xlabel('Warping (b)')
ylabel('Performance (%correct)')
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['overall: rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);

title('Days 1-5 (31 data points per animal)')



%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! --- Warp -> performance correlation ---

F = default_figure([-20,0.5,10,3.5]); hold on;

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
plot(this_data_x(:),this_data_y(:)*100,'ko')
plot(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'o','Color',p.col.nonrunner)
plot(this_data_x_runners(:),this_data_y_runners(:)*100,'o','Color',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k');
%xlim([0,0.25])
ylim([0,100])
xlabel('Exponential coefficient a (initial proportion)')
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
plot(this_data_x(:),this_data_y(:)*100,'ko')
plot(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'o','Color',p.col.nonrunner)
plot(this_data_x_runners(:),this_data_y_runners(:)*100,'o','Color',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k');
%xlim([0,0.25])
ylim([0,100])
xlabel('Exponential coefficient b (decay rate)')
ylabel('Performance')
ytickformat('percentage')
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['overall: rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2),newline,...
    'non-runners: rho = ',num2str(this_corr_r_nonrunners,2),', p = ',num2str(this_corr_p_nonrunners,2),newline,...
    'runners: rho = ',num2str(this_corr_r_runners,2),', p = ',num2str(this_corr_p_runners,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);

drawnow;
savefig(F,[path.root_summary,'plots\warping\warping_corrWithBehaviour_',these_animals,'.fig']);
saveas(F,[path.root_summary,'plots\warping\warping_corrWithBehaviour_',these_animals,'.png']);









%% --- SUMMARY FIGURE --- Proportion of cells as a function of trial block - for SfN

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

nrows = 1; ncols = numBins; m=0;
F = default_figure([-20,0.5,10,3.5]);

early_middle_late_runner = {};
early_middle_late_runner_corr = {};
for n=1:numBins
%    m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if runners_switch_d12(idx)
            temp = nan(1,numBlocks_switch);
            temp(find(these_blocks_switch(idx,:))) = this_data_switch(n,find(these_blocks_switch(idx,:)),idx);
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.runner);
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:)*100,'k');
    xline(blocks_switch(1)+0.5,'k:');
    xlim([0,numBlocks_switch])
    xlabel('Trial block')
    ylim([0,15])
    ytickformat('percentage')
    ylabel('Proportion of cells')
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),...
        ['n = ',num2str(size(temp_data,1)),', rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title(binLabels{n})
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
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.nonrunner);
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:)*100,'k');
    xline(blocks_switch(1)+0.5,'k:');
    xlim([0,numBlocks_switch])
    xlabel('Trial block')
    ylim([0,15])
    ytickformat('percentage')
    ylabel('Proportion of cells')
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),...
        ['n = ',num2str(size(temp_data,1)),', rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title(binLabels{n})
    early_middle_late_nonrunner_corr{n} = [size(temp_data,1),this_corr_r,this_corr_p];
    early_middle_late_nonrunner{n} = temp_data;
end

for n=1:numBins
    m = m+1; subplot(nrows,ncols,m); hold on;
    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_nonrunner{n},1)*100,nansem(early_middle_late_nonrunner{n},1)*100,'lineProps',p.col.nonrunner);
    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_runner{n},1)*100,nansem(early_middle_late_runner{n},1)*100,'lineProps',p.col.runner);
    xline(blocks_switch(1)+0.5,'k:');
    xlim([0,numBlocks_switch])
    xlabel('Trial block')
    ylim([0,15])
    ytickformat('percentage')
    ylabel('Proportion of cells')
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    %NE = [max(xlim) max(ylim)-max(ylim)/5]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
%     text(NE(1),NE(2),...
%         ['non-runners: n = ',num2str(early_middle_late_nonrunner_corr{n}(1)),', rho = ',num2str(early_middle_late_nonrunner_corr{n}(2),2),', p = ',num2str(early_middle_late_nonrunner_corr{n}(3),2),newline,...
%         'runners: n = ',num2str(early_middle_late_runner_corr{n}(1)),', rho = ',num2str(early_middle_late_runner_corr{n}(2),2),', p = ',num2str(early_middle_late_runner_corr{n}(3),2)],...
%         'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    text(SE(1),SE(2),...
        ['n = ',num2str(early_middle_late_nonrunner_corr{n}(1)),', rho = ',num2str(early_middle_late_nonrunner_corr{n}(2),2),', p = ',num2str(early_middle_late_nonrunner_corr{n}(3),2),newline,...
        'n = ',num2str(early_middle_late_runner_corr{n}(1)),', rho = ',num2str(early_middle_late_runner_corr{n}(2),2),', p = ',num2str(early_middle_late_runner_corr{n}(3),2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title(binLabels{n})
end

%suptitle(['Sequence development across days'])
drawnow;
savefig(F,[path.root_summary,'plots\sequenceDevelopment\sequenceDevelopment_',num2str(numBins),'.fig']);
saveas(F,[path.root_summary,'plots\sequenceDevelopment\sequenceDevelopment_',num2str(numBins),'.png']);


%% First try with mixed-effect model

% early_middle_late_runner{3}
% early_middle_late_nonrunner{3}

this_data = warp_exp1_a(:,blocks(1)+1:end); % b
%this_data = squeeze(seqCellsByTime_prop(1,blocks(1)+1:end,:))'; % early sequence cells
%this_data = squeeze(seqCellsByTime_prop(3,blocks(1)+1:end,:))'; % late sequence cells

this_data = squeeze(seqCellsByTime_prop(1,blocks(1)+1:end,:))' - squeeze(seqCellsByTime_prop(3,blocks(1)+1:end,:))';

tbl_in.y = this_data(:);
temp = repmat([1:13,1:13],d_info.numAnimals,1);
tbl_in.block = temp(:);
temp = repmat([1:8,1:5,1:8,1:5],d_info.numAnimals,1);
tbl_in.block_sess = temp(:);
temp = repmat([0:1/(8-1):1,0:1/(5-1):1,0:1/(8-1):1,0:1/(5-1):1],d_info.numAnimals,1);
tbl_in.block_sess_norm = temp(:);
tbl_in.animal = nominal(repmat((1:d_info.numAnimals)',size(this_data,2),1));
temp = repmat([1*ones(1,13),2*ones(1,13)],d_info.numAnimals,1);
tbl_in.switch = nominal(temp(:));
temp = repmat([1*ones(1,8),2*ones(1,5),1*ones(1,8),2*ones(1,5)],d_info.numAnimals,1);
tbl_in.switchday = nominal(temp(:));
temp = nan(size(this_data));
for j=2:5
    if j==2
        temp(:,1:8) = repmat(d_info.running(:,j),1,8);
    elseif j==3
        temp(:,9:13) = repmat(d_info.running(:,j),1,5);
    elseif j==4
        temp(:,14:21) = repmat(d_info.running(:,j),1,8);
    elseif j==5
        temp(:,22:26) = repmat(d_info.running(:,j),1,5);
    end
end
tbl_in.running = nominal(temp(:));
tbl = table(tbl_in.animal,tbl_in.switch,tbl_in.switchday,tbl_in.block,tbl_in.block_sess,tbl_in.block_sess_norm,tbl_in.running,tbl_in.y,...
    'VariableNames',{'animal','switch','switchday','block','block_sess','block_sess_norm','running','y'});

lme = fitlme(tbl,'y ~ running*block + (1|animal) + (1|switch)')

% lme = fitlme(tbl,'y ~ running*switchday*block_sess + (1|animal) + (1|switch)')


%lme = fitlme(tbl,'y ~ running*block + (block|animal)')
%lme = fitlme(tbl,'y ~ running*block + (block|animal)')

%lme = fitlme(tbl,'y ~ running*block + (block|animal) + (block|switch)')


% var 1: individual mixed-effect model for each trial block (Rolotti)

% y: proportion of (early, middle, late) sequence cells / warping (a, b) / number of sequence cells
% fixed effects: trial block, running, trial block*running (stim sessions: , stim, stim*running, stim*trial block, stim*running*trial block)
% random effects: mouse, switch (1 or 2)


%% !!!!!!!!!!!!!!!!!!!! --- SUMMARY FIGURE --- Warp as a function of trial block

% select data
this_data_switch = warp_exp1_a_switch;
these_sessions = ~d_info.excl; %& d_info.running==0;
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

nrows = 1; ncols = 2; m=0;
F = default_figure([-20,0.5,7.5,3.5]);

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
xlabel('Trial block')
ylabel('Coefficient')
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['n = ',num2str(exponential_nonrunner_corr(1)),', rho = ',num2str(exponential_nonrunner_corr(2),2),', p = ',num2str(exponential_nonrunner_corr(3),2),newline,...
    'n = ',num2str(exponential_runner_corr(1)),', rho = ',num2str(exponential_runner_corr(2),2),', p = ',num2str(exponential_runner_corr(3),2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
% text(SE(1),SE(2),...
%     ['non-runner: n = ',num2str(exponential_nonrunner_corr(1)),', rho = ',num2str(exponential_nonrunner_corr(2),2),', p = ',num2str(exponential_nonrunner_corr(3),2),newline,...
%     'runner: n = ',num2str(exponential_runner_corr(1)),', rho = ',num2str(exponential_runner_corr(2),2),', p = ',num2str(exponential_runner_corr(3),2)],...
%     'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
ylabel('Exponential coefficient a (initial proportion)')

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
ylim([-1.4,0.2])
xlabel('Trial block')
ylabel('Coefficient')
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['n = ',num2str(exponential_nonrunner_corr(1)),', rho = ',num2str(exponential_nonrunner_corr(2),2),', p = ',num2str(exponential_nonrunner_corr(3),2),newline,...
    'n = ',num2str(exponential_runner_corr(1)),', rho = ',num2str(exponential_runner_corr(2),2),', p = ',num2str(exponential_runner_corr(3),2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
% text(SE(1),SE(2),...
%     ['non-runner: n = ',num2str(exponential_nonrunner_corr(1)),', rho = ',num2str(exponential_nonrunner_corr(2),2),', p = ',num2str(exponential_nonrunner_corr(3),2),newline,...
%     'runner: n = ',num2str(exponential_runner_corr(1)),', rho = ',num2str(exponential_runner_corr(2),2),', p = ',num2str(exponential_runner_corr(3),2)],...
%     'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
ylabel('Exponential coefficient b (decay rate)')

% suptitle(['Sequence development across days (',these_animals,')'])
drawnow;
savefig(F,[path.root_summary,'plots\warping\warping_',these_animals,'.fig']);
saveas(F,[path.root_summary,'plots\warping\warping_',these_animals,'.png']);


%% !!!!!!!!!!!!!!!!!!!! --- SUMMARY FIGURE --- NumSeqCells as a function of trial block

% select data
this_data_switch = numSeqCells_switch;
these_sessions = ~d_info.excl; %& d_info.running==0;
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

nrows = 1; ncols = 1; m=0;
F = default_figure([-20,0.5,5,3.5]);

yline(0,'k:')
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
ylim([0,200])
yticks([0:100:200])
xlabel('Trial block')
ylabel('Coefficient')
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['non-runner: n = ',num2str(exponential_nonrunner_corr(1)),', rho = ',num2str(exponential_nonrunner_corr(2),2),', p = ',num2str(exponential_nonrunner_corr(3),2),newline,...
    'runner 1: n = ',num2str(exponential_runner_corr(1)),', rho = ',num2str(exponential_runner_corr(2),2),', p = ',num2str(exponential_runner_corr(3),2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
ylabel('Number of sequence cells')

% suptitle(['Sequence development across days (',these_animals,')'])
drawnow;
savefig(F,'C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Supp\numSeqCells_switch.fig');
saveas(F,'C:\Users\Moritz\Documents\Illustrator\SfN2022\Material\Supp\numSeqCells_switch.png');



%% !!!!!!!!!!!!!!! Sequence fit across blocks - runners vs. non-runners

nrows = 1; ncols = numBlocks_switch; m=0;
F = default_figure([-20,0.5,20,2.5]);

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
    ylim([0,0.1])
    
    if m==1 || m==numBlocks_switch
        NE = [max(xlim) max(ylim)-max(ylim)/2]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
        text(NE(1),NE(2),...
            ['a = ',num2str(nanmean(temp_a_nonrunners),2),newline,'b = ',num2str(nanmean(temp_b_nonrunners),2),newline,...
            'a = ',num2str(nanmean(temp_a_runners),2),newline,'b = ',num2str(nanmean(temp_b_runners),2)],...
            'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    end

    xlabel('Time (s)')
    if m==1
        ylabel('Proportion of sequence cells')
    end
    if m~=1
        set(gca,'yticklabel',{[]})
    end
    title(['Block ',num2str(m)])
end
drawnow;
savefig(F,[path.root_summary,'plots\warping\warping_runnonrun_fitOverTime.fig']);
saveas(F,[path.root_summary,'plots\warping\warping_runnonrun_fitOverTime.png']);



