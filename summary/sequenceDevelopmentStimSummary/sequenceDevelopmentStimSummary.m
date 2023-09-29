%% Preparations

blocks = [5,8,5,8,5];
numBlocks = sum(blocks);

binEdges = [0,1.5,3,5.3]; %[0,1,2,3,4,5,6]; %[0:5/3:5.3]; %[0,1,3,5.4]; %[0:5.3/3:5.3]; %[0:1:5];
numBins = length(binEdges)-1;
binLabels = {'Early peak (0 s - 1.5 s)','Middle peak (1.5 s - 3 s)','Late peak (3 s - 5.3 s)'};
%binLabels = {'eSequence cells (0 s - 1 s)','eSequence cells (1 s - 3 s)','eSequence cells (3 s - 5.3 s)'};


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
        if isfield(d{i,j},'tng_100t_stimVersion') %~isempty(d{i,j})
            for k=1:length(d{i,j}.tng_100t_stimVersion)
                this_block = temp0+k;

                % get idcs
                this_numCells = length(find(d{i,j}.tng_100t_stimVersion{k}.prop.iscell==1));
                these_idcs_A = find(d{i,j}.tng_100t_stimVersion{k}.passed_catch.AW.Acatchonly==1);
                these_idcs_X = find(d{i,j}.tng_100t_stimVersion{k}.passed_catch.AW.Xcatchonly==1);

                % get peak locations
                these_peakTimes_A = d{i,j}.tng_100t_stimVersion{k}.firingField.Acatch_AW.peakLocation_s(these_idcs_A);
                these_peakTimes_X = d{i,j}.tng_100t_stimVersion{k}.firingField.Xcatch_AW.peakLocation_s(these_idcs_X);
                these_peakTimes = [these_peakTimes_A; these_peakTimes_X];

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
warp_other = nan(d_info.numAnimals,numBlocks,12); % [animal, block, metric]
warp_linIdx = nan(d_info.numAnimals,numBlocks,12); % [animal, block, metric]
warp_nonlinIdx = nan(d_info.numAnimals,numBlocks,12); % [animal, block, metric]
warp_poly1_rsquare = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_exp1_rsquare = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_exp1_a = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_exp1_b = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_power2_a = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_power2_b = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_power2_c = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_exp1const_a = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_exp1const_b = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_exp1const_c = nan(d_info.numAnimals,numBlocks); % [animal, block]
numSeqCells = nan(d_info.numAnimals,numBlocks);
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
                
                warp_exp1_a(i,this_block) = d{i,j}.(warp_struct){k}.(firingField_metric).exp1.a;
                warp_exp1_b(i,this_block) = d{i,j}.(warp_struct){k}.(firingField_metric).exp1.b;
                warp_power2_a(i,this_block) = d{i,j}.(warp_struct){k}.(firingField_metric).power2.a;
                warp_power2_b(i,this_block) = d{i,j}.(warp_struct){k}.(firingField_metric).power2.b;
                warp_power2_c(i,this_block) = d{i,j}.(warp_struct){k}.(firingField_metric).power2.c;
%                 warp_exp1const_a(i,this_block) = d{i,j}.(warp_struct){k}.exp1const.a;
%                 warp_exp1const_b(i,this_block) = d{i,j}.(warp_struct){k}.exp1const.b;
%                 warp_exp1const_c(i,this_block) = d{i,j}.(warp_struct){k}.exp1const.c;
                try
                    numSeqCells(i,this_block) = d{i,j}.(warp_struct){k}.input.numSequenceCells;
                catch
                end
                
                warp_poly1_rsquare(i,this_block) = d{i,j}.(warp_struct){k}.(firingField_metric).gof.poly1.rsquare;
                warp_exp1_rsquare(i,this_block) = d{i,j}.(warp_struct){k}.(firingField_metric).gof.exp1.rsquare;
                
                warp_other(i,this_block,1) = d{i,j}.(warp_struct){k}.(firingField_metric).gof.poly1.rsquare;
                warp_other(i,this_block,2) = d{i,j}.(warp_struct){k}.(firingField_metric).gof.exp1.rsquare;
                warp_other(i,this_block,3) = d{i,j}.(warp_struct){k}.(firingField_metric).gof.poly1.rmse.^2;
                warp_other(i,this_block,4) = d{i,j}.(warp_struct){k}.(firingField_metric).gof.exp1.rmse.^2;
                warp_other(i,this_block,5) = d{i,j}.(warp_struct){k}.(firingField_metric).poly1.p1;
                warp_other(i,this_block,6) = d{i,j}.(warp_struct){k}.(firingField_metric).poly1.p2;
                warp_other(i,this_block,7) = d{i,j}.(warp_struct){k}.(firingField_metric).exp1.a;
                warp_other(i,this_block,8) = d{i,j}.(warp_struct){k}.(firingField_metric).exp1.b;
                
                warp_linIdx(i,this_block,1) = d{i,j}.(warp_struct){k}.(firingField_metric).linIdx1;
                warp_linIdx(i,this_block,2) = d{i,j}.(warp_struct){k}.(firingField_metric).linIdx2;
                warp_linIdx(i,this_block,3) = d{i,j}.(warp_struct){k}.(firingField_metric).linIdx3;
                warp_linIdx(i,this_block,4) = d{i,j}.(warp_struct){k}.(firingField_metric).linIdx4;
                warp_linIdx(i,this_block,5) = d{i,j}.(warp_struct){k}.(firingField_metric).linIdx5;
                warp_linIdx(i,this_block,6) = d{i,j}.(warp_struct){k}.(firingField_metric).linIdx6;
                warp_linIdx(i,this_block,7) = d{i,j}.(warp_struct){k}.(firingField_metric).linIdx7;
                warp_linIdx(i,this_block,8) = d{i,j}.(warp_struct){k}.(firingField_metric).linIdx8;
                warp_linIdx(i,this_block,9) = d{i,j}.(warp_struct){k}.(firingField_metric).linIdx9;
                warp_linIdx(i,this_block,10) = d{i,j}.(warp_struct){k}.(firingField_metric).linIdx10;
                warp_linIdx(i,this_block,11) = d{i,j}.(warp_struct){k}.(firingField_metric).linIdx11;
                warp_linIdx(i,this_block,12) = d{i,j}.(warp_struct){k}.(firingField_metric).linIdx12;
                warp_nonlinIdx(i,this_block,1) = d{i,j}.(warp_struct){k}.(firingField_metric).nonLinIdx1;
                warp_nonlinIdx(i,this_block,2) = d{i,j}.(warp_struct){k}.(firingField_metric).nonLinIdx2;
                warp_nonlinIdx(i,this_block,3) = d{i,j}.(warp_struct){k}.(firingField_metric).nonLinIdx3;
                warp_nonlinIdx(i,this_block,4) = d{i,j}.(warp_struct){k}.(firingField_metric).nonLinIdx4;
                warp_nonlinIdx(i,this_block,5) = d{i,j}.(warp_struct){k}.(firingField_metric).nonLinIdx5;
                warp_nonlinIdx(i,this_block,6) = d{i,j}.(warp_struct){k}.(firingField_metric).nonLinIdx6;
                warp_nonlinIdx(i,this_block,7) = d{i,j}.(warp_struct){k}.(firingField_metric).nonLinIdx7;
                warp_nonlinIdx(i,this_block,8) = d{i,j}.(warp_struct){k}.(firingField_metric).nonLinIdx8;
                warp_nonlinIdx(i,this_block,9) = d{i,j}.(warp_struct){k}.(firingField_metric).nonLinIdx9;
                warp_nonlinIdx(i,this_block,10) = d{i,j}.(warp_struct){k}.(firingField_metric).nonLinIdx10;
                warp_nonlinIdx(i,this_block,11) = d{i,j}.(warp_struct){k}.(firingField_metric).nonLinIdx11;
                warp_nonlinIdx(i,this_block,12) = d{i,j}.(warp_struct){k}.(firingField_metric).nonLinIdx12;
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

corr_correct_catch_exp1_b_all = nan(size(correct_catch_all,1),1);
for i=1:length(corr_correct_catch_exp1_b_all)
    corr_correct_catch_exp1_b_all(i) = corr(correct_catch_all(i,:)',warp_exp1_b_all(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_b_switch1 = nan(size(correct_catch_switch1,1),1);
for i=1:length(corr_correct_catch_exp1_b_switch1)
    corr_correct_catch_exp1_b_switch1(i) = corr(correct_catch_switch1(i,:)',warp_exp1_b_switch1(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_b_switch2 = nan(size(correct_catch_switch2,1),1);
for i=1:length(corr_correct_catch_exp1_b_switch2)
    corr_correct_catch_exp1_b_switch2(i) = corr(correct_catch_switch2(i,:)',warp_exp1_b_switch2(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_b_switch1_seq = nan(size(correct_catch_switch1_seq,1),1);
for i=1:length(corr_correct_catch_exp1_b_switch1_seq)
    corr_correct_catch_exp1_b_switch1_seq(i) = corr(correct_catch_switch1_seq(i,:)',warp_exp1_b_switch1_seq(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_b_switch1_ctrl = nan(size(correct_catch_switch1_ctrl,1),1);
for i=1:length(corr_correct_catch_exp1_b_switch1_ctrl)
    corr_correct_catch_exp1_b_switch1_ctrl(i) = corr(correct_catch_switch1_ctrl(i,:)',warp_exp1_b_switch1_ctrl(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_b_switch2_seq = nan(size(correct_catch_switch2_seq,1),1);
for i=1:length(corr_correct_catch_exp1_b_switch2_seq)
    corr_correct_catch_exp1_b_switch2_seq(i) = corr(correct_catch_switch2_seq(i,:)',warp_exp1_b_switch2_seq(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_b_switch2_ctrl = nan(size(correct_catch_switch2_ctrl,1),1);
for i=1:length(corr_correct_catch_exp1_b_switch2_ctrl)
    corr_correct_catch_exp1_b_switch2_ctrl(i) = corr(correct_catch_switch2_ctrl(i,:)',warp_exp1_b_switch2_ctrl(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_b_seq = nan(size(correct_catch_seq,1),1);
for i=1:length(corr_correct_catch_exp1_b_seq)
    corr_correct_catch_exp1_b_seq(i) = corr(correct_catch_seq(i,:)',warp_exp1_b_seq(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_b_ctrl = nan(size(correct_catch_ctrl,1),1);
for i=1:length(corr_correct_catch_exp1_b_ctrl)
    corr_correct_catch_exp1_b_ctrl(i) = corr(correct_catch_ctrl(i,:)',warp_exp1_b_ctrl(i,:)','Type','Pearson','Rows','Complete');
end

corr_correct_catch_exp1_a_all = nan(size(correct_catch_all,1),1);
for i=1:length(corr_correct_catch_exp1_a_all)
    corr_correct_catch_exp1_a_all(i) = corr(correct_catch_all(i,:)',warp_exp1_a_all(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_a_switch1 = nan(size(correct_catch_switch1,1),1);
for i=1:length(corr_correct_catch_exp1_a_switch1)
    corr_correct_catch_exp1_a_switch1(i) = corr(correct_catch_switch1(i,:)',warp_exp1_a_switch1(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_a_switch2 = nan(size(correct_catch_switch2,1),1);
for i=1:length(corr_correct_catch_exp1_a_switch2)
    corr_correct_catch_exp1_a_switch2(i) = corr(correct_catch_switch2(i,:)',warp_exp1_a_switch2(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_a_switch1_seq = nan(size(correct_catch_switch1_seq,1),1);
for i=1:length(corr_correct_catch_exp1_a_switch1_seq)
    corr_correct_catch_exp1_a_switch1_seq(i) = corr(correct_catch_switch1_seq(i,:)',warp_exp1_a_switch1_seq(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_a_switch1_ctrl = nan(size(correct_catch_switch1_ctrl,1),1);
for i=1:length(corr_correct_catch_exp1_a_switch1_ctrl)
    corr_correct_catch_exp1_a_switch1_ctrl(i) = corr(correct_catch_switch1_ctrl(i,:)',warp_exp1_a_switch1_ctrl(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_a_switch2_seq = nan(size(correct_catch_switch2_seq,1),1);
for i=1:length(corr_correct_catch_exp1_a_switch2_seq)
    corr_correct_catch_exp1_a_switch2_seq(i) = corr(correct_catch_switch2_seq(i,:)',warp_exp1_a_switch2_seq(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_a_switch2_ctrl = nan(size(correct_catch_switch2_ctrl,1),1);
for i=1:length(corr_correct_catch_exp1_a_switch2_ctrl)
    corr_correct_catch_exp1_a_switch2_ctrl(i) = corr(correct_catch_switch2_ctrl(i,:)',warp_exp1_a_switch2_ctrl(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_a_seq = nan(size(correct_catch_seq,1),1);
for i=1:length(corr_correct_catch_exp1_a_seq)
    corr_correct_catch_exp1_a_seq(i) = corr(correct_catch_seq(i,:)',warp_exp1_a_seq(i,:)','Type','Pearson','Rows','Complete');
end
corr_correct_catch_exp1_a_ctrl = nan(size(correct_catch_ctrl,1),1);
for i=1:length(corr_correct_catch_exp1_a_ctrl)
    corr_correct_catch_exp1_a_ctrl(i) = corr(correct_catch_ctrl(i,:)',warp_exp1_a_ctrl(i,:)','Type','Pearson','Rows','Complete');
end



%% Plot different types of correlation with performance, animal-by-animal
% 
% for idx=1:d_info.numAnimals
%     if any(~d_info.excl(idx,2:5))
%         
%         nrows = 1; ncols = 2; m=0;
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
%         suptitle(['Sequences across days (',d_info.animals{idx},', group ',num2str(d_info.group(idx)),')'])
%     end
% end


%% --- SUMMARY FIGURE --- Proportion of cells as a function of trial block

% select data
this_data = seqCellsByTime_prop; %seqCellsByTime_prop;

% select sessions
these_animals = 'stim'; these_sessions = ~d_info.excl; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1;

% preparations
these_rgbs = discretisedColourMap('winter',true,numBins);
this_n_lower = min(sum(these_sessions(:,2:5),1));
this_n_upper = max(sum(these_sessions(:,2:5),1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

nrows = 2; ncols = numBins+1; m=0;
F = default_figure();

early_middle_late_seq = {};
for n=1:numBins
    m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
            temp = this_data(n,6:18,idx);
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
            temp = this_data(n,19:31,idx);
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
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
        ['n = ',num2str(size(temp_data,1)),', \rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title([binLabels{n},', seq stim'])
    early_middle_late_seq{n} = temp_data;
end

m = m+1; subplot(nrows,ncols,m); hold on;
for n=1:numBins
    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_seq{n},1)*100,nansem(early_middle_late_seq{n},1)*100,'lineProps',these_rgbs(n,:));
end
xline(blocks_switch(1)+0.5,'k:');
xlim([0,numBlocks_switch])
xlabel('Trial block')
ylim([0,ceil(nanmax(this_data(:))*100)])
ytickformat('percentage')
ylabel('Proportion of cells')
title('Switches averaged, seq stim')

early_middle_late_ctrl = {};
for n=1:numBins
    m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
            temp = this_data(n,6:18,idx);
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
            temp = this_data(n,19:31,idx);
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
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
        ['n = ',num2str(size(temp_data,1)),', \rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title([binLabels{n},', ctrl stim'])
    early_middle_late_ctrl{n} = temp_data;
end

m = m+1; subplot(nrows,ncols,m); hold on;
for n=1:numBins
    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_ctrl{n},1)*100,nansem(early_middle_late_ctrl{n},1)*100,'lineProps',these_rgbs(n,:));
end
xline(blocks_switch(1)+0.5,'k:');
xlim([0,numBlocks_switch])
xlabel('Trial block')
ylim([0,ceil(nanmax(this_data(:))*100)])
ytickformat('percentage')
ylabel('Proportion of cells')
title('Switches averaged, ctrl stim')

suptitle(['Sequence development (catch trials) across days (',these_animals,')'])
drawnow;
savefig(F,[path.root_summary,'plots\sequenceDevelopmentStim\sequenceDevelopmentStim_',num2str(numBins),'_',these_animals,'.fig']);
saveas(F,[path.root_summary,'plots\sequenceDevelopmentStim\sequenceDevelopmentStim_',num2str(numBins),'_',these_animals,'.png']);


%% --- SUMMARY FIGURE --- Proportion of cells as a function of trial block - for SfN

% select data
this_data = seqCellsByTime_prop; %seqCellsByTime_prop;

% select sessions
these_animals = 'stim'; these_sessions = ~d_info.excl; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1;

% preparations
%these_rgbs = discretisedColourMap('winter',true,numBins);
this_n_lower = min(sum(these_sessions(:,2:5),1));
this_n_upper = max(sum(these_sessions(:,2:5),1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

nrows = 3; ncols = numBins; m=0;
F = default_figure([-20,0.5,10,9.9]);

early_middle_late_seq = {};
for n=1:numBins
    m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
            temp = this_data(n,6:18,idx);
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
            temp = this_data(n,19:31,idx);
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.seq);
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
        ['n = ',num2str(size(temp_data,1)),', \rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title([binLabels{n}])
    early_middle_late_seq{n} = temp_data;
end

early_middle_late_ctrl = {};
for n=1:numBins
    m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
            temp = this_data(n,6:18,idx);
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
            temp = this_data(n,19:31,idx);
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.ctrl);
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:)*100,'k');
    xline(blocks_switch(1)+0.5,'k:');
    xlim([0,numBlocks_switch])
    xlabel('Trial block')
    ylim([0,ceil(nanmax(this_data(:))*100)])
    ytickformat('percentage')
    ylabel('Proportion of cells')
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; % NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),...
        ['n = ',num2str(size(temp_data,1)),', \rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title([binLabels{n}])
    early_middle_late_ctrl{n} = temp_data;
end

for n=1:numBins
    m = m+1; subplot(nrows,ncols,m); hold on;
    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_ctrl{n},1)*100,nansem(early_middle_late_ctrl{n},1)*100,'lineProps',p.col.ctrl);
    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_seq{n},1)*100,nansem(early_middle_late_seq{n},1)*100,'lineProps',p.col.seq);
    xline(blocks_switch(1)+0.5,'k:');
    xlim([0,numBlocks_switch])
    xlabel('Trial block')
    ylim([0,ceil(nanmax(this_data(:))*100)])
    ytickformat('percentage')
    ylabel('Proportion of cells')
    title([binLabels{n}])
end

suptitle(['Sequence development (catch trials) across days (',these_animals,')'])
drawnow;
savefig(F,[path.root_summary,'plots\sequenceDevelopmentStim\sequenceDevelopmentStim_',num2str(numBins),'_',these_animals,'.fig']);
saveas(F,[path.root_summary,'plots\sequenceDevelopmentStim\sequenceDevelopmentStim_',num2str(numBins),'_',these_animals,'.png']);


%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! --- SUMMARY FIGURE --- Proportion of cells as a function of trial block - for SfN - NEW

% select data
this_data = seqCellsByTime_prop; %seqCellsByTime_prop;

% select sessions
these_animals = 'stim'; these_sessions = ~d_info.excl; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1;

% preparations
%these_rgbs = discretisedColourMap('winter',true,numBins);
this_n_lower = min(sum(these_sessions(:,2:5),1));
this_n_upper = max(sum(these_sessions(:,2:5),1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

nrows = 1; ncols = numBins; m=0;
F = default_figure([-20,0.5,15,5]);

early_middle_late_seq = {};
early_middle_late_seq_corr{n} = {};
for n=1:numBins
    %m = m+1; subplot(nrows,ncols,m); hold on;
    temp_data = [];
    for idx=1:d_info.numAnimals
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
            temp = this_data(n,6:18,idx);
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
            temp = this_data(n,19:31,idx);
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.seq);
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
        ['n = ',num2str(size(temp_data,1)),', \rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title([binLabels{n}])
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
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
        if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
            temp = this_data(n,19:31,idx);
            plot(1:numBlocks_switch,temp*100,'Color',p.col.gray);
            temp_data = [temp_data; temp];
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
    shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1)*100,nansem(temp_data,1)*100,'lineProps',p.col.ctrl);
    temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
    [this_corr_r,this_corr_p] = fitLine(temp(:),temp_data(:)*100,'k');
    xline(blocks_switch(1)+0.5,'k:');
    xlim([0,numBlocks_switch])
    xlabel('Trial block')
    ylim([0,ceil(nanmax(this_data(:))*100)])
    ytickformat('percentage')
    ylabel('Proportion of cells')
    SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; % NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(SE(1),SE(2),...
        ['n = ',num2str(size(temp_data,1)),', \rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title([binLabels{n}])
    early_middle_late_ctrl{n} = temp_data;
    early_middle_late_ctrl_corr{n} = [size(temp_data,1),this_corr_r,this_corr_p];
end

for n=1:numBins
    m = m+1; subplot(nrows,ncols,m); hold on;
    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_ctrl{n},1)*100,nansem(early_middle_late_ctrl{n},1)*100,'lineProps',p.col.ctrl);
    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_seq{n},1)*100,nansem(early_middle_late_seq{n},1)*100,'lineProps',p.col.seq);
    xline(blocks_switch(1)+0.5,'k:');
    xlim([0,numBlocks_switch])
    xlabel('Trial block')
    ylim([0,7])%ylim([0,ceil(nanmax(this_data(:))*100)])
    ytickformat('percentage')
    ylabel('Proportion of cells')
    NE = [max(xlim) max(ylim)-max(ylim)/8]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(NE(1),NE(2),...
        ['seq: n = ',num2str(early_middle_late_seq_corr{n}(1)),', \rho = ',num2str(early_middle_late_seq_corr{n}(2),2),', p = ',num2str(early_middle_late_seq_corr{n}(3),2),newline,...
        'ctrl: n = ',num2str(early_middle_late_ctrl_corr{n}(1)),', \rho = ',num2str(early_middle_late_ctrl_corr{n}(2),2),', p = ',num2str(early_middle_late_ctrl_corr{n}(3),2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
    title([binLabels{n}])
end

%suptitle(['Sequence development (catch trials) across days (',these_animals,')'])
drawnow;
savefig(F,[path.root_summary,'plots\sequenceDevelopmentStim\sequenceDevelopmentStim_',num2str(numBins),'_',these_animals,'_new.fig']);
saveas(F,[path.root_summary,'plots\sequenceDevelopmentStim\sequenceDevelopmentStim_',num2str(numBins),'_',these_animals,'_new.png']);


%% Compare early, middle and late sequence cells in seq stim vs ctrl stim

nrows = 1; ncols = 3; m=0;
F = default_figure();

for n=1:numBins
    m = m+1; subplot(nrows,ncols,m); hold on;

    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_seq{n},1)*100,nansem(early_middle_late_seq{n},1)*100,'lineProps',p.col.seq);
    shadedErrorBar(1:numBlocks_switch,nanmean(early_middle_late_ctrl{n},1)*100,nansem(early_middle_late_ctrl{n},1)*100,'lineProps',p.col.ctrl);
    xline(blocks_switch(1)+0.5,'k:');
    xlim([0,numBlocks_switch])
    xlabel('Trial block')
    ylim([0,ceil(nanmax(this_data(:))*100)])
    ytickformat('percentage')
    ylabel('Proportion of cells')
    title([binLabels{n}])
end

suptitle(['Sequence development across days (',these_animals,')'])
drawnow;
savefig(F,[path.root_summary,'plots\sequenceDevelopmentStim\sequenceDevelopmentStim_seqctrl_',num2str(numBins),'_',these_animals,'.fig']);
saveas(F,[path.root_summary,'plots\sequenceDevelopmentStim\sequenceDevelopmentStim_seqctrl_',num2str(numBins),'_',these_animals,'.png']);




%% --- SUMMARY FIGURE --- Warp as a function of trial block

% select data
this_data = warp_exp1_b; %seqCellsByTime_prop;
%this_data(15,:) = NaN;

% select sessions
these_animals = 'stim-exp1b'; these_sessions = ~d_info.excl; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1; & d_info.running==1

% select plotting range
% this_lower = -0.9; this_upper = 0.9; % power2a
% this_lower = -22; this_upper = 22; % power2b
% this_lower = -1; this_upper = 1; % power2c
%this_lower = 0; this_upper = 0.3; % exp1a
this_lower = -2; this_upper = 0.2; % exp1b

% preparations
this_n_lower = min(sum(these_sessions(:,2:5),1));
this_n_upper = max(sum(these_sessions(:,2:5),1));
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

nrows = 1; ncols = 2; m=0;
F = default_figure([-20,0.5,15,7.5]);

m = m+1; subplot(nrows,ncols,m); hold on;
temp_data = [];
for idx=1:d_info.numAnimals
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
        temp = this_data(idx,6:18);
        plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
        temp = this_data(idx,19:31);
        plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
end
temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.seq);
temp = repmat(1:numBlocks_switch,size(temp_data,1),1);
% temp0 = temp(:,1:8); temp_data0 = temp_data(:,1:8);
% [this_corr_r_1,this_corr_p_1] = corr(temp0(:),temp_data0(:),'Type','Pearson','Rows','Complete')
% temp0 = temp(:,9:13); temp_data0 = temp_data(:,9:13);
% [this_corr_r_2,this_corr_p_2] = corr(temp0(:),temp_data0(:),'Type','Pearson','Rows','Complete')
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
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['overall: n = ',num2str(size(temp_data,1)),', \rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2),newline,...
    'switch day 1: n = ',num2str(size(temp_data,1)),', \rho = ',num2str(this_corr_r_1,2),', p = ',num2str(this_corr_p_1,2),newline,...
    'switch day 2: n = ',num2str(size(temp_data,1)),', \rho = ',num2str(this_corr_r_2,2),', p = ',num2str(this_corr_p_2,2),...
    ],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
title('seq stim')

m = m+1; subplot(nrows,ncols,m); hold on;
temp_data = [];
for idx=1:d_info.numAnimals
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==8
        temp = this_data(idx,6:18);
        plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
    if sum(these_sessions(idx,:))>0 && d_info.group(idx)==7
        temp = this_data(idx,19:31);
        plot(1:numBlocks_switch,temp,'Color',p.col.gray);
        temp_data = [temp_data; temp];
    end
end
temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
shadedErrorBar(1:numBlocks_switch,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.ctrl);
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
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['overall: n = ',num2str(size(temp_data,1)),', \rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2),newline,...
    'switch day 1: n = ',num2str(size(temp_data,1)),', \rho = ',num2str(this_corr_r_1,2),', p = ',num2str(this_corr_p_1,2),newline,...
    'switch day 2: n = ',num2str(size(temp_data,1)),', \rho = ',num2str(this_corr_r_2,2),', p = ',num2str(this_corr_p_2,2),...
    ],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);
title('ctrl stim')

suptitle(['Sequence development (catch trials) across days (',these_animals,')'])
drawnow;
savefig(F,[path.root_summary,'plots\warping\warping_seqctrl_',these_animals,'_temp.fig']);
saveas(F,[path.root_summary,'plots\warping\warping_seqctrl_',these_animals,'_temp.png']);



%% !!!!!!!!!!!!!!!!! --- SUMMARY FIGURE --- Warp as a function of trial block

% select sessions
these_animals = 'stim'; these_sessions = ~d_info.excl; % & repmat(d_info.group,1,d_info.numDays)==7 % & d_info.presponsive==1; & d_info.running==1
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

% select data
this_data = warp_exp1_a; %seqCellsByTime_prop;

nrows = 1; ncols = 2; m=0;
F = default_figure([-20,0.5,15,5]);
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
xlabel('Trial block')
ylim([0,0.12])
ylabel('Exponential coefficient a (initial proportion)')
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
xlabel('Trial block')
ylim([-0.7,0])
ylabel('Exponential coefficient b (decay rate)')
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['seq: n = ',num2str(this_corr_n_seq),', \rho = ',num2str(this_corr_r_seq,2),', p = ',num2str(this_corr_p_seq,2),newline,...
    'ctrl: n = ',num2str(this_corr_n_ctrl),', \rho = ',num2str(this_corr_r_ctrl,2),', p = ',num2str(this_corr_p_ctrl,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);


% title('seq stim')

% suptitle(['Sequence development (catch trials) across days (',these_animals,')'])
drawnow;
savefig(F,[path.root_summary,'plots\warping\warping_seqctrl_exponential.fig']);
saveas(F,[path.root_summary,'plots\warping\warping_seqctrl_exponential.png']);







%%

% exp1 model

% a
coef.seq_a_start = 0.08;
coef.seq_a_end = 0.055;
coef.ctrl_a_start = 0.07;
coef.ctrl_a_end = 0.07;

% b
coef.seq_b_start = -0.37;
coef.seq_b_end = -0.14;
coef.ctrl_b_start = -0.27;
coef.ctrl_b_end = -0.27;

x = 0:0.01:5.3;
ypred.seq_start = coef.seq_a_start*exp(x*coef.seq_b_start);
ypred.seq_end = coef.seq_a_end*exp(x*coef.seq_b_end);
ypred.ctrl_start = coef.ctrl_a_start*exp(x*coef.ctrl_b_start);
ypred.ctrl_end = coef.ctrl_a_end*exp(x*coef.ctrl_b_end);

figure

subplot(2,2,1); hold on
plot(x,ypred.seq_start,'Color',p.col.seq_rainbow(1,:),'LineWidth',3)
plot(x,ypred.seq_end,'Color',p.col.seq_rainbow(end,:),'LineWidth',3)
xlim([0.2,5.2])
legend('seq stim (start)','seq stim (end)')
title('seq stim: start -> end')

subplot(2,2,2); hold on
plot(x,ypred.ctrl_start,'Color',p.col.ctrl_rainbow(1,:),'LineWidth',3)
plot(x,ypred.ctrl_end,'Color',p.col.ctrl_rainbow(end,:),'LineWidth',3)
xlim([0.2,5.2])
legend('ctrl stim (start)','ctrl stim (end)')
title('ctrl stim: start -> end')

subplot(2,2,3); hold on
plot(x,ypred.seq_start,'Color',p.col.seq_rainbow(1,:),'LineWidth',3)
plot(x,ypred.ctrl_start,'Color',p.col.ctrl_rainbow(1,:),'LineWidth',3)
xlim([0.2,5.2])
legend('seq stim (start)','ctrl stim (start)')
title('session start: seq stim vs ctrl stim')

subplot(2,2,4); hold on
plot(x,ypred.seq_end,'Color',p.col.seq_rainbow(end,:),'LineWidth',3)
plot(x,ypred.ctrl_end,'Color',p.col.ctrl_rainbow(end,:),'LineWidth',3)
xlim([0.2,5.2])
legend('seq stim (end)','ctrl stim (end)')
title('session end: seq stim vs ctrl stim')

suptitle('Sequence warping')


%% !!!!!!!!!!!!!!! Sequence fit across blocks - seq stim vs ctrl stim

nrows = 1; ncols = numBlocks_switch; m=0;
F = default_figure([-20,0.5,20,5]);

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
            %plot(x,this_ypred,'Color',mean([mean([p.col.ctrl;p.col.white]);p.col.white]));
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
            %plot(x,this_ypred,'Color',mean([mean([p.col.ctrl;p.col.white]);p.col.white]));
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
    shadedErrorBar(x,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.ctrl);
    
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
            %plot(x,this_ypred,'Color',mean([mean([p.col.seq;p.col.white]);p.col.white]));
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
            %plot(x,this_ypred,'Color',mean([mean([p.col.seq;p.col.white]);p.col.white]));
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
    shadedErrorBar(x,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',p.col.seq);
    xlim([0,5.3])
    ylim([0,0.1])
    
    NE = [max(xlim) max(ylim)-max(ylim)/4]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
    text(NE(1),NE(2),...
        ['a = ',num2str(nanmean(temp_a_seq),2),newline,'b = ',num2str(nanmean(temp_b_seq),2),newline,newline,...
        'a = ',num2str(nanmean(temp_a_ctrl),2),newline,'b = ',num2str(nanmean(temp_b_ctrl),2)],...
        'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);


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
savefig(F,[path.root_summary,'plots\warping\warping_seqctrl_fitOverTime_',these_animals,'.fig']);
saveas(F,[path.root_summary,'plots\warping\warping_seqctrl_fitOverTime_',these_animals,'.png']);



%% Sequence fit across blocks - seq stim

these_rgbs = discretisedColourMap('jet',true,8);

nrows = 1; ncols = 2; m=0;
F = default_figure();

x = 0:0.01:5.3;    
subplot(nrows,ncols,1); hold on;
for m=1:8%numBlocks_switch

    temp_data = [];
    for idx=1:d_info.numAnimals
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
            %plot(x,this_ypred,'Color',mean([mean([p.col.ctrl;p.col.white]);p.col.white]));
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
            %plot(x,this_ypred,'Color',mean([mean([p.col.ctrl;p.col.white]);p.col.white]));
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
    plot(x,nanmean(temp_data,1),'Color',these_rgbs(m,:))
    %shadedErrorBar(x,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',these_rgbs(m,:));
    
    xlim([0,5])
    ylim([0,0.1])
end

n = 0;
subplot(nrows,ncols,2); hold on;
for m=9:13 %numBlocks_switch
    n = n+1;

    temp_data = [];
    for idx=1:d_info.numAnimals
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
            %plot(x,this_ypred,'Color',mean([mean([p.col.ctrl;p.col.white]);p.col.white]));
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
            %plot(x,this_ypred,'Color',mean([mean([p.col.ctrl;p.col.white]);p.col.white]));
        end
    end
    temp_data = rmmissing(temp_data,'MinNumMissing',numBlocks_switch);
    plot(x,nanmean(temp_data,1),'Color',these_rgbs(n,:))
    %shadedErrorBar(x,nanmean(temp_data,1),nansem(temp_data,1),'lineProps',these_rgbs(m,:));
    
    xlim([0,5])
    ylim([0,0.1])
end



%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! --- Warp -> performance correlation ---

F = default_figure([-20,0.5,10,5]); hold on;

this_data_x = warp_exp1_a; %this_data_x = warp_exp1_b;
this_data_y = correct_catch;

this_data_x_seq = [warp_exp1_a(d_info.group==7,6:18); warp_exp1_a(d_info.group==8,19:31)];
this_data_y_seq = [correct_catch(d_info.group==7,6:18); correct_catch(d_info.group==8,19:31)];
this_data_x_ctrl = [warp_exp1_a(d_info.group==8,6:18); warp_exp1_a(d_info.group==7,19:31)];
this_data_y_ctrl = [correct_catch(d_info.group==8,6:18); correct_catch(d_info.group==7,19:31)];

subplot(1,2,1); hold on;
yline(50,'k:');
plot(this_data_x(:),this_data_y(:)*100,'ko')
plot(this_data_x_ctrl(:),this_data_y_ctrl(:)*100,'o','Color',p.col.ctrl)
plot(this_data_x_seq(:),this_data_y_seq(:)*100,'o','Color',p.col.seq)
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k');
%xlim([0,0.25])
ylim([0,100])
xlabel('Exponential coefficient a (initial proportion)')
ylabel('Performance')
ytickformat('percentage')
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['\rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);

this_data_x = warp_exp1_b; %this_data_x = warp_exp1_b;
this_data_y = correct_catch;

this_data_x_seq = [warp_exp1_b(d_info.group==7,6:18); warp_exp1_b(d_info.group==8,19:31)];
this_data_y_seq = [correct_catch(d_info.group==7,6:18); correct_catch(d_info.group==8,19:31)];
this_data_x_ctrl = [warp_exp1_b(d_info.group==8,6:18); warp_exp1_b(d_info.group==7,19:31)];
this_data_y_ctrl = [correct_catch(d_info.group==8,6:18); correct_catch(d_info.group==7,19:31)];

subplot(1,2,2); hold on;
yline(50,'k:');
plot(this_data_x(:),this_data_y(:)*100,'ko')
plot(this_data_x_ctrl(:),this_data_y_ctrl(:)*100,'o','Color',p.col.ctrl)
plot(this_data_x_seq(:),this_data_y_seq(:)*100,'o','Color',p.col.seq)
% [this_corr_r,this_corr_p] = fitLine(this_data_x_ctrl(:),this_data_y_ctrl(:)*100,p.col.ctrl);
% [this_corr_r,this_corr_p] = fitLine(this_data_x_seq(:),this_data_y_seq(:)*100,p.col.seq);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k');
%xlim([0,0.25])
ylim([0,100])
xlabel('Exponential coefficient b (decay rate)')
ylabel('Performance')
ytickformat('percentage')
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05; %SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
text(SE(1),SE(2),...
    ['\rho = ',num2str(this_corr_r,2),', p = ',num2str(this_corr_p,2)],...
    'Units','data','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10);

drawnow;
savefig(F,[path.root_summary,'plots\warping\warping_seqctrl_corrWithBehaviour_',these_animals,'.fig']);
saveas(F,[path.root_summary,'plots\warping\warping_seqctrl_corrWithBehaviour_',these_animals,'.png']);





%%

this_data_x_all = warp_exp1_b(:);
this_data_y_all = correct_catch(:);

binEdges = [-2,-0.5,0,0.5];
%binEdges = [-2,-0.3,-0.15,0.5];
numBins = length(binEdges)-1;
binLabels = {'Very start-heavy, b = [-2,-0.5]','Mildly start-heavy, b = [-0.5,0]','End-heavy, b = [0,0.5]'};
%binLabels = {'Very start-heavy, b = [-2,-0.3]','Mildly start-heavy, b = [-0.3,-0.15]','End-heavy, b = [-0.15,0.5]'};

% bin data
these_bins = discretize(this_data_x_all,binEdges);
this_data_y = {this_data_y_all(these_bins==1),this_data_y_all(these_bins==2),this_data_y_all(these_bins==3)};
% this_data_x_1 = this_data_x_all(these_bins==1);
% this_data_x_1 = this_data_y_all(these_bins==1);
% this_data_x_2 = this_data_x_all(these_bins==2);
% this_data_x_2 = this_data_y_all(these_bins==2);
% this_data_x_3 = this_data_x_all(these_bins==3);
% this_data_x_3 = this_data_y_all(these_bins==3);

F = default_figure();


ranksum(this_data_y{1},this_data_y{2})
ranksum(this_data_y{2},this_data_y{3})
ranksum(this_data_y{1},this_data_y{3})

this_data = nan(nanmax(cellfun(@length, this_data_y)),length(this_data_y));
this_data(1:length(this_data_y{1}),1) = this_data_y{1};
this_data(1:length(this_data_y{2}),2) = this_data_y{2};
this_data(1:length(this_data_y{3}),3) = this_data_y{3}; 
v = violinplot(this_data*100,binLabels);
ytickformat('percentage');


%%

% ranksum(corr_correct_catch_exp1_b_switch1_seq,corr_correct_catch_exp1_b_switch1_ctrl)
% ranksum(corr_correct_catch_exp1_b_switch2_seq,corr_correct_catch_exp1_b_switch2_ctrl)
% ranksum(corr_correct_catch_exp1_b_seq,corr_correct_catch_exp1_b_ctrl)

signrank(corr_correct_catch_exp1_a_all)
signrank(corr_correct_catch_exp1_b_all)



%% !!!!!!!!!!!!!!!!!!!! --- SUMMARY FIGURE --- NumSeqCells as a function of trial block

% select data
this_data_switch = numSeqCells_switch;
these_sessions = ~d_info.excl; %& d_info.running==0;
these_blocks = [repmat(these_sessions(:,1),1,blocks(1)),repmat(these_sessions(:,2),1,blocks(2)),repmat(these_sessions(:,3),1,blocks(3)),repmat(these_sessions(:,4),1,blocks(4)),repmat(these_sessions(:,5),1,blocks(5))];
these_blocks_switch = these_blocks(:,6:18) | these_blocks(:,19:31);

nrows = 1; ncols = 1; m=0;
F = default_figure([-20,0.5,3.5,3.5]);

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
% ylim([-1.4,0.2])
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




%% First try with mixed-effect model

this_data = warp_linIdx(:,blocks(1)+1:end,4); % b%warp_exp1_b(:,blocks(1)+1:end); % b
%this_data = squeeze(seqCellsByTime_prop(1,blocks(1)+1:end,:))'; % early sequence cells
%this_data = squeeze(seqCellsByTime_prop(3,blocks(1)+1:end,:))'; % late sequence cells

%this_data = expFits(:,:,3) - expFits(:,:,1); % tiling sequence
%this_data = squeeze(seqCellsByTime_prop(1,blocks(1)+1:end,:))' - squeeze(seqCellsByTime_prop(3,blocks(1)+1:end,:))'; % early cells - late sequence cells

these_sessions = ~d_info.excl & d_info.presponsive==1;

% only include sessions with all data points
% for i=1:d_info.numAnimals
%     if nansum(these_sessions(i,:))<4
%         these_sessions(i,:)=0;
%     end
% end

these_sessions_blocks = [repmat(these_sessions(:,2),1,8),repmat(these_sessions(:,3),1,5),repmat(these_sessions(:,4),1,8),repmat(these_sessions(:,5),1,5)];
this_data(~these_sessions_blocks)=NaN;

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
temp = nan(size(this_data));
for i=1:d_info.numAnimals
%     if d_info.group(i)==2
%         temp(i,:) = zeros(1,26);
%     end
    if d_info.group(i)==7
        temp(i,:) = [1*ones(1,13),2*ones(1,13)]; % stim type 1 = seq, stim type 2 = ctrl, 
    elseif d_info.group(i)==8
        temp(i,:) = [2*ones(1,13),1*ones(1,13)];
    end
end
tbl_in.stim = nominal(temp(:));
temp = nan(size(this_data));
for i=1:d_info.numAnimals
%     if d_info.group(i)==2
%         temp(i,:) = zeros(1,26);
%     end
    if d_info.group(i)==7
        temp(i,:) = [ones(1,13),zeros(1,13)]; % seq 1 = seq, seq 0 = ctrl, 
    elseif d_info.group(i)==8
        temp(i,:) = [zeros(1,13),ones(1,13)];
    end
end
tbl_in.seq = nominal(temp(:));

tbl = table(tbl_in.animal,tbl_in.switch,tbl_in.switchday,tbl_in.block,tbl_in.block_sess,tbl_in.block_sess_norm,tbl_in.running,tbl_in.stim,tbl_in.seq,tbl_in.y,...
    'VariableNames',{'animal','switch','switchday','block','block_sess','block_sess_norm','running','stim','seq','y'});

lme = fitlme(tbl,'y ~ running*block + seq*block + (1|animal) + (1|switch)') % looks not too bad with b
%lme = fitlme(tbl,'y ~ running*block_sess + seq*block_sess + switchday + (1|animal) + (1|switch)')
%   lme = fitlme(tbl,'y ~ running*block_sess + seq*block_sess + (1|switchday) + (1|animal) + (1|switch)')

% lme = fitlme(tbl,'y ~ running*seq*block + (1|animal) + (1|switch)') % looks not too bad with b


%lme = fitlme(tbl,'y ~ stim*block + (1|animal) + (1|switch)') % looks not too bad with b

 
% lme = fitlme(tbl,'y ~ running*block_sess_norm + stim*block_sess_norm + switchday + (1|animal) + (1|switch)')
% lme = fitlme(tbl,'y ~ running*block_sess + stim*block_sess + (1|animal) + (1|switch) + (1|switch)')
% lme = fitlme(tbl,'y ~ running*block_sess_norm + stim*block_sess_norm + switchday + (1|animal) + (1|switch)')
% 

% lme = fitlme(tbl,'y ~ running*block_sess + stim*block_sess + (1|animal) + (1|switch)')
% lme = fitlme(tbl,'y ~ running*block_sess_norm + stim*block_sess_norm + (1|animal) + (1|switch)')
% 
% lme = fitlme(tbl,'y ~ running*block_sess + stim*block_sess + (1|animal) + (1|switch) + (1|switchday)')
% lme = fitlme(tbl,'y ~ running*block_sess_norm + stim*block_sess_norm  + (1|animal) + (1|switch) + (1|switchday)')


% lme = fitlme(tbl,'y ~ stim*block + (1|animal) + (1|switch)')
% lme = fitlme(tbl,'y ~ stim*block_sess + (1|animal) + (1|switch)')


%%

x = 0:0.01:5.3;
expFits = nan(d_info.numAnimals,numBlocks_switch,3);
for m=1:numBlocks_switch*2
    for idx=1:d_info.numAnimals
        if d_info.group(idx)==7 || d_info.group(idx)==8     
            this_a = warp_exp1_a(idx,5+m);
            this_b = warp_exp1_b(idx,5+m);
            this_ypred = this_a*exp(this_b*x);
            
            temp = discretize(x,binEdges);
            for n=1:numBins
                expFits(idx,m,n) = nanmean(this_ypred(temp==n));
            end
        end
    end
end






