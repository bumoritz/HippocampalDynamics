function F = snakePlots_ov20220105(sca,p,info)

%% Parameters for 8-tile Plot

ops.smoothingSd = 0; % 'sd', 0 for no smoothing
ops.split = 'odd_even'; % 'odd_even' or 'stim_catch'
ops.sort = 'ipsi_1'; % 'ipsi_1' or 'seq'
ops.selectivitySplit = NaN; % 0.2; % Nan for no splitting
ops.normalisation = 'ref_AW'; % 'ref_AW' (!) or 'ref_all' or 'ref_AW_onlyDivide' or 'ref_all_onlyDivide' or 'within_AW' or 'two-sided'
ops.illustrator = 0;

%% Parameters

data_A = sca.traces.traces_A;
data_X = sca.traces.traces_X;

if ops.smoothingSd~=0
    data_A = smoothdata(data_A,2,'gaussian',ops.smoothingSd*5);
    data_X = smoothdata(data_X,2,'gaussian',ops.smoothingSd*5);
end

if strcmp(ops.split,'odd_even')
    avgTraces_A_1 = nanmean(data_A(:,:,1:2:end),3);
    avgTraces_A_2 = nanmean(data_A(:,:,2:2:end),3);
    avgTraces_X_1 = nanmean(data_X(:,:,1:2:end),3);
    avgTraces_X_2 = nanmean(data_X(:,:,2:2:end),3);
elseif strcmp(ops.split,'stim_catch')
    avgTraces_A_1 = nanmean(data_A(:,:, find(task.var(find(task.odour1=="A"))) ),3);
    avgTraces_A_2 = nanmean(data_A(:,:, find(~task.var(find(task.odour1=="A"))) ),3);
    avgTraces_X_1 = nanmean(data_X(:,:, find(task.var(find(task.odour1=="X"))) ),3);
    avgTraces_X_2 = nanmean(data_X(:,:, find(~task.var(find(task.odour1=="X"))) ),3);
end

if strcmp(ops.sort,'ipsi_1')

    temp_A = floor((sca.shuffling_A.significant + sca.firingField_A.reliability_A + sca.firingField_A.peakInAW + sca.firingField_A.diffBl_pos)/4);
    temp_X = floor((sca.shuffling_X.significant + sca.firingField_X.reliability_X + sca.firingField_X.peakInAW + sca.firingField_X.diffBl_pos)/4);

    SELECTIVITY_THRESHOLD = 0.1;
    
    idcs_unsorted_Aseq = intersect(find(sca.passed.Aonly==1),find(sca.firingField_A.selectivity_A>SELECTIVITY_THRESHOLD));%find(temp_A==1); %find(sca.passed.Aonly);
    idcs_unsorted_Xseq = intersect(find(sca.passed.Xonly==1),find(sca.firingField_X.selectivity_X>SELECTIVITY_THRESHOLD));%find(temp_X==1);%find(sca.passed.Xonly);
    idcs_unsorted_AXseq = find(sca.passed.AandX==1);   


    if isnan(ops.selectivitySplit)
        [~,temp] = nanmax(avgTraces_A_1(idcs_unsorted_Aseq,sca.prop.analysisWindow),[],2);
        [~,temp2] = sort(temp);
        idcs_sorted_Aseq = idcs_unsorted_Aseq(temp2);
        numSel_A = NaN;
        [~,temp] = nanmax(avgTraces_X_1(idcs_unsorted_Xseq,sca.prop.analysisWindow),[],2);
        [~,temp2] = sort(temp);
        idcs_sorted_Xseq = idcs_unsorted_Xseq(temp2);
        numSel_X = NaN;
        [~,temp] = nanmax(avgTraces_A_1(idcs_unsorted_AXseq,sca.prop.analysisWindow),[],2);
        [~,temp2] = sort(temp);
        idcs_sorted_AXseq = idcs_unsorted_AXseq(temp2);
        numSel_AX = NaN;
    else
        idcs_unsorted_Aseq_sel = intersect(idcs_unsorted_Aseq,find(sca.firingField_A.selectivity_A >= ops.selectivitySplit));
        [~,temp] = nanmax(avgTraces_A_1(idcs_unsorted_Aseq_sel,sca.prop.analysisWindow),[],2);
        [~,temp2] = sort(temp);
        idcs_sorted_Aseq_sel = idcs_unsorted_Aseq_sel(temp2);
        idcs_unsorted_Aseq_unsel = intersect(idcs_unsorted_Aseq,find(sca.firingField_A.selectivity_A < ops.selectivitySplit));
        [~,temp] = nanmax(avgTraces_A_1(idcs_unsorted_Aseq_unsel,sca.prop.analysisWindow),[],2);
        [~,temp2] = sort(temp);
        idcs_sorted_Aseq_unsel = idcs_unsorted_Aseq_unsel(temp2);
        idcs_sorted_Aseq = [idcs_sorted_Aseq_sel;idcs_sorted_Aseq_unsel];
        numSel_A = length(idcs_sorted_Aseq_sel);
        idcs_unsorted_Xseq_sel = intersect(idcs_unsorted_Xseq,find(sca.firingField_X.selectivity_X >= ops.selectivitySplit));
        [~,temp] = nanmax(avgTraces_X_1(idcs_unsorted_Xseq_sel,sca.prop.analysisWindow),[],2);
        [~,temp2] = sort(temp);
        idcs_sorted_Xseq_sel = idcs_unsorted_Xseq_sel(temp2);
        idcs_unsorted_Xseq_unsel = intersect(idcs_unsorted_Xseq,find(sca.firingField_X.selectivity_X < ops.selectivitySplit));
        [~,temp] = nanmax(avgTraces_X_1(idcs_unsorted_Xseq_unsel,sca.prop.analysisWindow),[],2);
        [~,temp2] = sort(temp);
        idcs_sorted_Xseq_unsel = idcs_unsorted_Xseq_unsel(temp2);
        idcs_sorted_Xseq = [idcs_sorted_Xseq_sel;idcs_sorted_Xseq_unsel];
        numSel_X = length(idcs_sorted_Xseq_sel);
    end
elseif strcmp(ops.sort,'seq')
    if isnan(ops.selectivitySplit)
        temp = trg.idcs_targetedCells(:,trg.sequenceClusters(:,1));
        idcs_sorted_Aseq = rmmissing(temp(:));
        numSel_A = NaN;
        temp = trg.idcs_targetedCells(:,trg.sequenceClusters(:,2));
        idcs_sorted_Xseq = rmmissing(temp(:));
        numSel_X = NaN;
    else
        error('Selectivity split for seq ordering not implemented.')
    end
end

if strcmp(ops.normalisation,'two-sided')
    temp = avgTraces_A_1;
    temp2 = temp(:,1:sca.prop.frames_pre_binned);
    this_zero = nanmedian(temp2,2);
    normAvgTraces_A_1_Aseq = nan(size(avgTraces_A_1,1),size(avgTraces_A_1,2)+3);
    normAvgTraces_A_2_Aseq = nan(size(avgTraces_A_2,1),size(avgTraces_A_2,2)+3);
    normAvgTraces_X_1_Aseq = nan(size(avgTraces_X_1,1),size(avgTraces_X_1,2)+3);
    normAvgTraces_X_2_Aseq = nan(size(avgTraces_X_2,1),size(avgTraces_X_2,2)+3);
    for i=1:size(temp,1)
        if abs(nanmin(temp(i,:))-this_zero(i)) > abs(nanmax(temp(i,:))-this_zero(i))
            this_lower = nanmin(temp(i,:));
            this_upper = abs(nanmin(temp(i,:)))+2*this_zero(i);
        else
            this_lower = -nanmax(temp(i,:))+2*this_zero(i);
            this_upper = nanmax(temp(i,:));            
        end
        normAvgTraces_A_1_Aseq(i,:) = rescale([avgTraces_A_1(i,:),[this_lower,this_zero(i),this_upper]],-1,1);
        normAvgTraces_A_2_Aseq(i,:) = rescale([avgTraces_A_2(i,:),[this_lower,this_zero(i),this_upper]],-1,1);
        normAvgTraces_X_1_Aseq(i,:) = rescale([avgTraces_X_1(i,:),[this_lower,this_zero(i),this_upper]],-1,1);
        normAvgTraces_X_2_Aseq(i,:) = rescale([avgTraces_X_2(i,:),[this_lower,this_zero(i),this_upper]],-1,1);
    end
    normAvgTraces_A_1_Aseq = normAvgTraces_A_1_Aseq(:,1:end-3);
    normAvgTraces_A_2_Aseq = normAvgTraces_A_2_Aseq(:,1:end-3);
    normAvgTraces_X_1_Aseq = normAvgTraces_X_1_Aseq(:,1:end-3);
    normAvgTraces_X_2_Aseq = normAvgTraces_X_2_Aseq(:,1:end-3);    
    %%%
    temp = avgTraces_X_1;
    temp2 = temp(:,1:sca.prop.frames_pre_binned);
    this_zero = nanmedian(temp2,2);
    normAvgTraces_A_1_Xseq = nan(size(avgTraces_A_1,1),size(avgTraces_A_1,2)+3);
    normAvgTraces_A_2_Xseq = nan(size(avgTraces_A_2,1),size(avgTraces_A_2,2)+3);
    normAvgTraces_X_1_Xseq = nan(size(avgTraces_X_1,1),size(avgTraces_X_1,2)+3);
    normAvgTraces_X_2_Xseq = nan(size(avgTraces_X_2,1),size(avgTraces_X_2,2)+3);
    for i=1:size(temp,1)
        if abs(nanmin(temp(i,:))-this_zero(i)) > abs(nanmax(temp(i,:))-this_zero(i))
            this_lower = nanmin(temp(i,:));
            this_upper = abs(nanmin(temp(i,:)))+2*this_zero(i);
        else
            this_lower = -nanmax(temp(i,:))+2*this_zero(i);
            this_upper = nanmax(temp(i,:));            
        end
        normAvgTraces_A_1_Xseq(i,:) = rescale([avgTraces_A_1(i,:),[this_lower,this_zero(i),this_upper]],-1,1);
        normAvgTraces_A_2_Xseq(i,:) = rescale([avgTraces_A_2(i,:),[this_lower,this_zero(i),this_upper]],-1,1);
        normAvgTraces_X_1_Xseq(i,:) = rescale([avgTraces_X_1(i,:),[this_lower,this_zero(i),this_upper]],-1,1);
        normAvgTraces_X_2_Xseq(i,:) = rescale([avgTraces_X_2(i,:),[this_lower,this_zero(i),this_upper]],-1,1);
    end
    normAvgTraces_A_1_Xseq = normAvgTraces_A_1_Xseq(:,1:end-3);
    normAvgTraces_A_2_Xseq = normAvgTraces_A_2_Xseq(:,1:end-3);
    normAvgTraces_X_1_Xseq = normAvgTraces_X_1_Xseq(:,1:end-3);
    normAvgTraces_X_2_Xseq = normAvgTraces_X_2_Xseq(:,1:end-3);  
elseif strcmp(ops.normalisation,'none')
    normAvgTraces_A_1_Aseq = avgTraces_A_1;
    normAvgTraces_A_2_Aseq = avgTraces_A_2;
    normAvgTraces_X_1_Aseq = avgTraces_X_1;
    normAvgTraces_X_2_Aseq = avgTraces_X_2;
    normAvgTraces_A_1_Xseq = avgTraces_A_1;
    normAvgTraces_A_2_Xseq = avgTraces_A_2;
    normAvgTraces_X_1_Xseq = avgTraces_X_1;
    normAvgTraces_X_2_Xseq = avgTraces_X_2;
    normAvgTraces_A_1_AXseq = avgTraces_A_1;
    normAvgTraces_A_2_AXseq = avgTraces_A_2;
    normAvgTraces_X_1_AXseq = avgTraces_X_1;
    normAvgTraces_X_2_AXseq = avgTraces_X_2;
elseif strcmp(ops.normalisation,'ref_AW')
    normAvgTraces_A_1_Aseq = (avgTraces_A_1-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_A_2_Aseq = (avgTraces_A_2-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_X_1_Aseq = (avgTraces_X_1-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_X_2_Aseq = (avgTraces_X_2-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_A_1_Xseq = (avgTraces_A_1-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_A_2_Xseq = (avgTraces_A_2-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_X_1_Xseq = (avgTraces_X_1-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_X_2_Xseq = (avgTraces_X_2-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_A_1_AXseq = (avgTraces_A_1-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_A_2_AXseq = (avgTraces_A_2-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_X_1_AXseq = (avgTraces_X_1-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_X_2_AXseq = (avgTraces_X_2-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
elseif strcmp(ops.normalisation,'ref_all')
    normAvgTraces_A_1_Aseq = (avgTraces_A_1-nanmin(avgTraces_A_1,[],2)) ./ (nanmax(avgTraces_A_1,[],2)-nanmin(avgTraces_A_1,[],2));
    normAvgTraces_A_2_Aseq = (avgTraces_A_2-nanmin(avgTraces_A_1,[],2)) ./ (nanmax(avgTraces_A_1,[],2)-nanmin(avgTraces_A_1,[],2));
    normAvgTraces_X_1_Aseq = (avgTraces_X_1-nanmin(avgTraces_A_1,[],2)) ./ (nanmax(avgTraces_A_1,[],2)-nanmin(avgTraces_A_1,[],2));
    normAvgTraces_X_2_Aseq = (avgTraces_X_2-nanmin(avgTraces_A_1,[],2)) ./ (nanmax(avgTraces_A_1,[],2)-nanmin(avgTraces_A_1,[],2));
    normAvgTraces_A_1_Xseq = (avgTraces_A_1-nanmin(avgTraces_X_1,[],2)) ./ (nanmax(avgTraces_X_1,[],2)-nanmin(avgTraces_X_1,[],2));
    normAvgTraces_A_2_Xseq = (avgTraces_A_2-nanmin(avgTraces_X_1,[],2)) ./ (nanmax(avgTraces_X_1,[],2)-nanmin(avgTraces_X_1,[],2));
    normAvgTraces_X_1_Xseq = (avgTraces_X_1-nanmin(avgTraces_X_1,[],2)) ./ (nanmax(avgTraces_X_1,[],2)-nanmin(avgTraces_X_1,[],2));
    normAvgTraces_X_2_Xseq = (avgTraces_X_2-nanmin(avgTraces_X_1,[],2)) ./ (nanmax(avgTraces_X_1,[],2)-nanmin(avgTraces_X_1,[],2));
    normAvgTraces_A_1_AXseq = (avgTraces_A_1-nanmin(avgTraces_A_1,[],2)) ./ (nanmax(avgTraces_A_1,[],2)-nanmin(avgTraces_A_1,[],2));
    normAvgTraces_A_2_AXseq = (avgTraces_A_2-nanmin(avgTraces_A_1,[],2)) ./ (nanmax(avgTraces_A_1,[],2)-nanmin(avgTraces_A_1,[],2));
    normAvgTraces_X_1_AXseq = (avgTraces_X_1-nanmin(avgTraces_A_1,[],2)) ./ (nanmax(avgTraces_A_1,[],2)-nanmin(avgTraces_A_1,[],2));
    normAvgTraces_X_2_AXseq = (avgTraces_X_2-nanmin(avgTraces_A_1,[],2)) ./ (nanmax(avgTraces_A_1,[],2)-nanmin(avgTraces_A_1,[],2));
elseif strcmp(ops.normalisation,'ref_AW_onlyDevide')
    normAvgTraces_A_1_Aseq = avgTraces_A_1 ./ nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2);
    normAvgTraces_A_2_Aseq = avgTraces_A_2 ./ nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2);
    normAvgTraces_X_1_Aseq = avgTraces_X_1 ./ nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2);
    normAvgTraces_X_2_Aseq = avgTraces_X_2 ./ nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2);
    normAvgTraces_A_1_Xseq = avgTraces_A_1 ./ nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2);
    normAvgTraces_A_2_Xseq = avgTraces_A_2 ./ nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2);
    normAvgTraces_X_1_Xseq = avgTraces_X_1 ./ nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2);
    normAvgTraces_X_2_Xseq = avgTraces_X_2 ./ nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2);
    normAvgTraces_A_1_AXseq = avgTraces_A_1 ./ nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2);
    normAvgTraces_A_2_AXseq = avgTraces_A_2 ./ nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2);
    normAvgTraces_X_1_AXseq = avgTraces_X_1 ./ nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2);
    normAvgTraces_X_2_AXseq = avgTraces_X_2 ./ nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2);
elseif strcmp(ops.normalisation,'ref_all_onlyDevide')
    normAvgTraces_A_1_Aseq = avgTraces_A_1 ./ nanmax(avgTraces_A_1,[],2);
    normAvgTraces_A_2_Aseq = avgTraces_A_2 ./ nanmax(avgTraces_A_1,[],2);
    normAvgTraces_X_1_Aseq = avgTraces_X_1 ./ nanmax(avgTraces_A_1,[],2);
    normAvgTraces_X_2_Aseq = avgTraces_X_2 ./ nanmax(avgTraces_A_1,[],2);
    normAvgTraces_A_1_Xseq = avgTraces_A_1 ./ nanmax(avgTraces_X_1,[],2);
    normAvgTraces_A_2_Xseq = avgTraces_A_2 ./ nanmax(avgTraces_X_1,[],2);
    normAvgTraces_X_1_Xseq = avgTraces_X_1 ./ nanmax(avgTraces_X_1,[],2);
    normAvgTraces_X_2_Xseq = avgTraces_X_2 ./ nanmax(avgTraces_X_1,[],2);
    normAvgTraces_A_1_AXseq = avgTraces_A_1 ./ nanmax(avgTraces_A_1,[],2);
    normAvgTraces_A_2_AXseq = avgTraces_A_2 ./ nanmax(avgTraces_A_1,[],2);
    normAvgTraces_X_1_AXseq = avgTraces_X_1 ./ nanmax(avgTraces_A_1,[],2);
    normAvgTraces_X_2_AXseq = avgTraces_X_2 ./ nanmax(avgTraces_A_1,[],2);
elseif strcmp(ops.normalisation,'within_AW')
    normAvgTraces_A_1_Aseq = (avgTraces_A_1-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_A_2_Aseq = (avgTraces_A_2-nanmin(avgTraces_A_2(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_2(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_2(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_X_1_Aseq = (avgTraces_X_1-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_X_2_Aseq = (avgTraces_X_2-nanmin(avgTraces_X_2(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_2(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_2(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_A_1_Xseq = (avgTraces_A_1-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_A_2_Xseq = (avgTraces_A_2-nanmin(avgTraces_A_2(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_2(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_2(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_X_1_Xseq = (avgTraces_X_1-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_X_2_Xseq = (avgTraces_X_2-nanmin(avgTraces_X_2(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_2(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_2(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_A_1_AXseq = (avgTraces_A_1-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_A_2_AXseq = (avgTraces_A_2-nanmin(avgTraces_A_2(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_A_2(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_A_2(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_X_1_AXseq = (avgTraces_X_1-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_1(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_1(:,sca.prop.analysisWindow),[],2));
    normAvgTraces_X_2_AXseq = (avgTraces_X_2-nanmin(avgTraces_X_2(:,sca.prop.analysisWindow),[],2)) ./ (nanmax(avgTraces_X_2(:,sca.prop.analysisWindow),[],2)-nanmin(avgTraces_X_2(:,sca.prop.analysisWindow),[],2));
end


%% 8 tile

default_figure();
close;
F = figure;

plt = struct(); plt.illustrator = false; 
if strcmp(ops.normalisation,'two-sided')
    plt.clim = [-1,1]; plt.colormap = redblue;
else
    plt.clim = [0,1]; plt.colormap = parula;
end
subplot(2,4,1)
heatMap_task_ov20220105(normAvgTraces_A_1_Aseq(idcs_sorted_Aseq,:),numSel_A,sca,p,info,plt);
subplot(2,4,2)
heatMap_task_ov20220105(normAvgTraces_A_2_Aseq(idcs_sorted_Aseq,:),numSel_A,sca,p,info,plt);
subplot(2,4,3)
heatMap_task_ov20220105(normAvgTraces_X_1_Aseq(idcs_sorted_Aseq,:),numSel_A,sca,p,info,plt);
subplot(2,4,4)
heatMap_task_ov20220105(normAvgTraces_X_2_Aseq(idcs_sorted_Aseq,:),numSel_A,sca,p,info,plt);
subplot(2,4,5)
heatMap_task_ov20220105(normAvgTraces_A_1_Xseq(idcs_sorted_Xseq,:),numSel_X,sca,p,info,plt);
subplot(2,4,6)
heatMap_task_ov20220105(normAvgTraces_A_2_Xseq(idcs_sorted_Xseq,:),numSel_X,sca,p,info,plt);
subplot(2,4,7)
heatMap_task_ov20220105(normAvgTraces_X_1_Xseq(idcs_sorted_Xseq,:),numSel_X,sca,p,info,plt);
subplot(2,4,8)
heatMap_task_ov20220105(normAvgTraces_X_2_Xseq(idcs_sorted_Xseq,:),numSel_X,sca,p,info,plt);

if ~ops.illustrator
    for i=1:8
        subplot(2,4,i)
        if strcmp(ops.sort,'ipsi_1')
            if i==1
                ylabel('Sequence A cell')
            end
            if i==5
                ylabel('Sequence X cell')
            end
        elseif strcmp(ops.sort,'seq')
            if i==1
                ylabel('Seq Stim A target')
            end
            if i==5
                ylabel('Seq Stim X target')
            end
        end
        if strcmp(ops.split,'odd_even')
            plotTitles = ["Odour A trials - odd","Odour A trials - even","Odour X trials - odd","Odour X trials - even","Odour A trials - odd","Odour A trials - even","Odour X trials - odd","Odour X trials - even"];
            title(plotTitles(i))
        elseif strcmp(ops.split,'stim_catch')
            plotTitles = ["Odour A trials - stim","Odour A trials - catch","Odour X trials - stim","Odour X trials - catch","Odour A trials - stim","Odour A trials - catch","Odour X trials - stim","Odour X trials - catch"];
            title(plotTitles(i))
        end
    end
end

end