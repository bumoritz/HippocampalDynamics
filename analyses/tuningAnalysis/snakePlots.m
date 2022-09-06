function F = snakePlots(this_tracesStruct,these_idcs_unsorted,these_idcs_sorted,these_refs_sort,these_refs_norm,p,info,plt)
%this_tracesStruct = traces_all; these_idcs_unsorted = {find(passed_all.AW.Aonly_sel==1),find(passed_all.AW.Xonly_sel==1)}; these_idcs_sorted = []; these_refs_sort = [1;4]; these_refs_norm = [];
%this_tracesStruct = traces_60t; these_idcs_unsorted = {find(passed_all.AW.Aonly_sel==1)}; these_idcs_sorted = []; these_refs_sort = [1;NaN]; these_refs_norm = [];

if ~exist('plt','var')
    plt = struct();
end
if ~isfield(plt,'trialTypes')
    plt.trialTypes          = "A_X";            % "A","X","A_X"
end
if ~isfield(plt,'sequences')
    plt.sequences           = plt.trialTypes;
end
if ~isfield(plt,'split')
    plt.split               = 'odd_even';       % 'none','odd_even','corr_incorr','corrO_corrM_incorr','corrTrain_corrTest_incorrTest','stim_catch','stimO_stimM_catch'
end
if ~isfield(plt,'split3_numTrials')
    plt.split3_numTrials    = 'same';           % 'max','same'
end
if ~isfield(plt,'illustrator')
    plt.illustrator         = false;
end
if ~isfield(plt,'sortingWindow')
    plt.sortingWindow        = 'AW';            % 'AW','all' ['AW' is p.general.bins_analysisWindow]
end
if ~isfield(plt,'normalisationType')
    plt.normalisationType    = 'one-sided';     % 'none','one-sided','two-sided' ['two-sided' needs p.general.bins_baselineWindow]
end
if ~isfield(plt,'normalisationWindow')
    plt.normalisationWindow  = 'AW';            % 'AW','all' ['AW' is p.general.bins_analysisWindow]
end
if ~isfield(plt,'normalisationRef')
    plt.normalisationRef     = 'rowwise_all';	% 'rowwise_all','rowwise_ref','plotwise','uniform'
end
if ~isfield(plt,'normalisationMetric')
    plt.normalisationMetric  = 'standard';      % 'standard' or 'onlyDivide'
end

if strcmp(plt.normalisationType,'two-sided')
    plt.clim = [-1,1]; plt.colormap = redblue;
elseif strcmp(plt.normalisationType,'one-sided')
    plt.clim = [0,1]; plt.colormap = parula;
elseif strcmp(plt.normalisationType,'none')
    plt.colormap = parula;
end


%% Prepare plot layout

numBlocks = size(this_tracesStruct,2);
if ~isempty(these_idcs_sorted)
    numSequences = size(these_idcs_sorted,2);
else
    numSequences = size(these_idcs_unsorted,2);
end
numTrialTypes = length(split(plt.trialTypes,'_'));
numSplits = length(split(plt.split,'_'));

if numBlocks==1   
    nrows = numSequences;
    ncols = numTrialTypes*numSplits;
else
    if numSequences~=1
        error('Blockwise plots with multiple sequences not implemented.')
    end
    if numSplits~=1 && numTrialTypes~=1
        error('Blockwise plots with both splits and multiple trial types not implemented.')
    end
    if numSplits~=1
        nrows = numSplits;
    elseif numTrialTypes~=1
        nrows = numTrialTypes;
    end
    ncols = numBlocks;
end

% Sanity checks
if strcmp(plt.normalisationType,'two-sided') && ~strcmp(plt.normalisationWindow,'all')
    error('plt.normalisationWindow needs to be all for two-sided normalisation')
end


%% Calculate average traces

these_avgTraces = cell(nrows,ncols);

numTrials = [];
if numBlocks==1
    for i=1:numSequences
        temp = split(plt.trialTypes,'_');
        for j=1:numTrialTypes
            if temp{j}=="A" && strcmp(plt.split,'none')
                these_avgTraces{i,j} = nanmean(this_tracesStruct.A,3);
                numTrials = [numTrials,size(this_tracesStruct.A,3)];
            elseif temp{j}=="X" && strcmp(plt.split,'none')
                these_avgTraces{i,j} = nanmean(this_tracesStruct.X,3);
                numTrials = [numTrials,size(this_tracesStruct.X,3)];
            end
            if temp{j}=="A" && strcmp(plt.split,'odd_even')
                these_avgTraces{i,j*2-1} = nanmean(this_tracesStruct.A(:,:,1:2:end),3);
                these_avgTraces{i,j*2} = nanmean(this_tracesStruct.A(:,:,2:2:end),3);
                numTrials = [numTrials,size(this_tracesStruct.A(:,:,1:2:end),3),size(this_tracesStruct.A(:,:,2:2:end),3)];
            elseif temp{j}=="X" && strcmp(plt.split,'odd_even')    
                these_avgTraces{i,j*2-1} = nanmean(this_tracesStruct.X(:,:,1:2:end),3);
                these_avgTraces{i,j*2} = nanmean(this_tracesStruct.X(:,:,2:2:end),3);
                numTrials = [numTrials,size(this_tracesStruct.X(:,:,1:2:end),3),size(this_tracesStruct.X(:,:,2:2:end),3)];
            end
            if temp{j}=="A" && strcmp(plt.split,'corr_incorr')
                these_avgTraces{i,j*2-1} = nanmean(this_tracesStruct.A_correct,3);
                these_avgTraces{i,j*2} = nanmean(this_tracesStruct.A_incorrect,3);
                numTrials = [numTrials,size(this_tracesStruct.A_correct,3),size(this_tracesStruct.A_incorrect,3)];
            elseif temp{j}=="X" && strcmp(plt.split,'corr_incorr')
                these_avgTraces{i,j*2-1} = nanmean(this_tracesStruct.X_correct,3);
                these_avgTraces{i,j*2} = nanmean(this_tracesStruct.X_incorrect,3);
                numTrials = [numTrials,size(this_tracesStruct.X_correct,3),size(this_tracesStruct.X_incorrect,3)];
            end
            if strcmp(plt.split,'corrO_corrM_incorr') && strcmp(plt.split3_numTrials,'max')
                numTrials_A_corr = size(this_tracesStruct.A_correct,3);
                numTrials_A_incorr = size(this_tracesStruct.A_incorrect,3);
                numTrials_X_corr = size(this_tracesStruct.X_correct,3);
                numTrials_X_incorr = size(this_tracesStruct.X_incorrect,3);
            elseif  strcmp(plt.split,'corrO_corrM_incorr') && strcmp(plt.split3_numTrials,'same') 
                numTrials_A_corr = size(this_tracesStruct.A_correct,3);
                numTrials_X_corr = size(this_tracesStruct.X_correct,3);
                numTrials_A_incorr = nanmin([size(this_tracesStruct.A_incorrect,3),size(this_tracesStruct.X_incorrect,3)]);
                numTrials_X_incorr = nanmin([size(this_tracesStruct.A_incorrect,3),size(this_tracesStruct.X_incorrect,3)]);
            end
            if temp{j}=="A" && strcmp(plt.split,'corrO_corrM_incorr')
                rng(p.general.rgnSeed);
                selection_A_corrm = datasample(1:numTrials_A_corr,numTrials_A_incorr,'Replace',false);
                selection_A_corro = setdiff(1:numTrials_A_corr,selection_A_corrm);
                these_avgTraces{i,j*3-2} = nanmean(this_tracesStruct.A_correct(:,:,selection_A_corro),3);                
                these_avgTraces{i,j*3-1} = nanmean(this_tracesStruct.A_correct(:,:,selection_A_corrm),3);
                these_avgTraces{i,j*3} = nanmean(this_tracesStruct.A_incorrect,3);
                numTrials = [numTrials,numTrials_A_corr-numTrials_A_incorr,numTrials_A_incorr,numTrials_A_incorr];
            elseif temp{j}=="X" && strcmp(plt.split,'corrO_corrM_incorr')
                rng(p.general.rgnSeed);
                selection_X_corrm = datasample(1:numTrials_X_corr,numTrials_X_incorr,'Replace',false);
                selection_X_corro = setdiff(1:numTrials_X_corr,selection_X_corrm);
                these_avgTraces{i,j*3-2} = nanmean(this_tracesStruct.X_correct(:,:,selection_X_corro),3);                
                these_avgTraces{i,j*3-1} = nanmean(this_tracesStruct.X_correct(:,:,selection_X_corrm),3);
                these_avgTraces{i,j*3} = nanmean(this_tracesStruct.X_incorrect,3);
                numTrials = [numTrials,numTrials_X_corr-numTrials_X_incorr,numTrials_X_incorr,numTrials_X_incorr];
            end
            if temp{j}=="A" && strcmp(plt.split,'corrTrain_corrTest_incorrTest')
                rng(p.general.rgnSeed);
                these_avgTraces{i,j*3-2} = nanmean(this_tracesStruct.A_correct_train,3);                
                these_avgTraces{i,j*3-1} = nanmean(this_tracesStruct.A_correct_test,3);
                these_avgTraces{i,j*3} = nanmean(this_tracesStruct.A_incorrect,3);
                numTrials = [numTrials,size(this_tracesStruct.A_correct_train,3),size(this_tracesStruct.A_correct_test,3),size(this_tracesStruct.A_incorrect,3)];
            elseif temp{j}=="X" && strcmp(plt.split,'corrTrain_corrTest_incorrTest')
                rng(p.general.rgnSeed);
                these_avgTraces{i,j*3-2} = nanmean(this_tracesStruct.X_correct_train,3);                
                these_avgTraces{i,j*3-1} = nanmean(this_tracesStruct.X_correct_test,3);
                these_avgTraces{i,j*3} = nanmean(this_tracesStruct.X_incorrect,3);
                numTrials = [numTrials,size(this_tracesStruct.X_correct_train,3),size(this_tracesStruct.X_correct_test,3),size(this_tracesStruct.X_incorrect,3)];
            end
            if temp{j}=="A" && strcmp(plt.split,'stim_catch')
                these_avgTraces{i,j*2-1} = nanmean(this_tracesStruct.A_var1,3);
                these_avgTraces{i,j*2} = nanmean(this_tracesStruct.A_var0,3);
                numTrials = [numTrials,size(this_tracesStruct.A_var1,3),size(this_tracesStruct.A_var0,3)];
            elseif temp{j}=="X" && strcmp(plt.split,'stim_catch')
                these_avgTraces{i,j*2-1} = nanmean(this_tracesStruct.X_var1,3);
                these_avgTraces{i,j*2} = nanmean(this_tracesStruct.X_var0,3);
                numTrials = [numTrials,size(this_tracesStruct.X_var1,3),size(this_tracesStruct.X_var0,3)];
            end
            if strcmp(plt.split,'stimO_stimM_catch') && strcmp(plt.split3_numTrials,'max')
                numTrials_A_stim = size(this_tracesStruct.A_var1,3);
                numTrials_A_catch = size(this_tracesStruct.A_var0,3);
                numTrials_A_stim = size(this_tracesStruct.X_var1,3);
                numTrials_X_catch = size(this_tracesStruct.X_var0,3);
            elseif  strcmp(plt.split,'stimO_stimM_catch') && strcmp(plt.split3_numTrials,'same') 
                numTrials_A_stim = size(this_tracesStruct.A_var1,3);
                numTrials_X_stim = size(this_tracesStruct.X_var1,3);
                numTrials_A_catch = nanmin([size(this_tracesStruct.A_var0,3),size(this_tracesStruct.X_var0,3)]);
                numTrials_X_catch = nanmin([size(this_tracesStruct.A_var0,3),size(this_tracesStruct.X_var0,3)]);
            end
            if temp{j}=="A" && strcmp(plt.split,'stimO_stimM_catch')
                rng(p.general.rgnSeed);
                selection_A_stimm = datasample(1:numTrials_A_stim,numTrials_A_catch,'Replace',false);
                selection_A_stimo = setdiff(1:numTrials_A_stim,selection_A_stimm);
                these_avgTraces{i,j*3-2} = nanmean(this_tracesStruct.A_var1(:,:,selection_A_stimo),3);                
                these_avgTraces{i,j*3-1} = nanmean(this_tracesStruct.A_var1(:,:,selection_A_stimm),3);
                these_avgTraces{i,j*3} = nanmean(this_tracesStruct.A_var0,3);
                numTrials = [numTrials,numTrials_A_stim-numTrials_A_catch,numTrials_A_catch,numTrials_A_catch];
            elseif temp{j}=="X" && strcmp(plt.split,'stimO_stimM_catch')
                rng(p.general.rgnSeed);
                selection_X_stimm = datasample(1:numTrials_A_stim,numTrials_X_catch,'Replace',false);
                selection_X_stimo = setdiff(1:numTrials_A_stim,selection_X_stimm);
                these_avgTraces{i,j*3-2} = nanmean(this_tracesStruct.X_var1(:,:,selection_X_stimo),3);                
                these_avgTraces{i,j*3-1} = nanmean(this_tracesStruct.X_var1(:,:,selection_X_stimm),3);
                these_avgTraces{i,j*3} = nanmean(this_tracesStruct.X_var0,3);
                numTrials = [numTrials,numTrials_A_stim-numTrials_X_catch,numTrials_X_catch,numTrials_X_catch];
            end
        end
    end
else
    for j=1:numBlocks
        temp = split(plt.trialTypes,'_');
        if numTrialTypes~=1
            for i=1:numTrialTypes
                if temp{i}=="A" && strcmp(plt.split,'none')
                    these_avgTraces{i,j} = nanmean(this_tracesStruct{j}.A,3);
                    numTrials = [numTrials,size(this_tracesStruct{j}.A,3)];
                elseif temp{i}=="X" && strcmp(plt.split,'none')
                    these_avgTraces{i,j} = nanmean(this_tracesStruct{j}.X,3);
                    numTrials = [numTrials,size(this_tracesStruct{j}.X,3)];
                end
            end
        elseif numSplits~=1
            if temp{1}=="A" && strcmp(plt.split,'stim_catch')
                these_avgTraces{1,j} = nanmean(this_tracesStruct{j}.A_var1,3);
                these_avgTraces{2,j} = nanmean(this_tracesStruct{j}.A_var0,3);
                numTrials = [numTrials,size(this_tracesStruct{j}.A_var1,3),size(this_tracesStruct{j}.A_var0,3)];                    
            elseif temp{1}=="X" && strcmp(plt.split,'stim_catch')
                these_avgTraces{1,j} = nanmean(this_tracesStruct{j}.X_var1,3);
                these_avgTraces{2,j} = nanmean(this_tracesStruct{j}.X_var0,3);
                numTrials = [numTrials,size(this_tracesStruct{j}.X_var1,3),size(this_tracesStruct{j}.X_var0,3)];
            end
        end
    end
end

            
%% Sort traces

if isempty(these_idcs_sorted)
    
    % define normalisation window
    if strcmp(plt.sortingWindow,'all')
        this_sortingWindow = 1:size(these_avgTraces{1,1},2);
    elseif strcmp(plt.sortingWindow,'AW')
        this_sortingWindow = p.general.bins_analysisWindow;
    end
    
    these_idcs_sorted = {};
    for i=1:length(these_refs_sort)
        if ~isnan(these_refs_sort(i))
            try %if i<=length(these_idcs_unsorted)
                [~,temp] = nanmax(these_avgTraces{i,these_refs_sort(i)}(these_idcs_unsorted{i},this_sortingWindow),[],2);
                [~,temp2] = sort(temp);
                these_idcs_sorted{i} = these_idcs_unsorted{i}(temp2);
            catch % else
                [~,temp] = nanmax(these_avgTraces{i,these_refs_sort(i)}(these_idcs_unsorted{1},this_sortingWindow),[],2);
                [~,temp2] = sort(temp);
                these_idcs_sorted{i} = these_idcs_unsorted{1}(temp2);
            end
        end
    end
    for i=1:length(these_refs_sort)
        if i==1 && isnan(these_refs_sort(i))
            these_idcs_sorted{i} = these_idcs_sorted{2};
        elseif i==2 && isnan(these_refs_sort(i))
            these_idcs_sorted{i} = these_idcs_sorted{1};
        elseif isnan(these_refs_sort(i))
            these_idcs_sorted{i} = these_idcs_sorted{1};
        end
    end
elseif length(these_idcs_sorted)<nrows
    for i=1:nrows
        these_idcs_sorted{i} = these_idcs_sorted{1};
    end
end


%% Normalise traces

these_normAvgTraces = cell(nrows,ncols);

if strcmp(plt.normalisationType,'none')
    for r=1:nrows
        for c=1:ncols
            these_normAvgTraces{r,c} = these_avgTraces{r,c};
        end
    end
else
    
    % define normalisation window
    if strcmp(plt.normalisationWindow,'all')
        this_normalisationWindow = 1:size(these_avgTraces{1,1},2);
    elseif strcmp(plt.normalisationWindow,'AW')
        this_normalisationWindow = p.general.bins_analysisWindow;
    end
    
    % get normalisation reference data
    if strcmp(plt.normalisationRef,'uniform')
        this_normalisationData = [];
        this_normalisationData_bl = [];
        for r=1:nrows
            for c=1:ncols
                this_normalisationData = [this_normalisationData,these_avgTraces{r,c}(:,this_normalisationWindow)];
                this_normalisationData_bl = [this_normalisationData_bl,these_avgTraces{r,c}(:,p.general.bins_baselineWindow)];
            end
        end
    end
    for r=1:nrows
        if strcmp(plt.normalisationRef,'rowwise_all')
            this_normalisationData = [];
            this_normalisationData_bl = [];
            for c=1:ncols
                this_normalisationData = [this_normalisationData,these_avgTraces{r,c}(:,this_normalisationWindow)];
                this_normalisationData_bl = [this_normalisationData_bl,these_avgTraces{r,c}(:,p.general.bins_baselineWindow)];
            end
        elseif strcmp(plt.normalisationRef,'rowwise_ref')
            this_normalisationData = these_avgTraces{r,these_refs_norm(r)}(:,this_normalisationWindow);
            this_normalisationData_bl = these_avgTraces{r,these_refs_norm(r)}(:,p.general.bins_baselineWindow);
        end
        for c=1:ncols
            if strcmp(plt.normalisationRef,'plotwise')
                this_normalisationData = these_avgTraces{r,c}(:,this_normalisationWindow);
                this_normalisationData_bl = these_avgTraces{r,c}(:,p.general.bins_baselineWindow);
            end            
            
            % get lower and upper normalisation values
            if strcmp(plt.normalisationType,'one-sided')
                this_lower = nanmin(this_normalisationData,[],2);
                this_upper = nanmax(this_normalisationData,[],2);
            elseif strcmp(plt.normalisationType,'two-sided')
                this_zero = nanmedian(this_normalisationData_bl(),2);
                this_lower = nan(size(this_normalisationData,1),1);
                this_upper = nan(size(this_normalisationData,1),1);
                for i=1:size(this_normalisationData,1)
                    if abs(nanmin(this_normalisationData(i,:))-this_zero(i)) > abs(nanmax(this_normalisationData(i,:))-this_zero(i))
                        this_lower(i) = nanmin(this_normalisationData(i,:));
                        this_upper(i) = abs(nanmin(this_normalisationData(i,:)))+2*this_zero(i);
                    else
                        this_lower(i) = -nanmax(this_normalisationData(i,:))+2*this_zero(i);
                        this_upper(i) = nanmax(this_normalisationData(i,:));            
                    end
                end
            end
            
            % do normalisation
            if strcmp(plt.normalisationType,'one-sided')
                if strcmp(plt.normalisationMetric,'standard')
                    these_normAvgTraces{r,c} = (these_avgTraces{r,c}-this_lower) ./ (this_upper-this_lower);
                elseif strcmp(plt.normalisationMetric,'onlyDivide')
                    these_normAvgTraces{r,c} = these_avgTraces{r,c} ./ this_upper;
                end
            elseif strcmp(plt.normalisationType,'two-sided')
                these_normAvgTraces{r,c} = nan(size(these_avgTraces{r,c},1),size(these_avgTraces{r,c},2)+3);
                for i=1:size(these_avgTraces{r,c},1)
                    these_normAvgTraces{r,c}(i,:) = rescale([these_avgTraces{r,c}(i,:),[this_lower(i),this_zero(i),this_upper(i)]],-1,1);
                end
                these_normAvgTraces{r,c} = these_normAvgTraces{r,c}(:,1:end-3);
            end
        end
    end
end


%% Make plots

F = default_figure([20,0.5,20,9.9]);

for r=1:nrows
    for c=1:ncols
        subplot(nrows,ncols,(r-1)*ncols+c)
        heatMap_task(these_normAvgTraces{r,c}(these_idcs_sorted{r},:),NaN,[],p,info,plt);
    end
end


%% Label plot

if ~plt.illustrator
    n = 0;
    temp = split(plt.trialTypes,'_');
    temp1 = split(plt.sequences,'_');
    temp2 = split(plt.split,'_');
    for r=1:nrows
        for c=1:ncols
            n = n+1;
            subplot(nrows,ncols,(r-1)*ncols+c)
            
            if numBlocks==1
                ylabel(['Sequence ',temp{r},' cell'])
                if strcmp(plt.split,'none')
                    title([temp{c},' trials (',num2str(numTrials(c)),'t)'])
                else
                    title([temp{ceil(c/numSplits)},' trials, ',temp2{n},' (',num2str(numTrials(c)),'t)'])
                end
                if mod(c,numSplits)==0
                    n = 0;
                end
                if ~isempty(these_refs_sort)
                    if these_refs_sort(r)==c
                        set(gca,'LineWidth',3);
                    end
                end
            else
                ylabel(['Sequence ',temp1{1},' cell'])
                if numSplits~=1
                    title([temp{1},' tr., ',temp2{r},', bl. ',num2str(c)])
                elseif numTrialTypes~=1
                    title([temp{r},' tr., bl. ',num2str(c)])
                end
                if mod(c,numSplits)==0
                    n = 0;
                end
                if ~isempty(these_refs_sort)
                    if these_refs_sort(r)==c
                        set(gca,'LineWidth',3);
                    end
                end
            end
        end
    end
end


%% Return

% A = these_normAvgTraces{1,2}(these_idcs_sorted{1},:);
% [x,y] = nanmax(A,[],2);
% 
% %subplot(nrows,ncols,(1-1)*ncols+2)
% subplot(nrows,ncols,(2-1)*ncols+2)
% 
% hold on
% scatter(y,1:length(y),'k.')

drawnow;
%end



