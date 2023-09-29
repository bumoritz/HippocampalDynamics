%% Fig2_SequenceAnalysis_IN

% import data using Summary_Master with ops.do_sequenceCellSummary = true;
% d_info = selectDataset(d_info,'-g02345678-d1-e01-r01-p01-l01-i00',sheet,path,ops);

save_root_fig = [path.root_summary,'figures\Fig2_fig\'];
save_root_png = [path.root_summary,'figures\Fig2_png\'];
save_root_pdf = [path.root_summary,'figures\Fig2_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig2_txt\'];

labels = categorical({'Non-runners','Runners'});
labels = reordercats(labels,{'Non-runners','Runners'});

labelsCellTypes = categorical({'earlyNR','earlyR','lateNR','lateR'});
labelsCellTypes = reordercats(labelsCellTypes,{'earlyNR','earlyR','lateNR','lateR'});

tng_struct = 'tng_all';
warp_struct = 'warp_all';


%% Get performance

correct = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~isempty(d{i})
            correct(i) = d{i}.perf.general.correct;
    end
end


%% Get cell structs

% idcs.iscells
temp = extractVariable(d,[tng_struct,'.prop.iscell'],'cell','all');
idcs.iscells = cell(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    try
        idcs.iscells{i} = find(temp{i}==1);
    catch
    end
end

% idcs.passed
these_conditions = {'Aonly','Xonly','pref'}; %{'A','X','Aonly','Xonly','AandX','AorX'};
for k=1:length(these_conditions)-1
    temp = extractVariable(d,[tng_struct,'.passed.AW.',these_conditions{k}],'cell','all');
    idcs.passed.(these_conditions{k}) = cell(d_info.numAnimals,1);
    for i=1:d_info.numAnimals
        try
            idcs.passed.(these_conditions{k}){i} = find(temp{i}==1);
        catch
        end
    end
end
idcs.passed.pref = cell(d_info.numAnimals,1);
for i=1:d_info.numAnimals
	try
        idcs.passed.pref{i} = [idcs.passed.('Aonly'){i};idcs.passed.('Xonly'){i}];
    catch
    end
end

% peakTime
temp_A = extractVariable(d,[tng_struct,'.firingField.A_AW.peakLocation_s'],'cell','all');
temp_X = extractVariable(d,[tng_struct,'.firingField.X_AW.peakLocation_s'],'cell','all');
peakTime.Atrials_Aonly = cell(d_info.numAnimals,1);
peakTime.Xtrials_Xonly = cell(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        peakTime.Atrials_Aonly{i} = temp_A{i}(idcs.passed.Aonly{i});
        peakTime.Xtrials_Xonly{i} = temp_X{i}(idcs.passed.Xonly{i});
    end
end
peakTime.pref = cell(d_info.numAnimals,1);
for i=1:d_info.numAnimals
	try
        peakTime.pref{i} = [peakTime.('Atrials_Aonly'){i};peakTime.('Xtrials_Xonly'){i}];
    catch
    end
end

% comTime
temp_A = extractVariable(d,[tng_struct,'.firingField.A_AW.comLocation_s'],'cell','all');
temp_X = extractVariable(d,[tng_struct,'.firingField.X_AW.comLocation_s'],'cell','all');
comTime.Atrials_Aonly = cell(d_info.numAnimals,1);
comTime.Xtrials_Xonly = cell(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        comTime.Atrials_Aonly{i} = temp_A{i}(idcs.passed.Aonly{i});
        comTime.Xtrials_Xonly{i} = temp_X{i}(idcs.passed.Xonly{i});
    end
end
comTime.pref = cell(d_info.numAnimals,1);
for i=1:d_info.numAnimals
	try
        comTime.pref{i} = [comTime.('Atrials_Aonly'){i};comTime.('Xtrials_Xonly'){i}];
    catch
    end
end


%% Get data arrays - number of cells

% numCells
numCells = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        numCells(i) = length(idcs.iscells{i});
    end
end

% numSeqCells, propSeqCells
these_fields = fields(idcs.passed);
for k=1:length(these_fields)
    numSeqCells.(these_fields{k}) = nan(d_info.numAnimals,1);
    propSeqCells.(these_fields{k}) = nan(d_info.numAnimals,1);
    for i=1:d_info.numAnimals
        if ~isempty(idcs.iscells{i})
            numSeqCells.(these_fields{k})(i) = length(idcs.passed.(these_fields{k}){i});
            propSeqCells.(these_fields{k})(i) = numSeqCells.(these_fields{k})(i) / numCells(i);
            if strcmp(these_fields{k},'pref')
                if numSeqCells.(these_fields{k})(i)~= d{i,1}.(warp_struct).input.numSequenceCells
                    warning('Number of sequence cells not consistent between structs')
                end
            end
        end
    end
end

% numByPeakTime, propByPeakTime_perIscell, propByPeakTime_perSeqCells
binEdges = [0:0.2:5];
these_fields = fields(peakTime);
for k=1:length(these_fields)
    numByPeakTime.(these_fields{k}) = nan(d_info.numAnimals,length(binEdges)-1);
    propByPeakTime_perIscell.(these_fields{k}) = nan(d_info.numAnimals,length(binEdges)-1);
    propByPeakTime_perSeqCells.(these_fields{k}) = nan(d_info.numAnimals,length(binEdges)-1);
    for i=1:d_info.numAnimals
        if ~isempty(idcs.iscells{i})
            temp = discretize(peakTime.(these_fields{k}){i},binEdges);
            for n=1:length(binEdges)-1
                numByPeakTime.(these_fields{k})(i,n) = nansum(temp==n);
                propByPeakTime_perIscell.(these_fields{k})(i,n) = nansum(temp==n) / length(idcs.iscells{i});
                propByPeakTime_perSeqCells.(these_fields{k})(i,n) = nansum(temp==n) / length(idcs.passed.(these_conditions{k}){i});
            end
        end
    end
end

% numByComTime, propByComTime_perIscell, propByComTime_perSeqCells
binEdges = [0:0.2:5];
these_fields = fields(comTime);
for k=1:length(these_fields)
    numByComTime.(these_fields{k}) = nan(d_info.numAnimals,length(binEdges)-1);
    propByComTime_perIscell.(these_fields{k}) = nan(d_info.numAnimals,length(binEdges)-1);
    propByComTime_perSeqCells.(these_fields{k}) = nan(d_info.numAnimals,length(binEdges)-1);
    for i=1:d_info.numAnimals
        if ~isempty(idcs.iscells{i})
            temp = discretize(comTime.(these_fields{k}){i},binEdges);
            for n=1:length(binEdges)-1
                numByComTime.(these_fields{k})(i,n) = nansum(temp==n);
                propByComTime_perIscell.(these_fields{k})(i,n) = nansum(temp==n) / length(idcs.iscells{i});
                propByComTime_perSeqCells.(these_fields{k})(i,n) = nansum(temp==n) / length(idcs.passed.(these_conditions{k}){i});
            end
        end
    end
end

% numByCellType, propByCellType_perIscell, propByCellType_perSeqCells
cellTypeEdges = [0:5/3:5]; %[0,1,5.3];
these_fields = fields(peakTime);
for k=1:length(these_fields)
    numByCellType.(these_fields{k}) = nan(d_info.numAnimals,length(cellTypeEdges)-1);
    propByCellType_perIscell.(these_fields{k}) = nan(d_info.numAnimals,length(cellTypeEdges)-1);
    propByCellType_perSeqCells.(these_fields{k}) = nan(d_info.numAnimals,length(cellTypeEdges)-1);
    for i=1:d_info.numAnimals
        if ~isempty(idcs.iscells{i})
            temp = discretize(peakTime.(these_fields{k}){i},cellTypeEdges);
            for n=1:length(cellTypeEdges)-1
                numByCellType.(these_fields{k})(i,n) = nansum(temp==n);
                propByCellType_perIscell.(these_fields{k})(i,n) = nansum(temp==n) / length(idcs.iscells{i});
                propByCellType_perSeqCells.(these_fields{k})(i,n) = nansum(temp==n) / length(idcs.passed.(these_conditions{k}){i});
            end
        end
    end
end


%% Get data arrays - place cell distribution

% numByPeakDistance, propByPeakDistance_perIscell, propByPeakDistance_perSeqCells
distanceBinEdges = [0:50:600];
numByPeakDistance = nan(d_info.numAnimals,length(distanceBinEdges)-1);
propByPeakDistance_perIscell = nan(d_info.numAnimals,length(distanceBinEdges)-1);
propByPeakDistance_perSeqCells = nan(d_info.numAnimals,length(distanceBinEdges)-1);
for i=1:d_info.numAnimals
    try
        temp = discretize(d{i}.pca_all.x_peakBin_cm,distanceBinEdges);
        for n=1:length(distanceBinEdges)-1
            numByPeakDistance(i,n) = nansum(temp==n);
            propByPeakDistance_perIscell(i,n) = nansum(temp==n) / length(idcs.iscells{i});
            propByPeakDistance_perSeqCells(i,n) = nansum(temp==n) / length(idcs.passed.pref{i});
        end
    catch
    end
end


%% Get data arrays - amplitude

% peakAmplitude, meanAmplitude, meanAmplitudeActive
peakAmplitude.Atrials_Aonly = nan(d_info.numAnimals,1);
peakAmplitude.Xtrials_Xonly = nan(d_info.numAnimals,1);
peakAmplitude.pref = nan(d_info.numAnimals,1);
meanAmplitude.Atrials_Aonly = nan(d_info.numAnimals,1);
meanAmplitude.Xtrials_Xonly = nan(d_info.numAnimals,1);
meanAmplitude.pref = nan(d_info.numAnimals,1);
meanAmplitudeActive.Atrials_Aonly = nan(d_info.numAnimals,1);
meanAmplitudeActive.Xtrials_Xonly = nan(d_info.numAnimals,1);
meanAmplitudeActive.pref = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        peakAmplitude.Atrials_Aonly(i) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.peakAmplitude_blSub(idcs.passed.Aonly{i}));
        peakAmplitude.Xtrials_Xonly(i) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.peakAmplitude_blSub(idcs.passed.Xonly{i}));
        peakAmplitude.pref(i) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.peakAmplitude_blSub(idcs.passed.Aonly{i});d{i,1}.(tng_struct).firingField.X_AW.peakAmplitude_blSub(idcs.passed.Xonly{i})]);
        meanAmplitude.Atrials_Aonly(i) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_ipsi(idcs.passed.Aonly{i}));
        meanAmplitude.Xtrials_Xonly(i) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_ipsi(idcs.passed.Xonly{i}));
        meanAmplitude.pref(i) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_ipsi(idcs.passed.Aonly{i});d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_ipsi(idcs.passed.Xonly{i})]);
        meanAmplitudeActive.Atrials_Aonly(i) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Aonly{i}));
        meanAmplitudeActive.Xtrials_Xonly(i) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Xonly{i}));
        meanAmplitudeActive.pref(i) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Aonly{i});d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Xonly{i})]);
    end
end

% peakAmplitudeByPeakTime, meanAmplitudeByPeakTime, meanAmplitudeActiveByPeakTime
peakAmplitudeByPeakTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
peakAmplitudeByPeakTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
peakAmplitudeByPeakTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
meanAmplitudeByPeakTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
meanAmplitudeByPeakTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
meanAmplitudeByPeakTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
meanAmplitudeActiveByPeakTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
meanAmplitudeActiveByPeakTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
meanAmplitudeActiveByPeakTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(peakTime.Atrials_Aonly{i},binEdges);
        temp_X = discretize(peakTime.Xtrials_Xonly{i},binEdges);
        for n=1:length(binEdges)-1
            peakAmplitudeByPeakTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.peakAmplitude_blSub(idcs.passed.Aonly{i}(temp_A==n)));
            peakAmplitudeByPeakTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.peakAmplitude_blSub(idcs.passed.Xonly{i}(temp_X==n)));
            peakAmplitudeByPeakTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.peakAmplitude_blSub(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.peakAmplitude_blSub(idcs.passed.Xonly{i}(temp_X==n))]);
            meanAmplitudeByPeakTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_ipsi(idcs.passed.Aonly{i}(temp_A==n)));
            meanAmplitudeByPeakTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_ipsi(idcs.passed.Xonly{i}(temp_X==n)));
            meanAmplitudeByPeakTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_ipsi(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_ipsi(idcs.passed.Xonly{i}(temp_X==n))]);
            meanAmplitudeActiveByPeakTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Aonly{i}(temp_A==n)));
            meanAmplitudeActiveByPeakTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Xonly{i}(temp_X==n)));
            meanAmplitudeActiveByPeakTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end

% peakAmplitudeByComTime, meanAmplitudeByComTime, meanAmplitudeActiveByComTime
peakAmplitudeByComTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
peakAmplitudeByComTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
peakAmplitudeByComTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
meanAmplitudeByComTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
meanAmplitudeByComTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
meanAmplitudeByComTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
meanAmplitudeActiveByComTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
meanAmplitudeActiveByComTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
meanAmplitudeActiveByComTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(comTime.Atrials_Aonly{i},binEdges);
        temp_X = discretize(comTime.Xtrials_Xonly{i},binEdges);
        for n=1:length(binEdges)-1
            peakAmplitudeByComTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.peakAmplitude_blSub(idcs.passed.Aonly{i}(temp_A==n)));
            peakAmplitudeByComTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.peakAmplitude_blSub(idcs.passed.Xonly{i}(temp_X==n)));
            peakAmplitudeByComTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.peakAmplitude_blSub(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.peakAmplitude_blSub(idcs.passed.Xonly{i}(temp_X==n))]);
            meanAmplitudeByComTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_ipsi(idcs.passed.Aonly{i}(temp_A==n)));
            meanAmplitudeByComTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_ipsi(idcs.passed.Xonly{i}(temp_X==n)));
            meanAmplitudeByComTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_ipsi(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_ipsi(idcs.passed.Xonly{i}(temp_X==n))]);
            meanAmplitudeActiveByComTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Aonly{i}(temp_A==n)));
            meanAmplitudeActiveByComTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Xonly{i}(temp_X==n)));
            meanAmplitudeActiveByComTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end

% peakAmplitudeByCellType, meanAmplitudeByCellType, meanAmplitudeActiveByCellType
peakAmplitudeByCellType.Atrials_Aonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
peakAmplitudeByCellType.Xtrials_Xonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
peakAmplitudeByCellType.pref = nan(d_info.numAnimals,length(cellTypeEdges)-1);
meanAmplitudeByCellType.Atrials_Aonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
meanAmplitudeByCellType.Xtrials_Xonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
meanAmplitudeByCellType.pref = nan(d_info.numAnimals,length(cellTypeEdges)-1);
meanAmplitudeActiveByCellType.Atrials_Aonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
meanAmplitudeActiveByCellType.Xtrials_Xonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
meanAmplitudeActiveByCellType.pref = nan(d_info.numAnimals,length(cellTypeEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(peakTime.Atrials_Aonly{i},cellTypeEdges);
        temp_X = discretize(peakTime.Xtrials_Xonly{i},cellTypeEdges);
        for n=1:length(cellTypeEdges)-1
            peakAmplitudeByCellType.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.peakAmplitude_blSub(idcs.passed.Aonly{i}(temp_A==n)));
            peakAmplitudeByCellType.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.peakAmplitude_blSub(idcs.passed.Xonly{i}(temp_X==n)));
            peakAmplitudeByCellType.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.peakAmplitude_blSub(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.peakAmplitude_blSub(idcs.passed.Xonly{i}(temp_X==n))]);
            meanAmplitudeByCellType.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_ipsi(idcs.passed.Aonly{i}(temp_A==n)));
            meanAmplitudeByCellType.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_ipsi(idcs.passed.Xonly{i}(temp_X==n)));
            meanAmplitudeByCellType.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_ipsi(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_ipsi(idcs.passed.Xonly{i}(temp_X==n))]);
            meanAmplitudeActiveByCellType.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Aonly{i}(temp_A==n)));
            meanAmplitudeActiveByCellType.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Xonly{i}(temp_X==n)));
            meanAmplitudeActiveByCellType.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.meanAmplitude_blSub_active_ipsi(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end


%% Get data arrays - reliability

% activationProbability
activationProbability.Atrials_Aonly = nan(d_info.numAnimals,1);
activationProbability.Xtrials_Xonly = nan(d_info.numAnimals,1);
activationProbability.pref = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        activationProbability.Atrials_Aonly(i) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.activationProbability_ipsi(idcs.passed.Aonly{i}));
        activationProbability.Xtrials_Xonly(i) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.activationProbability_ipsi(idcs.passed.Xonly{i}));
        activationProbability.pref(i) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.activationProbability_ipsi(idcs.passed.Aonly{i});d{i,1}.(tng_struct).firingField.X_AW.activationProbability_ipsi(idcs.passed.Xonly{i})]);
    end
end

% activationProbabilityByPeakTime
activationProbabilityByPeakTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
activationProbabilityByPeakTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
activationProbabilityByPeakTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(peakTime.Atrials_Aonly{i},binEdges);
        temp_X = discretize(peakTime.Xtrials_Xonly{i},binEdges);
        for n=1:length(binEdges)-1
            activationProbabilityByPeakTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.activationProbability_ipsi(idcs.passed.Aonly{i}(temp_A==n)));
            activationProbabilityByPeakTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.activationProbability_ipsi(idcs.passed.Xonly{i}(temp_X==n)));
            activationProbabilityByPeakTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.activationProbability_ipsi(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.activationProbability_ipsi(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end

% activationProbabilityByComTime
activationProbabilityByComTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
activationProbabilityByComTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
activationProbabilityByComTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(comTime.Atrials_Aonly{i},binEdges);
        temp_X = discretize(comTime.Xtrials_Xonly{i},binEdges);
        for n=1:length(binEdges)-1
            activationProbabilityByComTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.activationProbability_ipsi(idcs.passed.Aonly{i}(temp_A==n)));
            activationProbabilityByComTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.activationProbability_ipsi(idcs.passed.Xonly{i}(temp_X==n)));
            activationProbabilityByComTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.activationProbability_ipsi(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.activationProbability_ipsi(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end

% activationProbabilityByCellType
activationProbabilityByCellType.Atrials_Aonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
activationProbabilityByCellType.Xtrials_Xonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
activationProbabilityByCellType.pref = nan(d_info.numAnimals,length(cellTypeEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(peakTime.Atrials_Aonly{i},cellTypeEdges);
        temp_X = discretize(peakTime.Xtrials_Xonly{i},cellTypeEdges);
        for n=1:length(cellTypeEdges)-1
            activationProbabilityByCellType.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.activationProbability_ipsi(idcs.passed.Aonly{i}(temp_A==n)));
            activationProbabilityByCellType.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.activationProbability_ipsi(idcs.passed.Xonly{i}(temp_X==n)));
            activationProbabilityByCellType.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.activationProbability_ipsi(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.activationProbability_ipsi(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end


%% Get data arrays - baseline

% baseline
baseline.Atrials_Aonly = nan(d_info.numAnimals,1);
baseline.Xtrials_Xonly = nan(d_info.numAnimals,1);
baseline.pref = nan(d_info.numAnimals,1);
baseline.iscell = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        baseline.Atrials_Aonly(i) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.baseline_meanOfTrialwise(idcs.passed.Aonly{i}));
        baseline.Xtrials_Xonly(i) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.baseline_meanOfTrialwise(idcs.passed.Xonly{i}));
        baseline.pref(i) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.baseline_meanOfTrialwise(idcs.passed.Aonly{i});d{i,1}.(tng_struct).firingField.X_AW.baseline_meanOfTrialwise(idcs.passed.Xonly{i})]);
        baseline.iscell(i) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.baseline_meanOfTrialwise(idcs.iscells{i}));
    end
end

% baselineByPeakTime
baselineByPeakTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
baselineByPeakTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
baselineByPeakTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(peakTime.Atrials_Aonly{i},binEdges);
        temp_X = discretize(peakTime.Xtrials_Xonly{i},binEdges);
        for n=1:length(binEdges)-1
            baselineByPeakTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.baseline_meanOfTrialwise(idcs.passed.Aonly{i}(temp_A==n)));
            baselineByPeakTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.baseline_meanOfTrialwise(idcs.passed.Xonly{i}(temp_X==n)));
            baselineByPeakTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.baseline_meanOfTrialwise(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.baseline_meanOfTrialwise(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end

% baselineByComTime
baselineByComTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
baselineByComTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
baselineByComTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(comTime.Atrials_Aonly{i},binEdges);
        temp_X = discretize(comTime.Xtrials_Xonly{i},binEdges);
        for n=1:length(binEdges)-1
            baselineByComTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.baseline_meanOfTrialwise(idcs.passed.Aonly{i}(temp_A==n)));
            baselineByComTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.baseline_meanOfTrialwise(idcs.passed.Xonly{i}(temp_X==n)));
            baselineByComTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.baseline_meanOfTrialwise(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.baseline_meanOfTrialwise(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end

% baselineByCellType
baselineByCellType.Atrials_Aonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
baselineByCellType.Xtrials_Xonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
baselineByCellType.pref = nan(d_info.numAnimals,length(cellTypeEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(peakTime.Atrials_Aonly{i},cellTypeEdges);
        temp_X = discretize(peakTime.Xtrials_Xonly{i},cellTypeEdges);
        for n=1:length(cellTypeEdges)-1
            baselineByCellType.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.baseline_meanOfTrialwise(idcs.passed.Aonly{i}(temp_A==n)));
            baselineByCellType.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.baseline_meanOfTrialwise(idcs.passed.Xonly{i}(temp_X==n)));
            baselineByCellType.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.baseline_meanOfTrialwise(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.baseline_meanOfTrialwise(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end


%% Get data arrays - peak time sd

% peakTimeStd
peakTimeStd.Atrials_Aonly = nan(d_info.numAnimals,1);
peakTimeStd.Xtrials_Xonly = nan(d_info.numAnimals,1);
peakTimeStd.pref = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        peakTimeStd.Atrials_Aonly(i) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.peakTimeSd_s(idcs.passed.Aonly{i}));
        peakTimeStd.Xtrials_Xonly(i) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.peakTimeSd_s(idcs.passed.Xonly{i}));
        peakTimeStd.pref(i) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.peakTimeSd_s(idcs.passed.Aonly{i});d{i,1}.(tng_struct).firingField.X_AW.peakTimeSd_s(idcs.passed.Xonly{i})]);
    end
end

% peakTimeStdByPeakTime
peakTimeStdByPeakTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
peakTimeStdByPeakTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
peakTimeStdByPeakTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(peakTime.Atrials_Aonly{i},binEdges);
        temp_X = discretize(peakTime.Xtrials_Xonly{i},binEdges);
        for n=1:length(binEdges)-1
            peakTimeStdByPeakTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.peakTimeSd_s(idcs.passed.Aonly{i}(temp_A==n)));
            peakTimeStdByPeakTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.peakTimeSd_s(idcs.passed.Xonly{i}(temp_X==n)));
            peakTimeStdByPeakTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.peakTimeSd_s(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.peakTimeSd_s(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end

% peakTimeStdByComTime
peakTimeStdByComTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
peakTimeStdByComTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
peakTimeStdByComTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(comTime.Atrials_Aonly{i},binEdges);
        temp_X = discretize(comTime.Xtrials_Xonly{i},binEdges);
        for n=1:length(binEdges)-1
            peakTimeStdByComTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.peakTimeSd_s(idcs.passed.Aonly{i}(temp_A==n)));
            peakTimeStdByComTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.peakTimeSd_s(idcs.passed.Xonly{i}(temp_X==n)));
            peakTimeStdByComTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.peakTimeSd_s(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.peakTimeSd_s(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end

% peakTimeStdByCellType
peakTimeStdByCellType.Atrials_Aonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
peakTimeStdByCellType.Xtrials_Xonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
peakTimeStdByCellType.pref = nan(d_info.numAnimals,length(cellTypeEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(peakTime.Atrials_Aonly{i},cellTypeEdges);
        temp_X = discretize(peakTime.Xtrials_Xonly{i},cellTypeEdges);
        for n=1:length(cellTypeEdges)-1
            peakTimeStdByCellType.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.peakTimeSd_s(idcs.passed.Aonly{i}(temp_A==n)));
            peakTimeStdByCellType.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.peakTimeSd_s(idcs.passed.Xonly{i}(temp_X==n)));
            peakTimeStdByCellType.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.peakTimeSd_s(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.peakTimeSd_s(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end


%% Get data arrays - selectivity

% selectivity
selectivity.Atrials_Aonly = nan(d_info.numAnimals,1);
selectivity.Xtrials_Xonly = nan(d_info.numAnimals,1);
selectivity.pref = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        selectivity.Atrials_Aonly(i) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.selectivity_ipsi(idcs.passed.Aonly{i}));
        selectivity.Xtrials_Xonly(i) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.selectivity_ipsi(idcs.passed.Xonly{i}));
        selectivity.pref(i) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.selectivity_ipsi(idcs.passed.Aonly{i});d{i,1}.(tng_struct).firingField.X_AW.selectivity_ipsi(idcs.passed.Xonly{i})]);
    end
end

% selectivityByPeakTime
selectivityByPeakTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
selectivityByPeakTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
selectivityByPeakTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(peakTime.Atrials_Aonly{i},binEdges);
        temp_X = discretize(peakTime.Xtrials_Xonly{i},binEdges);
        for n=1:length(binEdges)-1
            selectivityByPeakTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.selectivity_ipsi(idcs.passed.Aonly{i}(temp_A==n)));
            selectivityByPeakTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.selectivity_ipsi(idcs.passed.Xonly{i}(temp_X==n)));
            selectivityByPeakTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.selectivity_ipsi(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.selectivity_ipsi(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end

% selectivityByComTime
selectivityByComTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
selectivityByComTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
selectivityByComTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(comTime.Atrials_Aonly{i},binEdges);
        temp_X = discretize(comTime.Xtrials_Xonly{i},binEdges);
        for n=1:length(binEdges)-1
            selectivityByComTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.selectivity_ipsi(idcs.passed.Aonly{i}(temp_A==n)));
            selectivityByComTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.selectivity_ipsi(idcs.passed.Xonly{i}(temp_X==n)));
            selectivityByComTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.selectivity_ipsi(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.selectivity_ipsi(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end

% selectivityByCellType
selectivityByCellType.Atrials_Aonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
selectivityByCellType.Xtrials_Xonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
selectivityByCellType.pref = nan(d_info.numAnimals,length(cellTypeEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(peakTime.Atrials_Aonly{i},cellTypeEdges);
        temp_X = discretize(peakTime.Xtrials_Xonly{i},cellTypeEdges);
        for n=1:length(cellTypeEdges)-1
            selectivityByCellType.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.selectivity_ipsi(idcs.passed.Aonly{i}(temp_A==n)));
            selectivityByCellType.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.selectivity_ipsi(idcs.passed.Xonly{i}(temp_X==n)));
            selectivityByCellType.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.selectivity_ipsi(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.selectivity_ipsi(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end


%% Get data arrays - width

% width
width.Atrials_Aonly = nan(d_info.numAnimals,1);
width.Xtrials_Xonly = nan(d_info.numAnimals,1);
width.pref = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        width.Atrials_Aonly(i) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.width_s(idcs.passed.Aonly{i}));
        width.Xtrials_Xonly(i) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.width_s(idcs.passed.Xonly{i}));
        width.pref(i) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.width_s(idcs.passed.Aonly{i});d{i,1}.(tng_struct).firingField.X_AW.width_s(idcs.passed.Xonly{i})]);
    end
end

% widthByPeakTime
widthByPeakTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
widthByPeakTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
widthByPeakTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(peakTime.Atrials_Aonly{i},binEdges);
        temp_X = discretize(peakTime.Xtrials_Xonly{i},binEdges);
        for n=1:length(binEdges)-1
            widthByPeakTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.width_s(idcs.passed.Aonly{i}(temp_A==n)));
            widthByPeakTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.width_s(idcs.passed.Xonly{i}(temp_X==n)));
            widthByPeakTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.width_s(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.width_s(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end

% widthByComTime
widthByComTime.Atrials_Aonly = nan(d_info.numAnimals,length(binEdges)-1);
widthByComTime.Xtrials_Xonly = nan(d_info.numAnimals,length(binEdges)-1);
widthByComTime.pref = nan(d_info.numAnimals,length(binEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(comTime.Atrials_Aonly{i},binEdges);
        temp_X = discretize(comTime.Xtrials_Xonly{i},binEdges);
        for n=1:length(binEdges)-1
            widthByComTime.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.width_s(idcs.passed.Aonly{i}(temp_A==n)));
            widthByComTime.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.width_s(idcs.passed.Xonly{i}(temp_X==n)));
            widthByComTime.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.width_s(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.width_s(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end

% widthByCellType
widthByCellType.Atrials_Aonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
widthByCellType.Xtrials_Xonly = nan(d_info.numAnimals,length(cellTypeEdges)-1);
widthByCellType.pref = nan(d_info.numAnimals,length(cellTypeEdges)-1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        temp_A = discretize(peakTime.Atrials_Aonly{i},cellTypeEdges);
        temp_X = discretize(peakTime.Xtrials_Xonly{i},cellTypeEdges);
        for n=1:length(cellTypeEdges)-1
            widthByCellType.Atrials_Aonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.A_AW.width_s(idcs.passed.Aonly{i}(temp_A==n)));
            widthByCellType.Xtrials_Xonly(i,n) = nanmean(d{i,1}.(tng_struct).firingField.X_AW.width_s(idcs.passed.Xonly{i}(temp_X==n)));
            widthByCellType.pref(i,n) = nanmean([d{i,1}.(tng_struct).firingField.A_AW.width_s(idcs.passed.Aonly{i}(temp_A==n));d{i,1}.(tng_struct).firingField.X_AW.width_s(idcs.passed.Xonly{i}(temp_X==n))]);
        end
	end
end


%% Get data arrays - Warping

% exp_a
exp_a = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        exp_a(i) = d{i,1}.(warp_struct).peak.exp1.a;
    end
end

% exp_b
exp_b = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    if ~isempty(idcs.iscells{i})
        exp_b(i) = d{i,1}.(warp_struct).peak.exp1.b;
    end
end


%% Get data arrays - Bayesian decoding

timeDecodingError_s = nan(d_info.numAnimals,25);
timeDecodingError_withinCat_s = nan(d_info.numAnimals,25);
typeDecodingCorrect = nan(d_info.numAnimals,25);
for i=1:d_info.numAnimals
    try
        timeDecodingError_s(i,:) = nanmean(d{i}.dec_all.analysis.trainingSet_seq.timeDecodingError_s,1);
        timeDecodingError_withinCat_s(i,:) = nanmean(d{i}.dec_all.analysis.trainingSet_seq.timeDecodingError_withinCat_s,1);
        typeDecodingCorrect(i,:) = nanmean(d{i}.dec_all.analysis.trainingSet_seq.typeDecodingCorrect,1);
    catch
    end
end

timeDecodingError_s_allCells = nan(d_info.numAnimals,25);
timeDecodingError_withinCat_s_allCells = nan(d_info.numAnimals,25);
typeDecodingCorrect_allCells = nan(d_info.numAnimals,25);
for i=1:d_info.numAnimals
    try
        timeDecodingError_s_allCells(i,:) = nanmean(d{i}.dec_all.analysis.trainingSet.timeDecodingError_s,1);
        timeDecodingError_withinCat_s_allCells(i,:) = nanmean(d{i}.dec_all.analysis.trainingSet.timeDecodingError_withinCat_s,1);
        typeDecodingCorrect_allCells(i,:) = nanmean(d{i}.dec_all.analysis.trainingSet.typeDecodingCorrect,1);
    catch
    end
end


%% Get bcon

numTimeBins = 67;

distance = nan(d_info.numAnimals,numTimeBins);
velocity = nan(d_info.numAnimals,numTimeBins);
acceleration = nan(d_info.numAnimals,numTimeBins);
licking = nan(d_info.numAnimals,numTimeBins);
distance_A = nan(d_info.numAnimals,numTimeBins);
velocity_A = nan(d_info.numAnimals,numTimeBins);
acceleration_A = nan(d_info.numAnimals,numTimeBins);
licking_A = nan(d_info.numAnimals,numTimeBins);
distance_X = nan(d_info.numAnimals,numTimeBins);
velocity_X = nan(d_info.numAnimals,numTimeBins);
acceleration_X = nan(d_info.numAnimals,numTimeBins);
licking_X = nan(d_info.numAnimals,numTimeBins);
for i=1:d_info.numAnimals
    try
        distance(i,:) = nanmean(d{i,1}.bcon.binwise_full.distance,1);
        velocity(i,:) = nanmean(d{i,1}.bcon.binwise_full.velocity,1);
        acceleration(i,:) = nanmean(d{i,1}.bcon.binwise_full.acceleration,1);
        licking(i,:) = nanmean(d{i,1}.bcon.binwise_full.licking,1);
        distance_A(i,:) = nanmean(d{i,1}.bcon.binwise_full.distance(find(d{i,1}.task.odour1=="A"),:),1);
        velocity_A(i,:) = nanmean(d{i,1}.bcon.binwise_full.velocity(find(d{i,1}.task.odour1=="A"),:),1);
        acceleration_A(i,:) = nanmean(d{i,1}.bcon.binwise_full.acceleration(find(d{i,1}.task.odour1=="A"),:),1);
        licking_A(i,:) = nanmean(d{i,1}.bcon.binwise_full.licking(find(d{i,1}.task.odour1=="A"),:),1);
        distance_X(i,:) = nanmean(d{i,1}.bcon.binwise_full.distance(find(d{i,1}.task.odour1=="X"),:),1);
        velocity_X(i,:) = nanmean(d{i,1}.bcon.binwise_full.velocity(find(d{i,1}.task.odour1=="X"),:),1);
        acceleration_X(i,:) = nanmean(d{i,1}.bcon.binwise_full.acceleration(find(d{i,1}.task.odour1=="X"),:),1);
        licking_X(i,:) = nanmean(d{i,1}.bcon.binwise_full.licking(find(d{i,1}.task.odour1=="X"),:),1);
    catch
    end
end


%% Calculate trial-by-trial correlations

tbtcorr_correct_typeDecodingCorrect_late = nan(d_info.numAnimals,1);
for i=1:d_info.numAnimals
    try
        this_correct = (d{i}.task.response=="H" | d{i}.task.response=="CR")';
        this_decoding = nanmean(d{i}.dec_all.analysis.trainingSet_seq.typeDecodingCorrect(:,17:24),2);
        
        tbtcorr_correct_typeDecodingCorrect_late(i) = corr(this_correct,this_decoding,'Type','Pearson','Rows','Complete');
    catch
    end
end

%figure;histogram(tbtcorr_correct_typeDecodingCorrect_late);
%signrank(tbtcorr_correct_typeDecodingCorrect_late)






