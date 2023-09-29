%% Fig3_WarpingAndBehaviourConfounds

% import data using Summary_Master with ops.do_sequenceDevelopmentSummary = true;
% d_info = selectDataset(d_info,'-g2-d12345-e01-r01-p01-l01-i00',sheet,path,ops);

info = get_info;
save_root_fig = [path.root_summary,'figures\Fig3_fig\'];
save_root_png = [path.root_summary,'figures\Fig3_png\'];
save_root_pdf = [path.root_summary,'figures\Fig3_pdf\'];
save_root_txt = [path.root_summary,'figures\Fig3_txt\'];

blocks = [5,8,5,8,5];
numBlocks = sum(blocks);
blocks_switch = [8,5];
numBlocks_switch = sum(blocks_switch);

numTimeBins = 67;
binEdges = 0:5/3:5; %0:5.3/3:5;
numBins = length(binEdges)-1;
smallBinEdges = 0:0.2:5;
numSmallBins = length(smallBinEdges)-1;


%% Get running

this_running_cleaned = d_info.running;
for i=1:d_info.numAnimals
    if this_running_cleaned(i,2) ~= this_running_cleaned(i,3)
        this_running_cleaned(i,2) = NaN;
        this_running_cleaned(i,3) = NaN;
    end
    if this_running_cleaned(i,4) ~= this_running_cleaned(i,5)
        this_running_cleaned(i,4) = NaN;
        this_running_cleaned(i,5) = NaN;
    end
end

this_running = nan(d_info.numAnimals,numBlocks);
for j=1:5
    if j==1
        this_running(:,1:5) = repmat(this_running_cleaned(:,j),1,5);
    elseif j==2
        this_running(:,6:13) = repmat(this_running_cleaned(:,j),1,8);
    elseif j==3
        this_running(:,14:18) = repmat(this_running_cleaned(:,j),1,5);
    elseif j==4
        this_running(:,19:26) = repmat(this_running_cleaned(:,j),1,8);
    elseif j==5
        this_running(:,27:31) = repmat(this_running_cleaned(:,j),1,5);
    end
end


%% Get performance

correct = nan(d_info.numAnimals,numBlocks);
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
                correct(i,this_block) = nanmean(d{i,j}.perf.blocks_general.correct(these_blocks_20t));
            end
        end
    end
end

% correct split by running
correct_runner = nan(d_info.numAnimals,numBlocks);
correct_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            correct_runner(i,k) = correct(i,k);
        elseif this_running(i,k)==0
            correct_nonrunner(i,k) = correct(i,k);
        end
    end
end
correct_switch_runner = nanmean(cat(3,correct_runner(:,6:18),correct_runner(:,19:31)),3);
correct_switch_nonrunner = nanmean(cat(3,correct_nonrunner(:,6:18),correct_nonrunner(:,19:31)),3);


%% Get bcon

% for each 100t block across all 5 days
distance = nan(numTimeBins,numBlocks,d_info.numAnimals); % [time bin, block, animal]
velocity = nan(numTimeBins,numBlocks,d_info.numAnimals); % [time bin, block, animal]
acceleration = nan(numTimeBins,numBlocks,d_info.numAnimals); % [time bin, block, animal]
licking = nan(numTimeBins,numBlocks,d_info.numAnimals); % [time bin, block, animal]
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
                
                distance(:,this_block,i) = nanmean(d{i,j}.bcon.binwise_full.distance((k-1)*100+1:k*100,:),1);
                velocity(:,this_block,i) = nanmean(d{i,j}.bcon.binwise_full.velocity((k-1)*100+1:k*100,:),1);
                acceleration(:,this_block,i) = nanmean(d{i,j}.bcon.binwise_full.acceleration((k-1)*100+1:k*100,:),1);
                licking(:,this_block,i) = nanmean(d{i,j}.bcon.binwise_full.licking((k-1)*100+1:k*100,:),1);
            end
        catch
        end
    end
end

% distance split by running
distance_runner = nan(numTimeBins,numBlocks,d_info.numAnimals);
distance_nonrunner = nan(numTimeBins,numBlocks,d_info.numAnimals);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            distance_runner(:,k,i) = distance(:,k,i);
        elseif this_running(i,k)==0
            distance_nonrunner(:,k,i) = distance(:,k,i);
        end
    end
end
distance_switch_runner = nanmean(cat(4,distance_runner(:,6:18,:),distance_runner(:,19:31,:)),4);
distance_switch_nonrunner = nanmean(cat(4,distance_nonrunner(:,6:18,:),distance_nonrunner(:,19:31,:)),4);

% velocity split by running
velocity_runner = nan(numTimeBins,numBlocks,d_info.numAnimals);
velocity_nonrunner = nan(numTimeBins,numBlocks,d_info.numAnimals);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            velocity_runner(:,k,i) = velocity(:,k,i);
        elseif this_running(i,k)==0
            velocity_nonrunner(:,k,i) = velocity(:,k,i);
        end
    end
end
velocity_switch_runner = nanmean(cat(4,velocity_runner(:,6:18,:),velocity_runner(:,19:31,:)),4);
velocity_switch_nonrunner = nanmean(cat(4,velocity_nonrunner(:,6:18,:),velocity_nonrunner(:,19:31,:)),4);

% acceleration split by running
acceleration_runner = nan(numTimeBins,numBlocks,d_info.numAnimals);
acceleration_nonrunner = nan(numTimeBins,numBlocks,d_info.numAnimals);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            acceleration_runner(:,k,i) = acceleration(:,k,i);
        elseif this_running(i,k)==0
            acceleration_nonrunner(:,k,i) = acceleration(:,k,i);
        end
    end
end
acceleration_switch_runner = nanmean(cat(4,acceleration_runner(:,6:18,:),acceleration_runner(:,19:31,:)),4);
acceleration_switch_nonrunner = nanmean(cat(4,acceleration_nonrunner(:,6:18,:),acceleration_nonrunner(:,19:31,:)),4);

% licking split by running
licking_runner = nan(numTimeBins,numBlocks,d_info.numAnimals);
licking_nonrunner = nan(numTimeBins,numBlocks,d_info.numAnimals);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            licking_runner(:,k,i) = licking(:,k,i);
        elseif this_running(i,k)==0
            licking_nonrunner(:,k,i) = licking(:,k,i);
        end
    end
end
licking_switch_runner = nanmean(cat(4,licking_runner(:,6:18,:),licking_runner(:,19:31,:)),4);
licking_switch_nonrunner = nanmean(cat(4,licking_nonrunner(:,6:18,:),licking_nonrunner(:,19:31,:)),4);


%% Get data arrays - Bayesian decoding

% for each 100t block across all 5 days
typeDecodingCorrect_E = nan(d_info.numAnimals,numBlocks);
typeDecodingCorrect_I = nan(d_info.numAnimals,numBlocks);
typeDecodingCorrect_L = nan(d_info.numAnimals,numBlocks);
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
                
                typeDecodingCorrect_E(i,this_block) = nanmean(nanmean(d{i,j}.dec_100T{k}.analysis.trainingSet_seq.typeDecodingCorrect(:,1:8),2));
                typeDecodingCorrect_I(i,this_block) = nanmean(nanmean(d{i,j}.dec_100T{k}.analysis.trainingSet_seq.typeDecodingCorrect(:,9:16),2));
                typeDecodingCorrect_L(i,this_block) = nanmean(nanmean(d{i,j}.dec_100T{k}.analysis.trainingSet_seq.typeDecodingCorrect(:,17:24),2));
            end
        catch
        end
    end
end

% typeDecodingCorrect_E split by running
typeDecodingCorrect_E_runner = nan(d_info.numAnimals,numBlocks);
typeDecodingCorrect_E_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            typeDecodingCorrect_E_runner(i,k) = typeDecodingCorrect_E(i,k);
        elseif this_running(i,k)==0
            typeDecodingCorrect_E_nonrunner(i,k) = typeDecodingCorrect_E(i,k);
        end
    end
end
typeDecodingCorrect_E_switch_runner = nanmean(cat(3,typeDecodingCorrect_E_runner(:,6:18),typeDecodingCorrect_E_runner(:,19:31)),3);
typeDecodingCorrect_E_switch_nonrunner = nanmean(cat(3,typeDecodingCorrect_E_nonrunner(:,6:18),typeDecodingCorrect_E_nonrunner(:,19:31)),3);

% typeDecodingCorrect_I split by running
typeDecodingCorrect_I_runner = nan(d_info.numAnimals,numBlocks);
typeDecodingCorrect_I_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            typeDecodingCorrect_I_runner(i,k) = typeDecodingCorrect_I(i,k);
        elseif this_running(i,k)==0
            typeDecodingCorrect_I_nonrunner(i,k) = typeDecodingCorrect_I(i,k);
        end
    end
end
typeDecodingCorrect_I_switch_runner = nanmean(cat(3,typeDecodingCorrect_I_runner(:,6:18),typeDecodingCorrect_I_runner(:,19:31)),3);
typeDecodingCorrect_I_switch_nonrunner = nanmean(cat(3,typeDecodingCorrect_I_nonrunner(:,6:18),typeDecodingCorrect_I_nonrunner(:,19:31)),3);

% typeDecodingCorrect_L split by running
typeDecodingCorrect_L_runner = nan(d_info.numAnimals,numBlocks);
typeDecodingCorrect_L_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            typeDecodingCorrect_L_runner(i,k) = typeDecodingCorrect_E(i,k);
        elseif this_running(i,k)==0
            typeDecodingCorrect_L_nonrunner(i,k) = typeDecodingCorrect_E(i,k);
        end
    end
end
typeDecodingCorrect_L_switch_runner = nanmean(cat(3,typeDecodingCorrect_L_runner(:,6:18),typeDecodingCorrect_L_runner(:,19:31)),3);
typeDecodingCorrect_L_switch_nonrunner = nanmean(cat(3,typeDecodingCorrect_L_nonrunner(:,6:18),typeDecodingCorrect_L_nonrunner(:,19:31)),3);

% for each 100t block across all 5 days
typeDecodingCorrect_E_allCells = nan(d_info.numAnimals,numBlocks);
typeDecodingCorrect_I_allCells = nan(d_info.numAnimals,numBlocks);
typeDecodingCorrect_L_allCells = nan(d_info.numAnimals,numBlocks);
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
                
                typeDecodingCorrect_E_allCells(i,this_block) = nanmean(nanmean(d{i,j}.dec_100T{k}.analysis.trainingSet.typeDecodingCorrect(:,1:8),2));
                typeDecodingCorrect_I_allCells(i,this_block) = nanmean(nanmean(d{i,j}.dec_100T{k}.analysis.trainingSet.typeDecodingCorrect(:,9:16),2));
                typeDecodingCorrect_L_allCells(i,this_block) = nanmean(nanmean(d{i,j}.dec_100T{k}.analysis.trainingSet.typeDecodingCorrect(:,17:24),2));
            end
        catch
        end
    end
end

% typeDecodingCorrect_E_allCells split by running
typeDecodingCorrect_E_allCells_runner = nan(d_info.numAnimals,numBlocks);
typeDecodingCorrect_E_allCells_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            typeDecodingCorrect_E_allCells_runner(i,k) = typeDecodingCorrect_E_allCells(i,k);
        elseif this_running(i,k)==0
            typeDecodingCorrect_E_allCells_nonrunner(i,k) = typeDecodingCorrect_E_allCells(i,k);
        end
    end
end
typeDecodingCorrect_E_allCells_switch_runner = nanmean(cat(3,typeDecodingCorrect_E_allCells_runner(:,6:18),typeDecodingCorrect_E_allCells_runner(:,19:31)),3);
typeDecodingCorrect_E_allCells_switch_nonrunner = nanmean(cat(3,typeDecodingCorrect_E_allCells_nonrunner(:,6:18),typeDecodingCorrect_E_allCells_nonrunner(:,19:31)),3);

% typeDecodingCorrect_I_allCells split by running
typeDecodingCorrect_I_allCells_runner = nan(d_info.numAnimals,numBlocks);
typeDecodingCorrect_I_allCells_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            typeDecodingCorrect_I_allCells_runner(i,k) = typeDecodingCorrect_I_allCells(i,k);
        elseif this_running(i,k)==0
            typeDecodingCorrect_I_allCells_nonrunner(i,k) = typeDecodingCorrect_I_allCells(i,k);
        end
    end
end
typeDecodingCorrect_I_allCells_switch_runner = nanmean(cat(3,typeDecodingCorrect_I_allCells_runner(:,6:18),typeDecodingCorrect_I_allCells_runner(:,19:31)),3);
typeDecodingCorrect_I_allCells_switch_nonrunner = nanmean(cat(3,typeDecodingCorrect_I_allCells_nonrunner(:,6:18),typeDecodingCorrect_I_allCells_nonrunner(:,19:31)),3);

% typeDecodingCorrect_L_allCells split by running
typeDecodingCorrect_L_allCells_runner = nan(d_info.numAnimals,numBlocks);
typeDecodingCorrect_L_allCells_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            typeDecodingCorrect_L_allCells_runner(i,k) = typeDecodingCorrect_E_allCells(i,k);
        elseif this_running(i,k)==0
            typeDecodingCorrect_L_allCells_nonrunner(i,k) = typeDecodingCorrect_E_allCells(i,k);
        end
    end
end
typeDecodingCorrect_L_allCells_switch_runner = nanmean(cat(3,typeDecodingCorrect_L_allCells_runner(:,6:18),typeDecodingCorrect_L_allCells_runner(:,19:31)),3);
typeDecodingCorrect_L_allCells_switch_nonrunner = nanmean(cat(3,typeDecodingCorrect_L_allCells_nonrunner(:,6:18),typeDecodingCorrect_L_allCells_nonrunner(:,19:31)),3);


% for each 100t block across all 5 days
timeDecodingError_s_E = nan(d_info.numAnimals,numBlocks);
timeDecodingError_s_I = nan(d_info.numAnimals,numBlocks);
timeDecodingError_s_L = nan(d_info.numAnimals,numBlocks);
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
                
                timeDecodingError_s_E(i,this_block) = nanmean(nanmean(d{i,j}.dec_100T{k}.analysis.trainingSet_seq.timeDecodingError_s(:,1:8),2));
                timeDecodingError_s_I(i,this_block) = nanmean(nanmean(d{i,j}.dec_100T{k}.analysis.trainingSet_seq.timeDecodingError_s(:,9:16),2));
                timeDecodingError_s_L(i,this_block) = nanmean(nanmean(d{i,j}.dec_100T{k}.analysis.trainingSet_seq.timeDecodingError_s(:,17:24),2));
            end
        catch
        end
    end
end

% timeDecodingError_s_E split by running
timeDecodingError_s_E_runner = nan(d_info.numAnimals,numBlocks);
timeDecodingError_s_E_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            timeDecodingError_s_E_runner(i,k) = timeDecodingError_s_E(i,k);
        elseif this_running(i,k)==0
            timeDecodingError_s_E_nonrunner(i,k) = timeDecodingError_s_E(i,k);
        end
    end
end
timeDecodingError_s_E_switch_runner = nanmean(cat(3,timeDecodingError_s_E_runner(:,6:18),timeDecodingError_s_E_runner(:,19:31)),3);
timeDecodingError_s_E_switch_nonrunner = nanmean(cat(3,timeDecodingError_s_E_nonrunner(:,6:18),timeDecodingError_s_E_nonrunner(:,19:31)),3);

% timeDecodingError_s_I split by running
timeDecodingError_s_I_runner = nan(d_info.numAnimals,numBlocks);
timeDecodingError_s_I_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            timeDecodingError_s_I_runner(i,k) = timeDecodingError_s_I(i,k);
        elseif this_running(i,k)==0
            timeDecodingError_s_I_nonrunner(i,k) = timeDecodingError_s_I(i,k);
        end
    end
end
timeDecodingError_s_I_switch_runner = nanmean(cat(3,timeDecodingError_s_I_runner(:,6:18),timeDecodingError_s_I_runner(:,19:31)),3);
timeDecodingError_s_I_switch_nonrunner = nanmean(cat(3,timeDecodingError_s_I_nonrunner(:,6:18),timeDecodingError_s_I_nonrunner(:,19:31)),3);

% timeDecodingError_s_L split by running
timeDecodingError_s_L_runner = nan(d_info.numAnimals,numBlocks);
timeDecodingError_s_L_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            timeDecodingError_s_L_runner(i,k) = timeDecodingError_s_E(i,k);
        elseif this_running(i,k)==0
            timeDecodingError_s_L_nonrunner(i,k) = timeDecodingError_s_E(i,k);
        end
    end
end
timeDecodingError_s_L_switch_runner = nanmean(cat(3,timeDecodingError_s_L_runner(:,6:18),timeDecodingError_s_L_runner(:,19:31)),3);
timeDecodingError_s_L_switch_nonrunner = nanmean(cat(3,timeDecodingError_s_L_nonrunner(:,6:18),timeDecodingError_s_L_nonrunner(:,19:31)),3);

% for each 100t block across all 5 days
timeDecodingError_s_E_allCells = nan(d_info.numAnimals,numBlocks);
timeDecodingError_s_I_allCells = nan(d_info.numAnimals,numBlocks);
timeDecodingError_s_L_allCells = nan(d_info.numAnimals,numBlocks);
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
                
                timeDecodingError_s_E_allCells(i,this_block) = nanmean(nanmean(d{i,j}.dec_100T{k}.analysis.trainingSet.timeDecodingError_s(:,1:8),2));
                timeDecodingError_s_I_allCells(i,this_block) = nanmean(nanmean(d{i,j}.dec_100T{k}.analysis.trainingSet.timeDecodingError_s(:,9:16),2));
                timeDecodingError_s_L_allCells(i,this_block) = nanmean(nanmean(d{i,j}.dec_100T{k}.analysis.trainingSet.timeDecodingError_s(:,17:24),2));
            end
        catch
        end
    end
end

% timeDecodingError_s_E_allCells split by running
timeDecodingError_s_E_allCells_runner = nan(d_info.numAnimals,numBlocks);
timeDecodingError_s_E_allCells_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            timeDecodingError_s_E_allCells_runner(i,k) = timeDecodingError_s_E_allCells(i,k);
        elseif this_running(i,k)==0
            timeDecodingError_s_E_allCells_nonrunner(i,k) = timeDecodingError_s_E_allCells(i,k);
        end
    end
end
timeDecodingError_s_E_allCells_switch_runner = nanmean(cat(3,timeDecodingError_s_E_allCells_runner(:,6:18),timeDecodingError_s_E_allCells_runner(:,19:31)),3);
timeDecodingError_s_E_allCells_switch_nonrunner = nanmean(cat(3,timeDecodingError_s_E_allCells_nonrunner(:,6:18),timeDecodingError_s_E_allCells_nonrunner(:,19:31)),3);

% timeDecodingError_s_I_allCells split by running
timeDecodingError_s_I_allCells_runner = nan(d_info.numAnimals,numBlocks);
timeDecodingError_s_I_allCells_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            timeDecodingError_s_I_allCells_runner(i,k) = timeDecodingError_s_I_allCells(i,k);
        elseif this_running(i,k)==0
            timeDecodingError_s_I_allCells_nonrunner(i,k) = timeDecodingError_s_I_allCells(i,k);
        end
    end
end
timeDecodingError_s_I_allCells_switch_runner = nanmean(cat(3,timeDecodingError_s_I_allCells_runner(:,6:18),timeDecodingError_s_I_allCells_runner(:,19:31)),3);
timeDecodingError_s_I_allCells_switch_nonrunner = nanmean(cat(3,timeDecodingError_s_I_allCells_nonrunner(:,6:18),timeDecodingError_s_I_allCells_nonrunner(:,19:31)),3);

% timeDecodingError_s_L_allCells split by running
timeDecodingError_s_L_allCells_runner = nan(d_info.numAnimals,numBlocks);
timeDecodingError_s_L_allCells_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            timeDecodingError_s_L_allCells_runner(i,k) = timeDecodingError_s_E_allCells(i,k);
        elseif this_running(i,k)==0
            timeDecodingError_s_L_allCells_nonrunner(i,k) = timeDecodingError_s_E_allCells(i,k);
        end
    end
end
timeDecodingError_s_L_allCells_switch_runner = nanmean(cat(3,timeDecodingError_s_L_allCells_runner(:,6:18),timeDecodingError_s_L_allCells_runner(:,19:31)),3);
timeDecodingError_s_L_allCells_switch_nonrunner = nanmean(cat(3,timeDecodingError_s_L_allCells_nonrunner(:,6:18),timeDecodingError_s_L_allCells_nonrunner(:,19:31)),3);


%% Get sequence cells by peak time big bin

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

% seqCellsByTime_prop split by running
seqCellsByTime_prop_runner = nan(3,numBlocks,d_info.numAnimals);
seqCellsByTime_prop_nonrunner = nan(3,numBlocks,d_info.numAnimals);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            seqCellsByTime_prop_runner(:,k,i) = seqCellsByTime_prop(:,k,i);
        elseif this_running(i,k)==0
            seqCellsByTime_prop_nonrunner(:,k,i) = seqCellsByTime_prop(:,k,i);
        end
    end
end
seqCellsByTime_prop_switch_runner = nanmean(cat(4,seqCellsByTime_prop_runner(:,6:18,:),seqCellsByTime_prop_runner(:,19:31,:)),4);
seqCellsByTime_prop_switch_nonrunner = nanmean(cat(4,seqCellsByTime_prop_nonrunner(:,6:18,:),seqCellsByTime_prop_nonrunner(:,19:31,:)),4);

% seqCellsByTime_num split by running
seqCellsByTime_num_runner = nan(3,numBlocks,d_info.numAnimals);
seqCellsByTime_num_nonrunner = nan(3,numBlocks,d_info.numAnimals);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            seqCellsByTime_num_runner(:,k,i) = seqCellsByTime_num(:,k,i);
        elseif this_running(i,k)==0
            seqCellsByTime_num_nonrunner(:,k,i) = seqCellsByTime_num(:,k,i);
        end
    end
end
seqCellsByTime_num_switch_runner = nanmean(cat(4,seqCellsByTime_num_runner(:,6:18,:),seqCellsByTime_num_runner(:,19:31,:)),4);
seqCellsByTime_num_switch_nonrunner = nanmean(cat(4,seqCellsByTime_num_nonrunner(:,6:18,:),seqCellsByTime_num_nonrunner(:,19:31,:)),4);


%% Get sequence cells by peak time small bin

% for each 100t block across all 5 days
seqCellsByTime_sb_num = nan(numSmallBins,numBlocks,d_info.numAnimals); % [small time bin, block, animal]
seqCellsByTime_sb_prop = nan(numSmallBins,numBlocks,d_info.numAnimals); % [small time bin, block, animal]
seqCellsByTime_sb_propseq = nan(numSmallBins,numBlocks,d_info.numAnimals); % [small time bin, block, animal]
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
                temp = discretize(these_peakTimes,smallBinEdges);
                for n=1:numSmallBins
                    seqCellsByTime_sb_num(n,this_block,i) = nansum(temp==n);
                    seqCellsByTime_sb_prop(n,this_block,i) = nansum(temp==n) / this_numCells;
                    seqCellsByTime_sb_propseq(n,this_block,i) = nansum(temp==n) / length(these_peakTimes);
                end
            end
        catch
        end
    end
end
sequenceCells_num_avg = nanmean(seqCellsByTime_sb_num,3);
sequenceCells_prop_avg = nanmean(seqCellsByTime_sb_prop,3);
sequenceCells_propseq_avg = nanmean(seqCellsByTime_sb_propseq,3);

% seqCellsByTime_sb_prop split by running
seqCellsByTime_sb_prop_runner = nan(numSmallBins,numBlocks,d_info.numAnimals);
seqCellsByTime_sb_prop_nonrunner = nan(numSmallBins,numBlocks,d_info.numAnimals);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            seqCellsByTime_sb_prop_runner(:,k,i) = seqCellsByTime_sb_prop(:,k,i);
        elseif this_running(i,k)==0
            seqCellsByTime_sb_prop_nonrunner(:,k,i) = seqCellsByTime_sb_prop(:,k,i);
        end
    end
end
seqCellsByTime_sb_prop_switch_runner = nanmean(cat(4,seqCellsByTime_sb_prop_runner(:,6:18,:),seqCellsByTime_sb_prop_runner(:,19:31,:)),4);
seqCellsByTime_sb_prop_switch_nonrunner = nanmean(cat(4,seqCellsByTime_sb_prop_nonrunner(:,6:18,:),seqCellsByTime_sb_prop_nonrunner(:,19:31,:)),4);

% seqCellsByTime_sb_propseq split by running
seqCellsByTime_sb_propseq_runner = nan(numSmallBins,numBlocks,d_info.numAnimals);
seqCellsByTime_sb_propseq_nonrunner = nan(numSmallBins,numBlocks,d_info.numAnimals);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            seqCellsByTime_sb_propseq_runner(:,k,i) = seqCellsByTime_sb_propseq(:,k,i);
        elseif this_running(i,k)==0
            seqCellsByTime_sb_propseq_nonrunner(:,k,i) = seqCellsByTime_sb_propseq(:,k,i);
        end
    end
end
seqCellsByTime_sb_propseq_switch_runner = nanmean(cat(4,seqCellsByTime_sb_propseq_runner(:,6:18,:),seqCellsByTime_sb_propseq_runner(:,19:31,:)),4);
seqCellsByTime_sb_propseq_switch_nonrunner = nanmean(cat(4,seqCellsByTime_sb_propseq_nonrunner(:,6:18,:),seqCellsByTime_sb_propseq_nonrunner(:,19:31,:)),4);


%% Get sequence warp

% for each 100t block across all 5 days
warp_exp1_a = nan(d_info.numAnimals,numBlocks); % [animal, block]
warp_exp1_b = nan(d_info.numAnimals,numBlocks); % [animal, block]
numSeqCells = nan(d_info.numAnimals,numBlocks); % [animal, block]
propSeqCells = nan(d_info.numAnimals,numBlocks); % [animal, block]
for i=1:d_info.numAnimals
    for j=1:d_info.numDays
        if j==1
            temp0=0;
        else
            temp0 = cumsum(blocks);
            temp0 = temp0(j-1);
        end
        if isfield(d{i,j},'warp_100t')
            for k=1:length(d{i,j}.warp_100t)
                this_block = temp0+k;
                warp_exp1_a(i,this_block) = d{i,j}.warp_100t{k}.peak.exp1.a;
                warp_exp1_b(i,this_block) = d{i,j}.warp_100t{k}.peak.exp1.b;
                try
                    this_numCells = length(find(d{i,j}.tng_100t{k}.prop.iscell==1));
                    numSeqCells(i,this_block) = d{i,j}.warp_100t{k}.input.numSequenceCells;
                    propSeqCells(i,this_block) = numSeqCells(i,this_block) / this_numCells;
                catch
                end
            end
        end
    end
end

% numSeqCells split by running
numSeqCells_runner = nan(d_info.numAnimals,numBlocks);
numSeqCells_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            numSeqCells_runner(i,k) = numSeqCells(i,k);
        elseif this_running(i,k)==0
            numSeqCells_nonrunner(i,k) = numSeqCells(i,k);
        end
    end
end
numSeqCells_switch_runner = nanmean(cat(3,numSeqCells_runner(:,6:18),numSeqCells_runner(:,19:31)),3);
numSeqCells_switch_nonrunner = nanmean(cat(3,numSeqCells_nonrunner(:,6:18),numSeqCells_nonrunner(:,19:31)),3);

% propSeqCells split by running
propSeqCells_runner = nan(d_info.numAnimals,numBlocks);
propSeqCells_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            propSeqCells_runner(i,k) = propSeqCells(i,k);
        elseif this_running(i,k)==0
            propSeqCells_nonrunner(i,k) = propSeqCells(i,k);
        end
    end
end
propSeqCells_switch_runner = nanmean(cat(3,propSeqCells_runner(:,6:18),propSeqCells_runner(:,19:31)),3);
propSeqCells_switch_nonrunner = nanmean(cat(3,propSeqCells_nonrunner(:,6:18),propSeqCells_nonrunner(:,19:31)),3);

% warp_exp1_a split by running
warp_exp1_a_runner = nan(d_info.numAnimals,numBlocks);
warp_exp1_a_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            warp_exp1_a_runner(i,k) = warp_exp1_a(i,k);
        elseif this_running(i,k)==0
            warp_exp1_a_nonrunner(i,k) = warp_exp1_a(i,k);
        end
    end
end
warp_exp1_a_switch_runner = nanmean(cat(3,warp_exp1_a_runner(:,6:18),warp_exp1_a_runner(:,19:31)),3);
warp_exp1_a_switch_nonrunner = nanmean(cat(3,warp_exp1_a_nonrunner(:,6:18),warp_exp1_a_nonrunner(:,19:31)),3);

% warp_exp1_b split by running
warp_exp1_b_runner = nan(d_info.numAnimals,numBlocks);
warp_exp1_b_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            warp_exp1_b_runner(i,k) = warp_exp1_b(i,k);
        elseif this_running(i,k)==0
            warp_exp1_b_nonrunner(i,k) = warp_exp1_b(i,k);
        end
    end
end
warp_exp1_b_switch_runner = nanmean(cat(3,warp_exp1_b_runner(:,6:18),warp_exp1_b_runner(:,19:31)),3);
warp_exp1_b_switch_nonrunner = nanmean(cat(3,warp_exp1_b_nonrunner(:,6:18),warp_exp1_b_nonrunner(:,19:31)),3);


%% Calculate running-sequence correlations

corr_seq_v = nan(d_info.numAnimals,numBlocks);
corr_seq_a = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        corr_seq_v(i,k) = corr(seqCellsByTime_sb_propseq(:,k,i),velocity(p.general.bins_analysisWindow(1:end-1),k,i),'Type','Pearson','Rows','Complete');
        corr_seq_a(i,k) = corr(seqCellsByTime_sb_propseq(:,k,i),acceleration(p.general.bins_analysisWindow(1:end-1),k,i),'Type','Pearson','Rows','Complete');
    end
end

% corr_seq_v split by running
corr_seq_v_runner = nan(d_info.numAnimals,numBlocks);
corr_seq_v_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            corr_seq_v_runner(i,k) = corr_seq_v(i,k);
        elseif this_running(i,k)==0
            corr_seq_v_nonrunner(i,k) = corr_seq_v(i,k);
        end
    end
end
corr_seq_v_switch_runner = nanmean(cat(3,corr_seq_v_runner(:,6:18),corr_seq_v_runner(:,19:31)),3);
corr_seq_v_switch_nonrunner = nanmean(cat(3,corr_seq_v_nonrunner(:,6:18),corr_seq_v_nonrunner(:,19:31)),3);

% corr_seq_a split by running
corr_seq_a_runner = nan(d_info.numAnimals,numBlocks);
corr_seq_a_nonrunner = nan(d_info.numAnimals,numBlocks);
for i=1:d_info.numAnimals
    for k=1:numBlocks
        if this_running(i,k)==1
            corr_seq_a_runner(i,k) = corr_seq_a(i,k);
        elseif this_running(i,k)==0
            corr_seq_a_nonrunner(i,k) = corr_seq_a(i,k);
        end
    end
end
corr_seq_a_switch_runner = nanmean(cat(3,corr_seq_a_runner(:,6:18),corr_seq_a_runner(:,19:31)),3);
corr_seq_a_switch_nonrunner = nanmean(cat(3,corr_seq_a_nonrunner(:,6:18),corr_seq_a_nonrunner(:,19:31)),3);


%% Fig3_NumSeqCells

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;

shadedErrorBar(1:numBlocks_switch,nanmean(numSeqCells_switch_nonrunner,1),nansem(numSeqCells_switch_nonrunner,1),'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(numSeqCells_switch_runner,1),nansem(numSeqCells_switch_runner,1),'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
ylim([0,200])
yticks([0:100:200])
ylabel({'Number of Sequence Cells'})

savefig(F,[save_root_fig,'\Fig3_NumSeqCells.fig']);
saveas(F,[save_root_png,'\Fig3_NumSeqCells.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_NumSeqCells.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = numSeqCells_switch_nonrunner;
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = numSeqCells_switch_runner;
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_PropSeqCells

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;

shadedErrorBar(1:numBlocks_switch,nanmean(propSeqCells_switch_nonrunner,1)*100,nansem(propSeqCells_switch_nonrunner,1)*100,'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(propSeqCells_switch_runner,1)*100,nansem(propSeqCells_switch_runner,1)*100,'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
ytickformat('percentage')
ylim([0,30])
yticks([0:10:30])
ylabel({'Sequence participation','probability'})

savefig(F,[save_root_fig,'\Fig3_PropSeqCells.fig']);
saveas(F,[save_root_png,'\Fig3_PropSeqCells.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_PropSeqCells.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = propSeqCells_switch_nonrunner;
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = propSeqCells_switch_runner;
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_SequenceCells_perIscell

F = paper_figure([0,0.5,mm2inch(6.3*34),mm2inch(0.7*34)]); hold on;

for m=1:numBlocks_switch
    subplot(1,13,m); hold on;

    shadedErrorBar(smallBinEdges(1:end-1),nanmean(seqCellsByTime_sb_prop_switch_nonrunner(:,m,:),3),nansem(seqCellsByTime_sb_prop_switch_nonrunner(:,m,:),3),'lineProps',p.col.nonrunner);
    shadedErrorBar(smallBinEdges(1:end-1),nanmean(seqCellsByTime_sb_prop_switch_runner(:,m,:),3),nansem(seqCellsByTime_sb_prop_switch_runner(:,m,:),3),'lineProps',p.col.runner);
 
    xlim([0,5]); xticks([0,5]); xticklabels({});
    ylim([0,0.1]); yticks([0,0.1]); yticklabels({});
end

savefig(F,[save_root_fig,'\Fig3_SequenceCells_perIscell.fig']);
saveas(F,[save_root_png,'\Fig3_SequenceCells_perIscell.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_SequenceCells_perIscell.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_SequenceCells_perSeqCell

F = paper_figure([0,0.5,mm2inch(6.3*34),mm2inch(0.7*34)]); hold on;

for m=1:numBlocks_switch
    subplot(1,13,m); hold on;

    shadedErrorBar(smallBinEdges(1:end-1),nanmean(seqCellsByTime_sb_propseq_switch_nonrunner(:,m,:),3),nansem(seqCellsByTime_sb_propseq_switch_nonrunner(:,m,:),3),'lineProps',p.col.nonrunner);
    shadedErrorBar(smallBinEdges(1:end-1),nanmean(seqCellsByTime_sb_propseq_switch_runner(:,m,:),3),nansem(seqCellsByTime_sb_propseq_switch_runner(:,m,:),3),'lineProps',p.col.runner);
 
    xlim([0,5]); xticks([0,5]); xticklabels({});
    ylim([0,0.1]); yticks([0,0.1]); yticklabels({});   
end

savefig(F,[save_root_fig,'\Fig3_SequenceCells_perSeqCell.fig']);
saveas(F,[save_root_png,'\Fig3_SequenceCells_perSeqCell.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_SequenceCells_perSeqCell.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_SequenceFit

F = paper_figure([0,0.5,mm2inch(6.3*34),mm2inch(0.7*34)]); hold on;

this_x = 0:0.01:5;%.3;
for m=1:numBlocks_switch
    subplot(1,13,m); hold on;

    these_data = [];
    for idx=1:d_info.numAnimals
        if ~isnan(warp_exp1_b_switch_nonrunner(idx,1))
            these_blocks_1 = 6:18; these_blocks_2 = 19:31;
            this_a = nanmean([warp_exp1_a_nonrunner(idx,these_blocks_1(m)),warp_exp1_a_nonrunner(idx,these_blocks_2(m))]);
            this_b = nanmean([warp_exp1_b_nonrunner(idx,these_blocks_1(m)),warp_exp1_b_nonrunner(idx,these_blocks_2(m))]);
            this_ypred = this_a*exp(this_b*this_x);
            these_data = [these_data; this_ypred];
        end
    end
    shadedErrorBar(this_x,nanmean(these_data,1),nansem(these_data,1),'lineProps',p.col.nonrunner);
    
    these_data = [];
    for idx=1:d_info.numAnimals
        if ~isnan(warp_exp1_b_switch_runner(idx,1))
            these_blocks_1 = 6:18; these_blocks_2 = 19:31;
            this_a = nanmean([warp_exp1_a_runner(idx,these_blocks_1(m)),warp_exp1_a_runner(idx,these_blocks_2(m))]);
            this_b = nanmean([warp_exp1_b_runner(idx,these_blocks_1(m)),warp_exp1_b_runner(idx,these_blocks_2(m))]);
            this_ypred = this_a*exp(this_b*this_x);
            these_data = [these_data; this_ypred];
        end
    end
    shadedErrorBar(this_x,nanmean(these_data,1),nansem(these_data,1),'lineProps',p.col.runner);
    
    xlim([0,5]); xticks([0,5]); xticklabels({});
    ylim([0,0.1]); yticks([0,0.1]); yticklabels({});
end

savefig(F,[save_root_fig,'\Fig3_SequenceFit.fig']);
saveas(F,[save_root_png,'\Fig3_SequenceFit.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_SequenceFit.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_SequenceShape_a

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;

shadedErrorBar(1:numBlocks_switch,nanmean(warp_exp1_a_switch_nonrunner,1),nansem(warp_exp1_a_switch_nonrunner,1),'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(warp_exp1_a_switch_runner,1),nansem(warp_exp1_a_switch_runner,1),'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
ylim([0,0.2])
yticks([0:0.1:0.2])
ylabel({'Sequence shape','(initial proportion)'})

savefig(F,[save_root_fig,'\Fig3_SequenceShape_a.fig']);
saveas(F,[save_root_png,'\Fig3_SequenceShape_a.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_SequenceShape_a.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = warp_exp1_a_switch_nonrunner;
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = warp_exp1_a_switch_runner;
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_SequenceShape_b

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on; % F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;

%yline(0,'k-');
shadedErrorBar(1:numBlocks_switch,nanmean(warp_exp1_b_switch_nonrunner,1),nansem(warp_exp1_b_switch_nonrunner,1),'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(warp_exp1_b_switch_runner,1),nansem(warp_exp1_b_switch_runner,1),'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
ylim([-1,0])
yticks([-1:0.5:0])
ylabel({'Sequence shape','(decay rate, s-1)'})

savefig(F,[save_root_fig,'\Fig3_SequenceShape_b.fig']);
saveas(F,[save_root_png,'\Fig3_SequenceShape_b.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_SequenceShape_b.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = warp_exp1_b_switch_nonrunner;
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = warp_exp1_b_switch_runner;
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_SequenceCells_early

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

%yline(0,'k-');
shadedErrorBar(1:numBlocks_switch,nanmean(seqCellsByTime_prop_switch_nonrunner(1,:,:),3)*100,nansem(seqCellsByTime_prop_switch_nonrunner(1,:,:),3)*100,'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(seqCellsByTime_prop_switch_runner(1,:,:),3)*100,nansem(seqCellsByTime_prop_switch_runner(1,:,:),3)*100,'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
ytickformat('percentage')
ylim([0,20])
yticks([0:10:20])
ylabel({'Proportion of cells'})

savefig(F,[save_root_fig,'\Fig3_SequenceCells_early.fig']);
saveas(F,[save_root_png,'\Fig3_SequenceCells_early.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_SequenceCells_early.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(seqCellsByTime_prop_switch_nonrunner(1,:,:))';
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(seqCellsByTime_prop_switch_runner(1,:,:))';
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_SequenceCells_middle

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

shadedErrorBar(1:numBlocks_switch,nanmean(seqCellsByTime_prop_switch_nonrunner(2,:,:),3)*100,nansem(seqCellsByTime_prop_switch_nonrunner(2,:,:),3)*100,'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(seqCellsByTime_prop_switch_runner(2,:,:),3)*100,nansem(seqCellsByTime_prop_switch_runner(2,:,:),3)*100,'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
ytickformat('percentage')
ylim([0,15])
yticks([0:5:15])
ylabel({'Proportion of cells'})

savefig(F,[save_root_fig,'\Fig3_SequenceCells_middle.fig']);
saveas(F,[save_root_png,'\Fig3_SequenceCells_middle.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_SequenceCells_middle.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(seqCellsByTime_prop_switch_nonrunner(2,:,:))';
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(seqCellsByTime_prop_switch_runner(2,:,:))';
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_SequenceCells_late

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

shadedErrorBar(1:numBlocks_switch,nanmean(seqCellsByTime_prop_switch_nonrunner(3,:,:),3)*100,nansem(seqCellsByTime_prop_switch_nonrunner(3,:,:),3)*100,'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(seqCellsByTime_prop_switch_runner(3,:,:),3)*100,nansem(seqCellsByTime_prop_switch_runner(3,:,:),3)*100,'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
ytickformat('percentage')
ylim([0,8])
yticks([0:4:8])
ylabel({'Proportion of cells'})

savefig(F,[save_root_fig,'\Fig3_SequenceCells_late.fig']);
saveas(F,[save_root_png,'\Fig3_SequenceCells_late.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_SequenceCells_late.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(seqCellsByTime_prop_switch_nonrunner(3,:,:))';
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(seqCellsByTime_prop_switch_runner(3,:,:))';
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')



%% Fig3_SequenceCells_early_abs

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

%yline(0,'k-');
shadedErrorBar(1:numBlocks_switch,nanmean(seqCellsByTime_num_switch_nonrunner(1,:,:),3),nansem(seqCellsByTime_num_switch_nonrunner(1,:,:),3),'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(seqCellsByTime_num_switch_runner(1,:,:),3),nansem(seqCellsByTime_num_switch_runner(1,:,:),3),'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
% ytickformat('percentage')
% ylim([0,20])
% yticks([0:10:20])
ylabel({'Number of cells'})

savefig(F,[save_root_fig,'\Fig3_SequenceCells_early_abs.fig']);
saveas(F,[save_root_png,'\Fig3_SequenceCells_early_abs.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_SequenceCells_early_abs.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(seqCellsByTime_num_switch_nonrunner(1,:,:))';
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(seqCellsByTime_num_switch_runner(1,:,:))';
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_SequenceCells_middle_abs

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

shadedErrorBar(1:numBlocks_switch,nanmean(seqCellsByTime_num_switch_nonrunner(2,:,:),3),nansem(seqCellsByTime_num_switch_nonrunner(2,:,:),3),'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(seqCellsByTime_num_switch_runner(2,:,:),3),nansem(seqCellsByTime_num_switch_runner(2,:,:),3),'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
% ytickformat('percentage')
% ylim([0,15])
% yticks([0:5:15])
ylabel({'Number of cells'})

savefig(F,[save_root_fig,'\Fig3_SequenceCells_middle_abs.fig']);
saveas(F,[save_root_png,'\Fig3_SequenceCells_middle_abs.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_SequenceCells_middle_abs.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(seqCellsByTime_num_switch_nonrunner(2,:,:))';
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(seqCellsByTime_num_switch_runner(2,:,:))';
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_SequenceCells_late_abs

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

shadedErrorBar(1:numBlocks_switch,nanmean(seqCellsByTime_num_switch_nonrunner(3,:,:),3),nansem(seqCellsByTime_num_switch_nonrunner(3,:,:),3),'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(seqCellsByTime_num_switch_runner(3,:,:),3),nansem(seqCellsByTime_num_switch_runner(3,:,:),3),'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
% ytickformat('percentage')
% ylim([0,8])
% yticks([0:4:8])
ylabel({'Number of cells'})

savefig(F,[save_root_fig,'\Fig3_SequenceCells_late_abs.fig']);
saveas(F,[save_root_png,'\Fig3_SequenceCells_late_abs.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_SequenceCells_late_abs.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(seqCellsByTime_num_switch_nonrunner(3,:,:))';
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(seqCellsByTime_num_switch_runner(3,:,:))';
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_Bcon_Distance_Runners

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;
these_rgbs = discretisedColourMap('winter',false,13);

hold on;
taskLines(p,info);
yline(0,'k-');
for i=1:8
    plot(1:p.general.numBins,nanmean(distance_switch_runner(:,i,:),3),'Color',these_rgbs(i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(distance_switch_runner(:,i,:),3),nansem(distance_switch_runner(:,i,:),3),'lineProps',these_rgbs(i,:))
end
for i=1:5
    plot(1:p.general.numBins,nanmean(distance_switch_runner(:,i,:),3),'Color',these_rgbs(8+i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(distance_switch_runner(:,i,:),3),nansem(distance_switch_runner(:,i,:),3),'lineProps',these_rgbs(8+i,:))
end
ylim([-400,800])
yticks([-400:400:800])
ylabel('Distance traveled (cm)')

savefig(F,[save_root_fig,'\Fig3_Bcon_Distance_Runners.fig']);
saveas(F,[save_root_png,'\Fig3_Bcon_Distance_Runners.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_Bcon_Distance_Runners.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_Bcon_Distance_Nonrunners

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;
these_rgbs = discretisedColourMap('winter',false,13);

hold on;
taskLines(p,info);
yline(0,'k-');
for i=1:8
    plot(1:p.general.numBins,nanmean(distance_switch_nonrunner(:,i,:),3),'Color',these_rgbs(i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(distance_switch_nonrunner(:,i,:),3),nansem(distance_switch_nonrunner(:,i,:),3),'lineProps',these_rgbs(i,:))
end
for i=1:5
    plot(1:p.general.numBins,nanmean(distance_switch_nonrunner(:,i,:),3),'Color',these_rgbs(8+i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(distance_switch_nonrunner(:,i,:),3),nansem(distance_switch_nonrunner(:,i,:),3),'lineProps',these_rgbs(8+i,:))
end
ylim([-400,800])
yticks([-400:400:800])
ylabel('Distance traveled (cm)')

savefig(F,[save_root_fig,'\Fig3_Bcon_Distance_Nonrunners.fig']);
saveas(F,[save_root_png,'\Fig3_Bcon_Distance_Nonrunners.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_Bcon_Distance_Nonrunners.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_Bcon_Velocity_Runners

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;
these_rgbs = discretisedColourMap('winter',false,13);

hold on;
taskLines(p,info);
for i=1:8
    plot(1:p.general.numBins,nanmean(velocity_switch_runner(:,i,:),3),'Color',these_rgbs(i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(velocity_switch_runner(:,i,:),3),nansem(velocity_switch_runner(:,i,:),3),'lineProps',these_rgbs(i,:))
end
for i=1:5
    plot(1:p.general.numBins,nanmean(velocity_switch_runner(:,i,:),3),'Color',these_rgbs(8+i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(velocity_switch_runner(:,i,:),3),nansem(velocity_switch_runner(:,i,:),3),'lineProps',these_rgbs(8+i,:))
end
ylim([0,120])
yticks([0:40:120])
ylabel('Velocity (cm/s)')

savefig(F,[save_root_fig,'\Fig3_Bcon_Velocity_Runners.fig']);
saveas(F,[save_root_png,'\Fig3_Bcon_Velocity_Runners.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_Bcon_Velocity_Runners.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_Bcon_Velocity_Nonrunners

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;
these_rgbs = discretisedColourMap('winter',false,13);

hold on;
taskLines(p,info);
for i=1:8
    plot(1:p.general.numBins,nanmean(velocity_switch_nonrunner(:,i,:),3),'Color',these_rgbs(8+i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(velocity_switch_nonrunner(:,i,:),3),nansem(velocity_switch_nonrunner(:,i,:),3),'lineProps',these_rgbs(8+i,:))
end
for i=1:5
    plot(1:p.general.numBins,nanmean(velocity_switch_nonrunner(:,i,:),3),'Color',these_rgbs(8+i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(velocity_switch_nonrunner(:,i,:),3),nansem(velocity_switch_nonrunner(:,i,:),3),'lineProps',these_rgbs(8+i,:))
end
ylim([0,120])
yticks([0:40:120])
ylabel('Velocity (cm/s)')

savefig(F,[save_root_fig,'\Fig3_Bcon_Velocity_Nonrunners.fig']);
saveas(F,[save_root_png,'\Fig3_Bcon_Velocity_Nonrunners.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_Bcon_Velocity_Nonrunners.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_Bcon_Acceleration_Runners

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;
these_rgbs = discretisedColourMap('winter',false,13);

hold on;
taskLines(p,info);
yline(0,'k-');
for i=1:8
    plot(1:p.general.numBins,nanmean(acceleration_switch_runner(:,i,:),3),'Color',these_rgbs(i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(acceleration_switch_runner(:,i,:),3),nansem(acceleration_switch_runner(:,i,:),3),'lineProps',these_rgbs(i,:))
end
for i=1:5
    plot(1:p.general.numBins,nanmean(acceleration_switch_runner(:,i,:),3),'Color',these_rgbs(8+i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(acceleration_switch_runner(:,i,:),3),nansem(acceleration_switch_runner(:,i,:),3),'lineProps',these_rgbs(8+i,:))
end
ylim([-80,80])
yticks([-80:40:80])
ylabel('Acceleration (cm/s2)')

savefig(F,[save_root_fig,'\Fig3_Bcon_Acceleration_Runners.fig']);
saveas(F,[save_root_png,'\Fig3_Bcon_Acceleration_Runners.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_Bcon_Acceleration_Runners.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_Bcon_Acceleration_Nonrunners

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;
these_rgbs = discretisedColourMap('winter',false,13);

hold on;
taskLines(p,info);
yline(0,'k-');
for i=1:8
    plot(1:p.general.numBins,nanmean(acceleration_switch_nonrunner(:,i,:),3),'Color',these_rgbs(i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(acceleration_switch_nonrunner(:,i,:),3),nansem(acceleration_switch_nonrunner(:,i,:),3),'lineProps',these_rgbs(i,:))
end
for i=1:5
    plot(1:p.general.numBins,nanmean(acceleration_switch_nonrunner(:,i,:),3),'Color',these_rgbs(8+i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(acceleration_switch_nonrunner(:,i,:),3),nansem(acceleration_switch_nonrunner(:,i,:),3),'lineProps',these_rgbs(8+i,:))
end
ylim([-80,80])
yticks([-80:40:80])
ylabel('Acceleration (cm/s2)')

savefig(F,[save_root_fig,'\Fig3_Bcon_Acceleration_Nonrunners.fig']);
saveas(F,[save_root_png,'\Fig3_Bcon_Acceleration_Nonrunners.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_Bcon_Acceleration_Nonrunners.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_Bcon_Licking_Runners

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;
these_rgbs = discretisedColourMap('winter',false,13);

hold on;
taskLines(p,info);
for i=1:8
    plot(1:p.general.numBins,nanmean(licking_switch_runner(:,i,:),3)*100,'Color',these_rgbs(i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(licking_switch_runner(:,i,:),3)*100,nansem(licking_switch_runner(:,i,:),3)*100,'lineProps',these_rgbs(i,:))
end
for i=1:5
    plot(1:p.general.numBins,nanmean(licking_switch_runner(:,i,:),3)*100,'Color',these_rgbs(8+i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(licking_switch_runner(:,i,:),3)*100,nansem(licking_switch_runner(:,i,:),3)*100,'lineProps',these_rgbs(8+i,:))
end
ytickformat('percentage')
ylim([0,100])
yticks([0,50,100])
ylabel('Lick probability')

savefig(F,[save_root_fig,'\Fig3_Bcon_Licking_Runners.fig']);
saveas(F,[save_root_png,'\Fig3_Bcon_Licking_Runners.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_Bcon_Licking_Runners.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_Bcon_Licking_Nonrunners

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;
these_rgbs = discretisedColourMap('winter',false,13);

hold on;
taskLines(p,info);
for i=1:8
    plot(1:p.general.numBins,nanmean(licking_switch_nonrunner(:,i,:),3)*100,'Color',these_rgbs(i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(licking_switch_nonrunner(:,i,:),3)*100,nansem(licking_switch_nonrunner(:,i,:),3)*100,'lineProps',these_rgbs(i,:))
end
for i=1:5
    plot(1:p.general.numBins,nanmean(licking_switch_nonrunner(:,i,:),3)*100,'Color',these_rgbs(8+i,:))
    %shadedErrorBar(1:p.general.numBins,nanmean(licking_switch_nonrunner(:,i,:),3)*100,nansem(licking_switch_nonrunner(:,i,:),3)*100,'lineProps',these_rgbs(8+i,:))
end
ytickformat('percentage')
ylim([0,100])
yticks([0,50,100])
ylabel('Lick probability')

savefig(F,[save_root_fig,'\Fig3_Bcon_Licking_Nonrunners.fig']);
saveas(F,[save_root_png,'\Fig3_Bcon_Licking_Nonrunners.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_Bcon_Licking_Nonrunners.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_PerformanceCorr_SequenceShape

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_x_nonrunners = warp_exp1_b_switch_nonrunner;
this_data_y_nonrunners = correct_switch_nonrunner;
this_data_x_runners = warp_exp1_b_switch_runner;
this_data_y_runners = correct_switch_runner;
this_data_x = [this_data_x_nonrunners;this_data_x_runners];
this_data_y = [this_data_y_nonrunners;this_data_y_runners];

yline(50,'k:');
scatter(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.nonrunner)
scatter(this_data_x_runners(:),this_data_y_runners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k')

xlim([-2,0.5])
xticks([-2,0,0.5])
xlabel({'Sequence shape','(decay rate, s-1)'})
ytickformat('percentage')
ylim([35,100]) %ylim([0,100])
yticks([50,100]) %yticks([0,50,100])
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig3_PerformanceCorr_SequenceShape.fig']);
saveas(F,[save_root_png,'\Fig3_PerformanceCorr_SequenceShape.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_PerformanceCorr_SequenceShape.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_PerformanceCorr_Velocity

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_x_nonrunners = squeeze(nanmean(velocity_switch_nonrunner(p.general.bins_analysisWindow,:,:),1))';
this_data_y_nonrunners = correct_switch_nonrunner;
this_data_x_runners = squeeze(nanmean(velocity_switch_runner(p.general.bins_analysisWindow,:,:),1))';
this_data_y_runners = correct_switch_runner;
this_data_x = [this_data_x_nonrunners;this_data_x_runners];
this_data_y = [this_data_y_nonrunners;this_data_y_runners];

yline(50,'k:');
scatter(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.nonrunner)
scatter(this_data_x_runners(:),this_data_y_runners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k')

xlim([0,80])
xticks([0:40:80])
xlabel({'Velocity','(cm/s)'})
ytickformat('percentage')
ylim([35,100])
yticks([50,100])
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig3_PerformanceCorr_Velocity.fig']);
saveas(F,[save_root_png,'\Fig3_PerformanceCorr_Velocity.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_PerformanceCorr_Velocity.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_PerformanceCorr_Acceleration

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_x_nonrunners = squeeze(nanmean(acceleration_switch_nonrunner(p.general.bins_analysisWindow,:,:),1))';
this_data_y_nonrunners = correct_switch_nonrunner;
this_data_x_runners = squeeze(nanmean(acceleration_switch_runner(p.general.bins_analysisWindow,:,:),1))';
this_data_y_runners = correct_switch_runner;
this_data_x = [this_data_x_nonrunners;this_data_x_runners];
this_data_y = [this_data_y_nonrunners;this_data_y_runners];

yline(50,'k:');
scatter(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.nonrunner)
scatter(this_data_x_runners(:),this_data_y_runners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k')

xlim([0,3])
xticks([0:1.5:3])
xlabel({'Acceleration','(cm/s2)'})
ytickformat('percentage')
ylim([35,100])
yticks([50,100])
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig3_PerformanceCorr_Acceleration.fig']);
saveas(F,[save_root_png,'\Fig3_PerformanceCorr_Acceleration.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_PerformanceCorr_Acceleration.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_PerformanceCorr_CorrV

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_x_nonrunners = corr_seq_v_switch_nonrunner;
this_data_y_nonrunners = correct_switch_nonrunner;
this_data_x_runners = corr_seq_v_switch_runner;
this_data_y_runners = correct_switch_runner;
this_data_x = [this_data_x_nonrunners;this_data_x_runners];
this_data_y = [this_data_y_nonrunners;this_data_y_runners];

yline(50,'k:');
scatter(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.nonrunner)
scatter(this_data_x_runners(:),this_data_y_runners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k')

xlim([-1,1])
xticks([-1,0,1])
xlabel({'Correlation','(peak times x velocity)'})
ytickformat('percentage')
ylim([35,100])
yticks([50,100])
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig3_PerformanceCorr_CorrV.fig']);
saveas(F,[save_root_png,'\Fig3_PerformanceCorr_CorrV.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_PerformanceCorr_CorrV.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_PerformanceCorr_CorrA

F = paper_figure([0,0.5,mm2inch(34),mm2inch(34)]); hold on;

this_data_x_nonrunners = corr_seq_a_switch_nonrunner;
this_data_y_nonrunners = correct_switch_nonrunner;
this_data_x_runners = corr_seq_a_switch_runner;
this_data_y_runners = correct_switch_runner;
this_data_x = [this_data_x_nonrunners;this_data_x_runners];
this_data_y = [this_data_y_nonrunners;this_data_y_runners];

yline(50,'k:');
scatter(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.nonrunner)
scatter(this_data_x_runners(:),this_data_y_runners(:)*100,'.','SizeData',30,'MarkerEdgeColor',p.col.runner)
[this_corr_r_nonrunners,this_corr_p_nonrunners] = fitLine(this_data_x_nonrunners(:),this_data_y_nonrunners(:)*100,p.col.nonrunner);
[this_corr_r_runners,this_corr_p_runners] = fitLine(this_data_x_runners(:),this_data_y_runners(:)*100,p.col.runner);
[this_corr_r,this_corr_p] = fitLine(this_data_x(:),this_data_y(:)*100,'k')

xlim([-1,1])
xticks([-1,0,1])
xlabel({'Correlation','(peak times x acceleration)'})
ytickformat('percentage')
ylim([35,100])
yticks([50,100])
ylabel('Performance')

savefig(F,[save_root_fig,'\Fig3_PerformanceCorr_CorrA.fig']);
saveas(F,[save_root_png,'\Fig3_PerformanceCorr_CorrA.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_PerformanceCorr_CorrA.pdf']); set(gcf,'Color',[1,1,1])


%% Fig3_CorrHistograms_full

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;

subplot(2,2,1); hold on;
histogram(corr_seq_v_switch_nonrunner(:),-1:0.33:1)
xline(0,'r');
xlabel('non-runners: corr(seq,v)')
[temp1,~,temp2] = signrank(corr_seq_v_switch_nonrunner(:))

subplot(2,2,2); hold on;
histogram(corr_seq_v_switch_runner(:),-1:0.33:1)
xline(0,'r');
xlabel('runners: corr(seq,v)')
[temp1,~,temp2] = signrank(corr_seq_v_switch_runner(:))

subplot(2,2,3); hold on;
histogram(corr_seq_a_switch_nonrunner(:),-1:0.33:1)
xline(0,'r');
xlabel('non-runners: corr(seq,a)')
[temp1,~,temp2] = signrank(corr_seq_a_switch_nonrunner(:))

subplot(2,2,4); hold on;
histogram(corr_seq_a_switch_runner(:),-1:0.33:1)
xline(0,'r');
xlabel('runners: corr(seq,a)')
[temp1,~,temp2] = signrank(corr_seq_a_switch_runner(:))


%% Fig3_CorrHistograms

F = paper_figure([0,0.5,mm2inch(1.25*34),mm2inch(34)]); hold on;

subplot(2,1,1); hold on;
histogram(rmmissing(corr_seq_v_switch_runner(:)),-1:0.2:1,'Normalization','probability','FaceColor',p.col.runner,'FaceAlpha',1,'EdgeColor','none');
xline(0,'k:','Color',p.col.black,'LineWidth',1,'Alpha',1);
xline(nanmean(corr_seq_v_switch_runner(:)),'k-','Color',p.col.black,'LineWidth',1,'Alpha',1);
xlim([-1,1])
xticks([-1:1:1])
xlabel('Correlation of firing field peak and velocity')
ylim([0,0.3])
yticks([0,0.3])
yticklabels({'0%','30%'})
ylabel('Proportion')


subplot(2,1,2); hold on;
histogram(rmmissing(corr_seq_a_switch_runner(:)),-1:0.2:1,'Normalization','probability','FaceColor',p.col.runner,'FaceAlpha',1,'EdgeColor','none');
xline(0,'k:','Color',p.col.black,'LineWidth',1,'Alpha',1);
xline(nanmean(corr_seq_a_switch_runner(:)),'k-','Color',p.col.black,'LineWidth',1,'Alpha',1);
xlim([-1,1])
xticks([-1:1:1])
xlabel('Correlation of firing field peak and acceleration')
ylim([0,0.3])
yticks([0,0.3])
yticklabels({'0%','30%'})
ylabel('Proportion')


savefig(F,[save_root_fig,'\Fig3_CorrHistograms.fig']);
saveas(F,[save_root_png,'\Fig3_CorrHistograms.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_CorrHistograms.pdf']); set(gcf,'Color',[1,1,1])



%% Fig3_Decoding_E

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

%yline(0,'k-');
shadedErrorBar(1:numBlocks_switch,nanmean(typeDecodingCorrect_E_switch_nonrunner,1)*100,nansem(typeDecodingCorrect_E_switch_nonrunner,1)*100,'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(typeDecodingCorrect_E_switch_runner,1)*100,nansem(typeDecodingCorrect_E_switch_runner,1)*100,'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
ytickformat('percentage')
ylim([50,100])
yticks([50,75,100])
ylabel({'1st Odor decoding accuracy'})

savefig(F,[save_root_fig,'\Fig3_Decoding_E.fig']);
saveas(F,[save_root_png,'\Fig3_Decoding_E.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_Decoding_E.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(typeDecodingCorrect_E_switch_nonrunner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(typeDecodingCorrect_E_switch_runner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_Decoding_I

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

%yline(0,'k-');
shadedErrorBar(1:numBlocks_switch,nanmean(typeDecodingCorrect_I_switch_nonrunner,1)*100,nansem(typeDecodingCorrect_I_switch_nonrunner,1)*100,'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(typeDecodingCorrect_I_switch_runner,1)*100,nansem(typeDecodingCorrect_I_switch_runner,1)*100,'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
ytickformat('percentage')
ylim([50,100])
yticks([50,75,100])
ylabel({'1st Odor decoding accuracy'})

savefig(F,[save_root_fig,'\Fig3_Decoding_I.fig']);
saveas(F,[save_root_png,'\Fig3_Decoding_I.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_Decoding_I.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(typeDecodingCorrect_I_switch_nonrunner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(typeDecodingCorrect_I_switch_runner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_Decoding_L

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

%yline(0,'k-');
shadedErrorBar(1:numBlocks_switch,nanmean(typeDecodingCorrect_L_switch_nonrunner,1)*100,nansem(typeDecodingCorrect_L_switch_nonrunner,1)*100,'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(typeDecodingCorrect_L_switch_runner,1)*100,nansem(typeDecodingCorrect_L_switch_runner,1)*100,'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
ytickformat('percentage')
ylim([50,100])
yticks([50,75,100])
ylabel({'1st Odor decoding accuracy'})

savefig(F,[save_root_fig,'\Fig3_Decoding_L.fig']);
saveas(F,[save_root_png,'\Fig3_Decoding_L.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_Decoding_L.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(typeDecodingCorrect_L_switch_nonrunner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(typeDecodingCorrect_L_switch_runner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_Decoding_E_allCells

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

%yline(0,'k-');
shadedErrorBar(1:numBlocks_switch,nanmean(typeDecodingCorrect_E_allCells_switch_nonrunner,1)*100,nansem(typeDecodingCorrect_E_allCells_switch_nonrunner,1)*100,'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(typeDecodingCorrect_E_allCells_switch_runner,1)*100,nansem(typeDecodingCorrect_E_allCells_switch_runner,1)*100,'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
ytickformat('percentage')
ylim([50,100])
yticks([50,75,100])
ylabel({'1st Odor decoding accuracy'})

savefig(F,[save_root_fig,'\Fig3_Decoding_E_allCells.fig']);
saveas(F,[save_root_png,'\Fig3_Decoding_E_allCells.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_Decoding_E_allCells.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(typeDecodingCorrect_E_allCells_switch_nonrunner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(typeDecodingCorrect_E_allCells_switch_runner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_Decoding_I_allCells

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

%yline(0,'k-');
shadedErrorBar(1:numBlocks_switch,nanmean(typeDecodingCorrect_I_allCells_switch_nonrunner,1)*100,nansem(typeDecodingCorrect_I_allCells_switch_nonrunner,1)*100,'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(typeDecodingCorrect_I_allCells_switch_runner,1)*100,nansem(typeDecodingCorrect_I_allCells_switch_runner,1)*100,'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
ytickformat('percentage')
ylim([50,100])
yticks([50,75,100])
ylabel({'1st Odor decoding accuracy'})

savefig(F,[save_root_fig,'\Fig3_Decoding_I_allCells.fig']);
saveas(F,[save_root_png,'\Fig3_Decoding_I_allCells.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_Decoding_I_allCells.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(typeDecodingCorrect_I_allCells_switch_nonrunner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(typeDecodingCorrect_I_allCells_switch_runner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_Decoding_L_allCells

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

%yline(0,'k-');
shadedErrorBar(1:numBlocks_switch,nanmean(typeDecodingCorrect_L_allCells_switch_nonrunner,1)*100,nansem(typeDecodingCorrect_L_allCells_switch_nonrunner,1)*100,'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(typeDecodingCorrect_L_allCells_switch_runner,1)*100,nansem(typeDecodingCorrect_L_allCells_switch_runner,1)*100,'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')
ytickformat('percentage')
ylim([50,100])
yticks([50,75,100])
ylabel({'1st Odor decoding accuracy'})

savefig(F,[save_root_fig,'\Fig3_Decoding_L_allCells.fig']);
saveas(F,[save_root_png,'\Fig3_Decoding_L_allCells.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_Decoding_L_allCells.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(typeDecodingCorrect_L_allCells_switch_nonrunner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(typeDecodingCorrect_L_allCells_switch_runner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')




%% Fig3_DecodingT_E

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

%yline(0,'k-');
shadedErrorBar(1:numBlocks_switch,nanmean(timeDecodingError_s_E_switch_nonrunner,1),nansem(timeDecodingError_s_E_switch_nonrunner,1),'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(timeDecodingError_s_E_switch_runner,1),nansem(timeDecodingError_s_E_switch_runner,1),'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')

ylim([0,4]) 

ylabel({'Time decoding error (s)'})

savefig(F,[save_root_fig,'\Fig3_DecodingT_E.fig']);
saveas(F,[save_root_png,'\Fig3_DecodingT_E.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_DecodingT_E.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(timeDecodingError_s_E_switch_nonrunner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(timeDecodingError_s_E_switch_runner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_DecodingT_I

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

%yline(0,'k-');
shadedErrorBar(1:numBlocks_switch,nanmean(timeDecodingError_s_I_switch_nonrunner,1),nansem(timeDecodingError_s_I_switch_nonrunner,1),'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(timeDecodingError_s_I_switch_runner,1),nansem(timeDecodingError_s_I_switch_runner,1),'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')

ylim([0,4]) 

ylabel({'Time decoding error (s)'})

savefig(F,[save_root_fig,'\Fig3_DecodingT_I.fig']);
saveas(F,[save_root_png,'\Fig3_DecodingT_I.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_DecodingT_I.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(timeDecodingError_s_I_switch_nonrunner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(timeDecodingError_s_I_switch_runner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_DecodingT_L

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

%yline(0,'k-');
shadedErrorBar(1:numBlocks_switch,nanmean(timeDecodingError_s_L_switch_nonrunner,1),nansem(timeDecodingError_s_L_switch_nonrunner,1),'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(timeDecodingError_s_L_switch_runner,1),nansem(timeDecodingError_s_L_switch_runner,1),'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')

ylim([0,4]) 

ylabel({'Time decoding error (s)'})

savefig(F,[save_root_fig,'\Fig3_DecodingT_L.fig']);
saveas(F,[save_root_png,'\Fig3_DecodingT_L.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_DecodingT_L.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(timeDecodingError_s_L_switch_nonrunner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(timeDecodingError_s_L_switch_runner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_DecodingT_E_allCells

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

%yline(0,'k-');
shadedErrorBar(1:numBlocks_switch,nanmean(timeDecodingError_s_E_allCells_switch_nonrunner,1),nansem(timeDecodingError_s_E_allCells_switch_nonrunner,1),'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(timeDecodingError_s_E_allCells_switch_runner,1),nansem(timeDecodingError_s_E_allCells_switch_runner,1),'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')

ylim([0,4]) 

ylabel({'Time decoding error (s)'})

savefig(F,[save_root_fig,'\Fig3_DecodingT_E_allCells.fig']);
saveas(F,[save_root_png,'\Fig3_DecodingT_E_allCells.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_DecodingT_E_allCells.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(timeDecodingError_s_E_allCells_switch_nonrunner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(timeDecodingError_s_E_allCells_switch_runner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_DecodingT_I_allCells

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

%yline(0,'k-');
shadedErrorBar(1:numBlocks_switch,nanmean(timeDecodingError_s_I_allCells_switch_nonrunner,1),nansem(timeDecodingError_s_I_allCells_switch_nonrunner,1),'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(timeDecodingError_s_I_allCells_switch_runner,1),nansem(timeDecodingError_s_I_allCells_switch_runner,1),'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')

ylim([0,4]) 

ylabel({'Time decoding error (s)'})

savefig(F,[save_root_fig,'\Fig3_DecodingT_I_allCells.fig']);
saveas(F,[save_root_png,'\Fig3_DecodingT_I_allCells.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_DecodingT_I_allCells.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(timeDecodingError_s_I_allCells_switch_nonrunner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(timeDecodingError_s_I_allCells_switch_runner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')


%% Fig3_DecodingT_L_allCells

F = paper_figure([0,0.5,mm2inch(1.05*34),mm2inch(34)]); hold on;

%yline(0,'k-');
shadedErrorBar(1:numBlocks_switch,nanmean(timeDecodingError_s_L_allCells_switch_nonrunner,1),nansem(timeDecodingError_s_L_allCells_switch_nonrunner,1),'lineProps',p.col.nonrunner)
shadedErrorBar(1:numBlocks_switch,nanmean(timeDecodingError_s_L_allCells_switch_runner,1),nansem(timeDecodingError_s_L_allCells_switch_runner,1),'lineProps',p.col.runner)

xlim([0,13]); xticks([1:13])
xticklabels({'1','','','','','','','8','','','','','13'})
xlabel('Block of 100 trials')

ylim([0,4]) 

ylabel({'Time decoding error (s)'})

savefig(F,[save_root_fig,'\Fig3_DecodingT_L_allCells.fig']);
saveas(F,[save_root_png,'\Fig3_DecodingT_L_allCells.png']);
set(gcf,'Color','none'); saveas(F,[save_root_pdf,'\Fig3_DecodingT_L_allCells.pdf']); set(gcf,'Color',[1,1,1])

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(timeDecodingError_s_L_allCells_switch_nonrunner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')

temp1 = repmat(1:numBlocks_switch,d_info.numAnimals,1);
temp2 = squeeze(timeDecodingError_s_L_allCells_switch_runner);
[this_corr_r,this_corr_p] = corr(temp1(:),temp2(:),'Type','Pearson','Rows','Complete')












