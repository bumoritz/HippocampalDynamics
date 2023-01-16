%% Select data and parameters

clc;
clear;
tic;

p.stimType              = 'ExpertStim1';
info.animal             = 'Python'; %'Mufasa';
info.date               = '20211210'; %'20200520';
info.recording          = 'Python-20211210-base1'; %'Mufasa-20200520-5min';
info.rootPath           = 'N:\Data\SniffinHippo\'; % 'N:\Data\SniffinHippo\';
info.calibrationFile    = 'N:\Data\SniffinHippo\PowerCalibration_old.csv'; % 'N:\Data\SniffinHippo\';



%%% --------------------------- %%%

% photostimulation power
p.powerPerCell_1cell    = 2.5; % [mW] % lowest calibration point has to be below that
p.powerPerCell_2cells   = 2; % [mW]
p.powerPerCell          = 1.5; % [mW]

% stimulation cluster assignment
p.discardFOVEdges       = 50; % [um]
p.minDistanceAcross     = 8; % [um] (above 10 takes long, even above 8)
p.maxClusterExtent      = 450; % [um](below 400-450 can take very long)
p.minDistanceWithin     = 20; % [um]
p.galvoDisplacement     = 50; % [um]
info.FOVsize            = 1000.78; % [um]
info.numPixels          = 512;

% experiment structure
info.trialsPerBlock                     = 20;
info.task.trialStructure.tOdour1        = 0.3; % [s]
info.task.trialStructure.tGap           = 5; % [s]
info.task.trialStructure.tOdour2        = 0.3; % [s]
info.task.trialStructure.tRespDelay     = 0.5; % [s]
info.task.trialStructure.tRespWindow    = 1; % [s]

% stimulation structure
p.numClusterSlots       = 20;
p.stimDuringOdour1      = false;

% colours
p.col.odour             = [252,193,2]/255; % yellow
p.col.reward            = [128,204,223]/255; % light blue
p.col.photostim         = [229,49,18]/255; % red
p.col.AB_rainbow        = [linspace(231,147,20);linspace(126,82,20);linspace(33,23,20)]'/255; % orange (gradient to brown)
p.col.XY_rainbow        = [linspace(192,130,20);linspace(202,120,20);linspace(50,29,20)]'/255; % light green (gradient to olive green)

% xlm info
xml.name_root           = 'SLMPattern_';
xml.shape               = 'Ellipse';
xml.pxSpacing           = '1';
xml.durationMS          = '2'; % spiral duration [ms]
xml.iterations          = '30'; % spiral iterations
xml.prePatIdleMS        = '0';
xml.postPatIdleMS       = '160'; % time between patterns [ms]
xml.preIteIdleMS        = '0';
xml.postIteIdleMS       = '1';
xml.measurePowerMW      = '0';
xml.measurePowerMWPerUM2= '0';
xml.sequenceEpochCount  = '1';
xml.roiWidthUM          = 8; % spiral size [um]
xml.roiHeightUM         = 8; % spiral size [um]

%%% --- checks and preparations ---

info.pixelSize          = info.FOVsize / info.numPixels;
xml.roiWidthPx          = num2str(xml.roiWidthUM / info.pixelSize);
xml.roiHeightPx         = num2str(xml.roiHeightUM / info.pixelSize);
if p.stimDuringOdour1
    error('stimDuringOdour1 not yet implemented.')
else
    info.stimEdges      = info.task.trialStructure.tOdour1:info.task.trialStructure.tGap/p.numClusterSlots:info.task.trialStructure.tOdour1+info.task.trialStructure.tGap;
end

disp(['Target selection for ',p.stimType,' experiment.'])


%% Load data

disp(['- Loading data.'])

info.year = info.date(1:4);
info.month = info.date(5:6);
info.day = info.date(7:8);
path.homeFolder = [info.rootPath,info.year,'\',info.year,'-',info.month,'\',info.year,'-',info.month,'-',info.day,'\',info.animal,'\','Targeting','\'];
path.s2pPath = [path.homeFolder,info.recording,'\suite2p\plane0\Fall.mat'];
path.iscellPath = [path.homeFolder,info.recording,'\suite2p\plane0\iscell.npy']; % just as a control, is not used otherwise
path.scaPath = [path.homeFolder,info.animal,'_',info.date,'_sca.mat'];
temp = dir([path.homeFolder,info.animal,'_',info.date,'_*_SEQ.txt']);
path.seqPath = [path.homeFolder,temp.name];
path.savePath = [path.homeFolder,'TargetSelection\'];

% import files
sca = load(path.scaPath);
sca = sca.sca;
iscells = sca.prop.iscells;
s2p = load(path.s2pPath);
numRois = length(s2p.stat);
temp = readNPY(path.iscellPath);
if ~isequal(iscells,find(temp(:,1)))
    warning('iscells assignment differs between iscell.npy and sca.mat files. Using the one from the sca.mat file.')
end
numCells = length(iscells);
cal = readmatrix(info.calibrationFile)';

% import seq file
seq = readmatrix(path.seqPath);
stimTrials.idcs = find(seq(2,:)>0);
stimTrials.type = seq(2,find(seq(2,:)>0));
info.stimTypes = unique(stimTrials.type);
info.numStimTrials = length(stimTrials.idcs);
info.numTrialsPerGroup = info.numStimTrials / 2;
info.numStimTrialsPerBlock = length(find(seq(2,1:info.trialsPerBlock)>0));
if mod(info.numTrialsPerGroup*2,info.numStimTrialsPerBlock)~=0
    error('Number of stim trials is not a multiple of number of stim trials per block.')
end
info.numStimBlocks = info.numStimTrials/info.numStimTrialsPerBlock;


%% Process data: suite2p -> cmp

disp(['- Processing data.'])

% medians
med1 = zeros(numRois,1);
med2 = zeros(numRois,1);
for c=1:numRois
    med1(c,1) = s2p.stat{1,c}.med(1);
    med2(c,1) = s2p.stat{1,c}.med(2);
end

% pairwise distance compatibility matrices
pwd = squareform(pdist([med1,med2],'euclidean'));
cmp_within = ones(size(pwd));
cmp_within(1:numRois+1:end) = 0;
cmp_within(pwd>p.maxClusterExtent/info.pixelSize) = 0;
cmp_within(pwd<p.minDistanceWithin/info.pixelSize) = 0;
cmp_across = ones(size(pwd));
cmp_across(1:numRois+1:end) = 0;
cmp_across(pwd<p.minDistanceAcross/info.pixelSize) = 0;
cmp = cmp_within & cmp_across;


%% Process data: sca -> cand

% get candidate pool
cand.pool_A_all = find(sca.passed.Aonly==1);
cand.pool_X_all = find(sca.passed.Xonly==1);

% discard cells close to FOV edges from candidate pool
cand.pool_A = [];
if p.discardFOVEdges~=0
    for i=1:length(cand.pool_A_all)
        this_idx = cand.pool_A_all(i);
        if all([med1(this_idx)>ceil(p.discardFOVEdges/info.pixelSize),...
                med2(this_idx)>ceil(p.discardFOVEdges/info.pixelSize),...
                med1(this_idx)<=info.numPixels-ceil(p.discardFOVEdges/info.pixelSize),...
                med2(this_idx)<=info.numPixels-ceil(p.discardFOVEdges/info.pixelSize)])
            cand.pool_A = [cand.pool_A; this_idx];
        end
    end
end
cand.pool_X = [];
if p.discardFOVEdges~=0
    for i=1:length(cand.pool_X_all)
        this_idx = cand.pool_X_all(i);
        if all([med1(this_idx)>ceil(p.discardFOVEdges/info.pixelSize),...
                med2(this_idx)>ceil(p.discardFOVEdges/info.pixelSize),...
                med1(this_idx)<=info.numPixels-ceil(p.discardFOVEdges/info.pixelSize),...
                med2(this_idx)<=info.numPixels-ceil(p.discardFOVEdges/info.pixelSize)])
            cand.pool_X = [cand.pool_X; this_idx];
        end
    end
end

% bin candidate cells by firing field peak
cand.candidates_A = {};
for i=1:p.numClusterSlots
    cand.candidates_A{i} = [];
    for j=1:length(cand.pool_A)
        this_idx = cand.pool_A(j);
        if sca.firingField_A.peakLocation_s(this_idx)>=info.stimEdges(i) && sca.firingField_A.peakLocation_s(this_idx)<info.stimEdges(i+1)
            cand.candidates_A{i} = [cand.candidates_A{i}; this_idx];
        end
    end
end
cand.candidates_X = {};
for i=1:p.numClusterSlots
    cand.candidates_X{i} = [];
    for j=1:length(cand.pool_X)
        this_idx = cand.pool_X(j);
        if sca.firingField_X.peakLocation_s(this_idx)>=info.stimEdges(i) && sca.firingField_X.peakLocation_s(this_idx)<info.stimEdges(i+1)
            cand.candidates_X{i} = [cand.candidates_X{i}; this_idx];
        end
    end
end

% group candidate cells to compatible candidate cliques
cand.cliques_A = {};
cand.maximalCliques_A = {};
cand.maximalCliques_A_num = zeros(1,p.numClusterSlots);
for i=1:p.numClusterSlots
    if ~isempty(cand.candidates_A{i})
        cand.cliques_A{i} = maximalCliques( cmp(cand.candidates_A{i},cand.candidates_A{i}) );
        [~,temp] = max( sum(cand.cliques_A{i},1) );
        cand.maximalCliques_A{i} = cand.candidates_A{i}(find(cand.cliques_A{i}(:,temp)));
        cand.maximalCliques_A_num(i) = length(cand.maximalCliques_A{i});
    else
        cand.cliques_A{i} = [];
        cand.maximalCliques_A{i} = [];
    end
end
cand.cliques_X = {};
cand.maximalCliques_X = {};
cand.maximalCliques_X_num = zeros(1,p.numClusterSlots);
for i=1:p.numClusterSlots
    if ~isempty(cand.candidates_X{i})
        cand.cliques_X{i} = maximalCliques( cmp(cand.candidates_X{i},cand.candidates_X{i}) );
        [~,temp] = max( sum(cand.cliques_X{i},1) );
        cand.maximalCliques_X{i} = cand.candidates_X{i}(find(cand.cliques_X{i}(:,temp)));
        cand.maximalCliques_X_num(i) = length(cand.maximalCliques_X{i});
    else
        cand.cliques_X{i} = [];
        cand.maximalCliques_X{i} = [];
    end
end


%% Finalise output (clustering, grouping, selections, cluster_centres, galvo, sequences, oris, stim power)

disp(['- Finalising output.'])

% clustering
clustering = nan(nanmax([cand.maximalCliques_A_num,cand.maximalCliques_X_num]),2*p.numClusterSlots);
for i=1:p.numClusterSlots
    temp = cand.maximalCliques_A{i};
    if ~isempty(temp)
        clustering(1:length(temp),i) = temp;
    else
        clustering(1,i) = 0;
    end
end
for i=1:p.numClusterSlots
    temp = cand.maximalCliques_X{i};
    if ~isempty(temp)
        clustering(1:length(temp),i+p.numClusterSlots) = temp;
    else
        clustering(1,i+p.numClusterSlots) = 0;
    end
end

% grouping
grouping = [1:p.numClusterSlots;p.numClusterSlots+1:p.numClusterSlots*2]';

% selection_A and selection_X
selection_A = clustering(:,grouping(:,1));
selection_A = selection_A(:);
selection_A = sort(nonzeros(selection_A(~isnan(selection_A))));
selection_X = clustering(:,grouping(:,2));
selection_X = selection_X(:);
selection_X = sort(nonzeros(selection_X(~isnan(selection_X))));

% cluster_centres and galvo
cluster_centres = zeros(2,p.numClusterSlots*2);
galvo = zeros(2,p.numClusterSlots*2);
for i=1:p.numClusterSlots*2
    if clustering(1,i)==0
        cluster_centres(1,i) = info.numPixels/2;
        cluster_centres(2,i) = info.numPixels/2;
    else
        temp = clustering(~isnan(clustering(:,i)),i);
        cluster_centres(1,i) = (min(med1(temp)) + max(med1(temp))) / 2;
        cluster_centres(2,i) = (min(med2(temp)) + max(med2(temp))) / 2;
    end
    
    if clustering(1,i)==0
        galvo(1,i) = info.numPixels/2;
        galvo(2,i) = info.numPixels/2;
    else
        temp2 = zeros(info.numPixels,info.numPixels);
        temp = zeros(numel(temp2),2);
        [temp(:,1),temp(:,2)] = find(temp2==0);
        temp3 = clustering(~isnan(clustering(:,i)),i);
        temp2 = reshape(min(pdist2(temp,[med1(temp3),med2(temp3)]),[],2) > (p.galvoDisplacement/info.pixelSize),info.numPixels,info.numPixels);
        this_map = reshape(pdist2(temp,cluster_centres(:,i)'),info.numPixels,info.numPixels);
        this_map(~temp2) = NaN;
        temp = nanmin(this_map(:));
        temp2=[];
        [temp2(1,:),temp2(2,:)] = find(this_map==temp);
        galvo(:,i) = temp2(:,randsample(1:size(temp2,2),1));
    end
    
%     mat = this_map;
%     [r, c] = size(mat); 
%     figure; imagesc((1:c)+0.5, (1:r)+0.5, mat); hold on; plot(cluster_centres(2,i),cluster_centres(1,i),'ro'); plot(galvo(2,i),galvo(1,i),'r+'); hold off; colormap(parula); axis equal;
end

% sequenceClusters
sequenceClusters = zeros(p.numClusterSlots,2+info.numStimTrials/2);
sequenceClusters(:,1) = grouping(:,1);    
sequenceClusters(:,2) = grouping(:,2);
for i=1:info.numStimTrials/4
    sequenceClusters(:,2+i) = grouping(randperm(p.numClusterSlots),1);
end
for i=1:info.numStimTrials/4
    sequenceClusters(:,2+info.numStimTrials/4+i) = grouping(randperm(p.numClusterSlots),2);
end

% stimTrials.idcs

% - ispi_seq (stimTrials.type==1)
temp1 = stimTrials.idcs(stimTrials.type==1);
temp2 = sort([find(seq(1,:)==1),find(seq(1,:)==3),find(seq(1,:)==5),find(seq(1,:)==7)]);
[~,stimTrials.ispi_seq_A] = intersect(stimTrials.idcs,intersect(temp1,temp2));
temp2 = sort([find(seq(1,:)==2),find(seq(1,:)==4),find(seq(1,:)==6),find(seq(1,:)==8)]);
[~,stimTrials.ispi_seq_X] = intersect(stimTrials.idcs,intersect(temp1,temp2));

% - ispi_ctrl (stimTrials.type==2)
temp1 = stimTrials.idcs(stimTrials.type==2);
temp2 = sort([find(seq(1,:)==1),find(seq(1,:)==3),find(seq(1,:)==5),find(seq(1,:)==7)]);
[~,stimTrials.ispi_ctrl_A] = intersect(stimTrials.idcs,intersect(temp1,temp2));
temp2 = sort([find(seq(1,:)==2),find(seq(1,:)==4),find(seq(1,:)==6),find(seq(1,:)==8)]);
[~,stimTrials.ispi_ctrl_X] = intersect(stimTrials.idcs,intersect(temp1,temp2));

% - contra_seq (stimTrials.type==3)
temp1 = stimTrials.idcs(stimTrials.type==3);
temp2 = sort([find(seq(1,:)==1),find(seq(1,:)==3),find(seq(1,:)==5),find(seq(1,:)==7)]);
[~,stimTrials.contra_seq_X] = intersect(stimTrials.idcs,intersect(temp1,temp2));
temp2 = sort([find(seq(1,:)==2),find(seq(1,:)==4),find(seq(1,:)==6),find(seq(1,:)==8)]);
[~,stimTrials.contra_seq_A] = intersect(stimTrials.idcs,intersect(temp1,temp2));

% - contra_ctrl (stimTrials.type==4)
temp1 = stimTrials.idcs(stimTrials.type==4);
temp2 = sort([find(seq(1,:)==1),find(seq(1,:)==3),find(seq(1,:)==5),find(seq(1,:)==7)]);
[~,stimTrials.contra_ctrl_X] = intersect(stimTrials.idcs,intersect(temp1,temp2));
temp2 = sort([find(seq(1,:)==2),find(seq(1,:)==4),find(seq(1,:)==6),find(seq(1,:)==8)]);
[~,stimTrials.contra_ctrl_A] = intersect(stimTrials.idcs,intersect(temp1,temp2));

% sequenceOrder
sequenceOrder = zeros(1,info.numStimTrials);
sequenceOrder(stimTrials.ispi_seq_A) = 1;
sequenceOrder(stimTrials.contra_seq_A) = 1;
sequenceOrder(stimTrials.ispi_seq_X) = 2;
sequenceOrder(stimTrials.contra_seq_X) = 2;
sequenceOrder(stimTrials.ispi_ctrl_A) = 2+(1:info.numStimTrials/8);
sequenceOrder(stimTrials.contra_ctrl_A) = 2+1*info.numStimTrials/8+(1:info.numStimTrials/8);
sequenceOrder(stimTrials.ispi_ctrl_X) = 2+2*info.numStimTrials/8+(1:info.numStimTrials/8);
sequenceOrder(stimTrials.contra_ctrl_X) = 2+3*info.numStimTrials/8+(1:info.numStimTrials/8);

% sequences
sequences = zeros(p.numClusterSlots,info.numStimTrials);
for i=1:info.numStimTrials
    sequences(:,i) = sequenceClusters(:,sequenceOrder(i));
end

% oris (stim order for STA movie maker)
oris = NaN(1,info.numStimTrials*p.numClusterSlots);
for i=1:info.numStimTrials
    oris((i-1)*p.numClusterSlots+1) = sequenceOrder(i);
end

% assign stim powers
info.calibration = cal;
temp1 = 0:0.1:max(info.calibration(:,1));
temp2 = interp1(info.calibration(:,1),info.calibration(:,2),temp1);
xml.power = {}; % power per cluster [% of 1V output]
for i=1:p.numClusterSlots*2
    temp = clustering(:,i);
    temp = temp(~isnan(temp));
    if ~isempty(nonzeros(temp))
        this_num = length(temp);
        if this_num==1
            this_mW = this_num*p.powerPerCell_1cell;
        elseif this_num==2
            this_mW = this_num*p.powerPerCell_2cells;
        else
            this_mW = this_num*p.powerPerCell;
        end
        [~,temp] = min(abs(temp2-this_mW));
        temp = temp1(temp);
        xml.power{i} = num2str(temp,2);
        %plot(this_num,this_mW,'ro')
    else
        xml.power{i} = num2str(0,2);
    end
end


%% --- Plot results ---

disp(['- Plotting results.'])

default_figure();
close;


%% Plot - Target groups on FOV

F1 = figure;
hold on
scatter(med2(iscells),med1(iscells),'.','MarkerEdgeColor',[0.8,0.8,0.8])
temp2 = [p.col.AB_rainbow;p.col.XY_rainbow];
for i=1:p.numClusterSlots*2
    temp1 = clustering(:,i);
    temp1 = nonzeros(temp1(~isnan(temp1)));
    l(i) = scatter(med2(temp1),med1(temp1),'d','MarkerEdgeColor','k','MarkerFaceColor',temp2(i,:));
end
xline(ceil(p.discardFOVEdges/info.pixelSize),'LineStyle',':');
xline(info.numPixels-ceil(p.discardFOVEdges/info.pixelSize),'LineStyle',':');
xline(info.numPixels,'LineStyle','-');
yline(ceil(p.discardFOVEdges/info.pixelSize),'LineStyle',':');
yline(info.numPixels-ceil(p.discardFOVEdges/info.pixelSize),'LineStyle',':');
yline(info.numPixels,'LineStyle','-');
hold off
pbaspect([1 1 1]);
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlim([0,info.numPixels])
xticks([0:info.numPixels/4:info.numPixels])
xlabel('medial <-> lateral (pixels)')
ylim([0,info.numPixels])
yticks([0:info.numPixels/4:info.numPixels])
ylabel('posterior <-> anterior (pixels)')
title('Target groups on FOV')


%% Plot - Target clusters on FOV

nrows = 4;
ncols = 5;

temp1 = clustering(:,grouping(:,1));
temp2 = clustering(:,grouping(:,2));
temp1_2 = galvo(:,grouping(:,1));
temp2_2 = galvo(:,grouping(:,2));
temp_col = [p.col.AB_rainbow;p.col.XY_rainbow];

F2 = figure;
for i=1:p.numClusterSlots
    subplot(nrows,ncols,i)
    hold on
    scatter(med2(iscells),med1(iscells),'.','MarkerEdgeColor',[0.8,0.8,0.8])
    
    temp = clustering(:,i);
    temp = nonzeros(temp(~isnan(temp)));
    scatter(med2(temp),med1(temp),'d','MarkerEdgeColor','k','MarkerFaceColor',temp_col(i,:));
    temp = clustering(:,i+20);
    temp = nonzeros(temp(~isnan(temp)));
    scatter(med2(temp),med1(temp),'s','MarkerEdgeColor','k','MarkerFaceColor',temp_col(i+20,:));
    scatter(temp1_2(2,i),temp1_2(1,i),'x','MarkerEdgeColor',temp_col(i,:),'LineWidth',1)
    scatter(temp2_2(2,i),temp2_2(1,i),'x','MarkerEdgeColor',temp_col(i+20,:),'LineWidth',1)
    xline(ceil(p.discardFOVEdges/info.pixelSize),'LineStyle',':');
    xline(info.numPixels-ceil(p.discardFOVEdges/info.pixelSize),'LineStyle',':');
    xline(info.numPixels,'LineStyle','-');
    yline(ceil(p.discardFOVEdges/info.pixelSize),'LineStyle',':');
    yline(info.numPixels-ceil(p.discardFOVEdges/info.pixelSize),'LineStyle',':');
    yline(info.numPixels,'LineStyle','-');
    hold off
    pbaspect([1 1 1]);
    set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
    xlim([0,info.numPixels])
    xticks([])
    xticklabels([])
    ylim([0,info.numPixels])
    yticks([])
    yticklabels([])
    xlabel(['Cluster ',num2str(i),', t = ',num2str(info.stimEdges(i),2),' s'])
end
han=axes(F2,'visible','off'); han.Title.Visible='on'; han.XLabel.Visible='on'; han.YLabel.Visible='on';
xlabel(han,'medial <-> lateral');
ylabel(han,'posterior <-> anterior');
suptitle('Target clusters on FOV');


%% Plot - Target groups on FOV

F3 = figure;
subplot(1,2,1)
hold on
temp_col = [p.col.AB_rainbow;p.col.XY_rainbow];
for i=1:p.numClusterSlots
    temp = clustering(:,i);
    temp = nonzeros(temp(~isnan(temp)));
    for j=1:length(temp)
        this_data = permute(squeeze(sca.traces.traces_A(temp(j),:,:)),[2,1]);
        this_lower = nanmin(nanmean(this_data,1));
        this_upper = nanmax(nanmean(this_data,1));
        this_normData = (this_data-this_lower) ./ (this_upper-this_lower);      
        temp2=shadedErrorBar(1:size(this_normData,2),nanmean(this_normData,1)+p.numClusterSlots+1-i,nansem(this_normData,1),'lineProps',temp_col(i,:)); temp2.mainLine.LineWidth = 1;   
        hold on
    end
end
plt = struct(); traces_task(sca,p,info,plt); 
for i=1:p.numClusterSlots
    temp = clustering(:,i);
    temp = nonzeros(temp(~isnan(temp)));
    temp1 = interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.stimEdges(i));
    temp2 = interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.stimEdges(i)+(str2num(xml.durationMS)+str2num(xml.postIteIdleMS))*str2num(xml.iterations)/1000) - temp1;
    rectangle('Position',[temp1,p.numClusterSlots+1-i+0.3,temp2,0.7],'FaceColor',p.col.photostim,'EdgeColor','none','LineWidth',0.5) 
end
hold off
ylim([0,p.numClusterSlots+3])
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
title('seqA')
xlabel('time (s)')
ylabel('normalised activity (stacked by stim cluster)')

subplot(1,2,2)
hold on
temp_col = [p.col.AB_rainbow;p.col.XY_rainbow];
for i=1:p.numClusterSlots
    temp = clustering(:,i+p.numClusterSlots);
    temp = nonzeros(temp(~isnan(temp)));
    for j=1:length(temp)
        this_data = permute(squeeze(sca.traces.traces_X(temp(j),:,:)),[2,1]);
        this_lower = nanmin(nanmean(this_data,1));
        this_upper = nanmax(nanmean(this_data,1));
        this_normData = (this_data-this_lower) ./ (this_upper-this_lower);      
        temp2=shadedErrorBar(1:size(this_normData,2),nanmean(this_normData,1)+p.numClusterSlots+1-i,nansem(this_normData,1),'lineProps',temp_col(i+p.numClusterSlots,:)); temp2.mainLine.LineWidth = 1;   
        hold on
    end
end
plt = struct(); traces_task(sca,p,info,plt); 
for i=1:p.numClusterSlots
    temp = clustering(:,i);
    temp = nonzeros(temp(~isnan(temp)));
    temp1 = interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.stimEdges(i));
    temp2 = interp1(sca.prop.t_binned,1:length(sca.prop.t_binned),info.stimEdges(i)+(str2num(xml.durationMS)+str2num(xml.postIteIdleMS))*str2num(xml.iterations)/1000) - temp1;
    rectangle('Position',[temp1,p.numClusterSlots+1-i+0.3,temp2,0.7],'FaceColor',p.col.photostim,'EdgeColor','none','LineWidth',0.5) 
end
hold off
ylim([0,p.numClusterSlots+3])
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
title('seqX')
xlabel('time (s)')
ylabel('normalised activity (stacked by stim cluster)')

suptitle('Endogenous sequence targets');

    
%% Save results

disp(['- Saving results.'])

f.med1 = med1;
f.med2 = med2;
trg.info = info;
trg.path = path;
trg.p = p;
trg.xml = xml;
trg.s2p = s2p;
trg.iscells = iscells;
trg.seq = seq;
trg.stimTrials = stimTrials;
trg.f = f;
trg.pwd = pwd;
trg.cmp_within = cmp_within;
trg.cmp_across = cmp_across;
trg.cmp = cmp;
trg.clustering = clustering;
trg.grouping = grouping;
trg.selection_A = selection_A;
trg.selection_X = selection_X;
trg.cluster_centres = cluster_centres;
trg.galvo = galvo;
trg.sequenceClusters = sequenceClusters;
trg.sequenceOrder = sequenceOrder;
trg.sequences = sequences;
trg.oris = oris;
trg.sca = sca;
trg.cand = cand;

if ~exist(path.savePath,'dir')
    mkdir(path.savePath);
end

temp = datestr(now,'YYYYmmdd_HHMMSS');
writeTargetXML_ExpertStim(trg,[path.savePath,['Targeting_',info.animal,'_',temp,'.xml']]);
save([path.savePath,['Targeting_',info.animal,'_',temp,'.mat']],'trg','-v7.3');
save([path.savePath,['oris_',info.animal,'_',temp,'.mat']],'oris');
savefig(F1,[path.savePath,['Targeting_',info.animal,'_',temp,'_F1.fig']]);
savefig(F2,[path.savePath,['Targeting_',info.animal,'_',temp,'_F2.fig']]);
savefig(F3,[path.savePath,['Targeting_',info.animal,'_',temp,'_F3.fig']]);

temp = toc;
disp(['- Done. Elapsed time: ',num2str(floor(temp/60)),' minutes and ',num2str(floor(rem(temp,60))),' seconds.']);

