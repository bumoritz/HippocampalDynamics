%% Continued targeting

path.sourceFile = 'G:\Data\SniffinHippo\2021\2021-12\2021-12-13\Python\Targeting\TargetSelection\Targeting_Python_20211213_111507.mat'; % Here comes the .mat file from the target selection folder of the previous day

info.animal             = 'Python'; %'Mufasa';
info.date               = '20211214'; %'20200520';
info.rootPath           = 'N:\Data\SniffinHippo\'; % 'G:\Data\SniffinHippo\';
% will save stuff here: e.g. 'C:\Users\Moritz\Desktop\OnlineAnalysis\2020\2020-10\2020-10-01\MB040\Targeting\TargetSelection\'

info.trialsPerBlock     = 20;

disp(['Continued targeting.'])


%% Find paths and load data

disp(['- Loading data.'])

% import source file
trg_source = load(path.sourceFile);
trg_source = trg_source.trg;

% extract relevant source file components
grouping = trg_source.grouping;
p = trg_source.p;

% find paths
info.year = info.date(1:4);
info.month = info.date(5:6);
info.day = info.date(7:8);
path.homeFolder = [info.rootPath,info.year,'\',info.year,'-',info.month,'\',info.year,'-',info.month,'-',info.day,'\',info.animal,'\','Targeting','\'];
temp = dir([path.homeFolder,info.animal,'_',info.date,'_*_SEQ.txt']);
path.seqPath = [path.homeFolder,temp.name];
path.savePath = [path.homeFolder,'TargetSelection\'];

% import seq file
seq = readmatrix(path.seqPath);
seq_stimTrials = seq(1,seq(2,:)==1);
p.numStimTrials = sum(seq(2,:),2);
p.numTrialsPerGroup = p.numStimTrials / 2;
p.numStimTrialsPerBlock = sum(seq(2,1:info.trialsPerBlock),2);
if mod(p.numTrialsPerGroup*2,p.numRepsPerCtrl)~=0
    error('Number of stim trials is not a multiple of number of repetitions for each ctrl.')
end
if mod(p.numTrialsPerGroup*2,p.numStimTrialsPerBlock)~=0
    error('Number of stim trials is not a multiple of number of stim trials per block.')
end
p.numDifferentCtrlSeqs  = p.numTrialsPerGroup / p.numRepsPerCtrl;
p.numStimTrialsPerSuperBlock = p.numStimTrials / p.numRepsPerCtrl;
p.numBlocks = p.numStimTrials/p.numStimTrialsPerBlock;


%% Core

disp(['- Processing data.'])

%%% new from 20210607 %%%

if strcmp(p.stimType,'seq')
    
%     % define sequences clusters
%     sequenceClusters = zeros(p.numClustersPerGroup,2);
%     sequenceClusters(:,1) = grouping(:,1);    
%     sequenceClusters(:,2) = grouping(randperm(p.numClustersPerGroup),2);
%     
    sequenceClusters = trg_source.sequenceClusters; % MB20210813
    
    % assign order of sequences
    sequenceOrder = [];
    for j=1:p.numBlocks/p.numRepsPerCtrl
        this_seq = seq_stimTrials((j-1)*p.numStimTrialsPerBlock+1:j*p.numStimTrialsPerBlock);
        temp = zeros(1,p.numStimTrialsPerBlock);
        temp(this_seq==1|this_seq==3) = 1;
        temp(this_seq==2|this_seq==4) = 2;
        sequenceOrder = [sequenceOrder, temp];
    end


elseif strcmp(p.stimType,'ctrl')
    
    % define sequences clusters
    sequenceClusters = zeros(p.numClustersPerGroup,p.numDifferentCtrlSeqs*2);
    for i=1:p.numDifferentCtrlSeqs
        sequenceClusters(:,i) = grouping(randperm(p.numClustersPerGroup),1);
        sequenceClusters(:,p.numDifferentCtrlSeqs+i) = grouping(randperm(p.numClustersPerGroup),2);
    end
    
    % assign order of sequences
    sequenceOrder = [];
    this_pool_A = 1:p.numDifferentCtrlSeqs;
    this_pool_X = p.numDifferentCtrlSeqs+1:2*p.numDifferentCtrlSeqs;
    for j=1:p.numBlocks/p.numRepsPerCtrl
        this_seq = seq_stimTrials((j-1)*p.numStimTrialsPerBlock+1:j*p.numStimTrialsPerBlock);
        temp = zeros(1,p.numStimTrialsPerBlock); % MB20210701 to troubleshoot below error

        temp2 = datasample(this_pool_A,8,'Replace',false);
        temp(this_seq==1|this_seq==3) = temp2; % sometimes: Assignment between unlike types is not allowed.
        this_pool_A = setdiff(this_pool_A,temp2);
        
        temp2 = datasample(this_pool_X,8,'Replace',false);
        temp(this_seq==2|this_seq==4) = temp2; 
        this_pool_X = setdiff(this_pool_X,temp2);
        
        sequenceOrder = [sequenceOrder, temp];
    end

end

% prepare stim order for STA movie maker
oris = NaN(1,p.numStimTrials*p.numClustersPerGroup);
for i=1:p.numStimTrials
    oris((i-1)*p.numClustersPerGroup+1) = sequenceOrder(i);
end


%% Save results

disp(['- Saving results.'])

trg = trg_source;
trg.info_source = trg_source.info;
trg.info = info;
trg.p_source = trg_source.p;
trg.p = p;
trg.path_source = trg_source.path;
trg.path = path;
trg.seq_source = trg_source.seq;
trg.seq = seq;
trg.sequenceClusters_source = trg_source.sequenceClusters;
trg.sequenceClusters = sequenceClusters;
trg.sequenceOrder_source = trg_source.sequenceOrder;
trg.sequenceOrder = sequenceOrder;
trg.oris_source = trg_source.oris;
trg.oris = oris;

if ~exist(path.savePath,'dir')
    mkdir(path.savePath);
end

temp = datestr(now,'YYYYmmdd_HHMMSS');
writeTargetXML(trg,[path.savePath,['Targeting_',info.animal,'_',temp,'.xml']]);
save([path.savePath,['Targeting_',info.animal,'_',temp,'.mat']],'trg','-v7.3');
save([path.savePath,['oris_',info.animal,'_',temp,'.mat']],'oris');

