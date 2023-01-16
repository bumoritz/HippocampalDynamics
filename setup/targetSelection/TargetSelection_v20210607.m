%% Select data and parameters
%%% ___MOST_RECENT_VERSION___ %%%

% suite2p QuickCuration:
% 1) include all ROIs
% 2) compact: remove top n -> removes scanning artefacts (threshold: ~1.1)
% 3) skew: add top n back in, then remove bottom n -> selects for cells with sparse transients (threshold: ~1)
% 4) draw selection to remove big bad areas

clc;
clear;
tic;

info.animal             = 'Python'; %'Mufasa';
info.date               = '20211213'; %'20200520';
info.recording          = 'Python-20211213-spont1'; %'Mufasa-20200520-5min';
info.rootPath           = 'G:\Data\SniffinHippo\'; % 'G:\Data\SniffinHippo\';
% will save stuff here: e.g. 'BASEPATH\2020\2020-10\2020-10-01\MB040\Targeting\TargetSelection\'


%%% --- critical parameters --- %%%

p.stimType              = 'seq'; % seq or ctrl
p.groupNames            = ["seqA","seqX"]; % ["seqA","seqX"] or ["ctrlA","ctrlX"]
xml.power               = '5.2'; % power per cluster [% of 1V output]

%%% --------------------------- %%%


% group assignment
p.numCellsPerGroup      = 120;
p.discardFOVEdges       = 50; % [um]
p.minDistanceAcross     = 8; % [um] (above 10 takes long, even above 8)

% cluster assignment
p.numClustersPerGroup   = 20;
p.numBufferClusters     = 10;
p.maxClusterExtent      = 450; % [um](below 400-450 can take very long)
p.minDistanceWithin     = 50; % [um]
p.clustAssignmentIters  = 1000;
p.galvoDisplacement     = 50; % [um] %MB20210521: changed from 100 to 50

p.numRepsPerCtrl = 1; % added as hail mary hack attempt NR

% plot labels
p.groupCols             = [136,19,80;4,88,156]/255;

% experiment conditions
info.trialsPerBlock     = 20;
info.FOVsize            = 1000.78; % [um]
info.numPixels          = 512;
p.smoothingSD           = 2; % [frames]
p.baselineSdThreshold   = 3;
p.neuropilSubtraction   = true;

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

if mod(p.numCellsPerGroup,p.numClustersPerGroup)~=0
    error('Number of cells per group is not a multiple of number of cells per cluster.')
end
p.numCellsPerCluster    = p.numCellsPerGroup / p.numClustersPerGroup;
p.numClusters           = p.numClustersPerGroup*2 + p.numBufferClusters;
info.pixelSize          = info.FOVsize / info.numPixels;
xml.roiWidthPx          = num2str(xml.roiWidthUM / info.pixelSize);
xml.roiHeightPx         = num2str(xml.roiHeightUM / info.pixelSize);

disp(['Target selection.'])


%% Load and process data

disp(['- Loading data.'])

info.year = info.date(1:4);
info.month = info.date(5:6);
info.day = info.date(7:8);
path.homeFolder = [info.rootPath,info.year,'\',info.year,'-',info.month,'\',info.year,'-',info.month,'-',info.day,'\',info.animal,'\','Targeting','\'];
path.s2pPath = [path.homeFolder,info.recording,'\suite2p\plane0\Fall.mat'];
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

% import s2p file
s2p = load(path.s2pPath);
numROIs = length(s2p.stat);
iscells = find(s2p.iscell(:,1));
iscells_all = iscells;
if p.discardFOVEdges~=0
    temp = iscells;
    for i=1:length(iscells)
        if ((s2p.stat{1,temp(i)}.med(1)<=ceil(p.discardFOVEdges/info.pixelSize))|...
                (s2p.stat{1,temp(i)}.med(2)<=ceil(p.discardFOVEdges/info.pixelSize))|...
                (s2p.stat{1,temp(i)}.med(1)>info.numPixels-ceil(p.discardFOVEdges/info.pixelSize))|...
                (s2p.stat{1,temp(i)}.med(2)>info.numPixels-ceil(p.discardFOVEdges/info.pixelSize)))
            iscells = iscells(iscells~=temp(i));
        end
    end
end
numCells = length(iscells);

disp(['- Processing data.'])

if p.neuropilSubtraction
    F = s2p.F - 0.7*s2p.Fneu;
else
    F = s2p.F;
end
% was like this until 20210426... temp = F; temp(~(F < p.baselineSdThreshold*std(F,[],2)))=NaN; F0 = nanmean(temp,2);
F0 = nanmedian(F,2);
dFF = (F-F0)./F0;
if p.smoothingSD~=0
    dFF = smoothdata(dFF,2,'gaussian',5*p.smoothingSD);
end

f.med1 = zeros(numCells,1);
f.med2 = zeros(numCells,1);
f.npix = zeros(numCells,1);
f.compact = zeros(numCells,1);
f.footprint = zeros(numCells,1);
f.aspect_ratio = zeros(numCells,1);
f.bsl = zeros(numCells,1);
f.mean = zeros(numCells,1);
f.std = zeros(numCells,1);
f.skew = zeros(numCells,1);
for c=1:numCells
    
    % localisation features
    f.med1(c,1) = s2p.stat{1,iscells(c)}.med(1);
    f.med2(c,1) = s2p.stat{1,iscells(c)}.med(2);
    
    % morphology features
    f.npix(c,1) = s2p.stat{1,iscells(c)}.npix;
    f.compact(c,1) = s2p.stat{1,iscells(c)}.compact;
    f.footprint(c,1) = s2p.stat{1,iscells(c)}.footprint;
    f.aspect_ratio(c,1) = s2p.stat{1,iscells(c)}.aspect_ratio;
    
    % activity features
    f.bsl(c,1) = F0(iscells(c),:);
    f.mean(c,1) = nanmean(dFF(iscells(c),:)); %mean(s2p.F(iscells(c),:)-0.7*s2p.Fneu(iscells(c),:));
    f.std(c,1) = nanstd(dFF(iscells(c),:)); %s2p.stat{1,iscells(c)}.std; % is calculated like this std(s2p.F(1,:)-0.7*s2p.Fneu(1,:))
    f.skew(c,1) = skewness(dFF(iscells(c),:)); %s2p.stat{1,iscells(c)}.skew; % is calculated like this skewness(s2p.F(1,:)-0.7*s2p.Fneu(1,:))
end
f_labels = fields(f);

f_vec = [];
f_z_vec = [];
for i=1:length(f_labels)
    f_z.(f_labels{i}) = nanzscore(f.(f_labels{i}));
    f_vec = [f_vec, f.(f_labels{i})];
    f_z_vec = [f_z_vec, nanzscore(f.(f_labels{i}))];
end

% coordinates of all cells for plotting
for c=1:length(iscells_all)
    med1_all(c,1) = s2p.stat{1,iscells_all(c)}.med(1);
    med2_all(c,1) = s2p.stat{1,iscells_all(c)}.med(2);
end

% pairwise metrics
pwd = squareform(pdist([f.med1,f.med2],'euclidean'));
pwc = 1-squareform(pdist(dFF(iscells,:),'correlation')); % it's Pearson

pwd_z = squareform(nanzscore(squareform(pwd)));
temp = pwc;
temp(1:numCells+1:end) = 0;
pwc_z = squareform(nanzscore(squareform(temp)));

numFeatures = length(f_labels)+2;

% pairwise compatibility matrices
cmp_within = ones(size(pwd));
cmp_within(1:numCells+1:end) = 0;
cmp_within(pwd>p.maxClusterExtent/info.pixelSize) = 0;
cmp_within(pwd<p.minDistanceWithin/info.pixelSize) = 0;

cmp_across = ones(size(pwd));
cmp_across(1:numCells+1:end) = 0;
cmp_across(pwd<p.minDistanceAcross/info.pixelSize) = 0;

cmp = cmp_within & cmp_across;


%% Assign cells into tentative clusters

disp(['- Identifying allowed clusterings.'])

n=0;
clusterings = zeros(p.numCellsPerCluster,p.numClusters,p.clustAssignmentIters);
while n<p.clustAssignmentIters
    try

        % initialise
        this_pool = 1:numCells;
        this_clustering = zeros(p.numCellsPerCluster,p.numClusters);

        % randomly draw clusters from cmp matrices
        for j=1:p.numClusters
            for k=1:p.numCellsPerCluster
                
                % draw seeds
                if j==1 && k==1 % MB20210513: corrected from j==k==1
                    this_clustering(k,j) = randsample(this_pool,1);
                elseif k==1
                    [~,~,temp] = find(this_clustering(:));
                    these_compatible_across = find( all(cmp_across(temp,:),1) );
                    this_clustering(k,j) = randsample( intersect(these_compatible_across,this_pool) ,1);
                    
                % draw compatible cells
                else
                    [~,~,temp] = find(this_clustering(:));
                    these_compatible_across = find( all(cmp_across(temp,:),1) );     
                    [~,~,temp] = find(this_clustering(:,j));
                    these_compatible_within = find( all(cmp_within(temp,:),1) );                  
                    this_clustering(k,j) = randsample( intersect(intersect(these_compatible_across,these_compatible_within),this_pool) ,1);
                    if ~all(squareform(cmp(nonzeros(this_clustering(:,j)),nonzeros(this_clustering(:,j)))))
                        error(['Just drew a cluster that doesnt fulfill criteria. j, k: ',num2str(j),',',num2str(k)])
                    end
                end
                
                % remove selected cells from pool
                this_pool = setdiff(this_pool,this_clustering(:));
            end
            
            % sort
            this_clustering(:,j) = sort(this_clustering(:,j));
        end
        
        n=n+1;
        if mod(n,100)==0
            disp(['--- ',num2str(n),' allowed clusterings identified.'])
        end
        clusterings(:,:,n) = this_clustering;     
    catch
%        warning('Jumped into catch. Probably too little cells.')
    end
end


%% Assign tentative clusters into two groups

disp(['- Separating clusters into two groups.'])

groupings = zeros(p.numClustersPerGroup,2,p.clustAssignmentIters);
for i=1:p.clustAssignmentIters
    this_clustering = clusterings(:,:,i);
    
    % collect average features
    these_avg_features = zeros(p.numClusters,numFeatures);
    for j=1:p.numClusters
        these_avg_features(j,1:end-2) = nanmean(f_z_vec(this_clustering(:,j),:),1);
        these_avg_features(j,end-1) = nanmean(squareform(pwd_z(this_clustering(:,j),this_clustering(:,j))));
        these_avg_features(j,end) = nanmean(squareform(pwc_z(this_clustering(:,j),this_clustering(:,j))));
    end
    
    % identify cluster pairs
    these_pwfd = squareform(pdist(these_avg_features,'euclidean'));
    for j=1:p.numClustersPerGroup
        [this_pair,~]=find(these_pwfd==nanmin(nanmin(squareform(these_pwfd))));
        
        % assign pair elements to different groups
        if j==1
            groupings(j,:,i) = this_pair';
        else
            these_avg_features_G1 = these_avg_features(groupings(find(groupings(:,1,i)),1,i),:);
            these_avg_features_G2 = these_avg_features(groupings(find(groupings(:,2,i)),2,i),:);
            this_dist_E1_G1 = pdist2(these_avg_features(this_pair(1),:),nanmean(these_avg_features_G1,1),'euclidean');
            this_dist_E1_G2 = pdist2(these_avg_features(this_pair(1),:),nanmean(these_avg_features_G2,1),'euclidean');
            this_dist_E2_G1 = pdist2(these_avg_features(this_pair(2),:),nanmean(these_avg_features_G1,1),'euclidean');
            this_dist_E2_G2 = pdist2(these_avg_features(this_pair(2),:),nanmean(these_avg_features_G2,1),'euclidean');
            if nanmean([this_dist_E1_G1,this_dist_E2_G2]) > nanmean([this_dist_E1_G2,this_dist_E2_G1])
                groupings(j,1,i) = this_pair(1);
                groupings(j,2,i) = this_pair(2);
            else
                groupings(j,1,i) = this_pair(2);
                groupings(j,2,i) = this_pair(1);
            end
        end
        
        % remove this pair from pool
        these_pwfd(this_pair,:) = NaN;
        these_pwfd(:,this_pair) = NaN;
        these_pwfd(1:p.numClusters+1:end) = 0;
    end
    groupings(:,:,i) = sortrows(groupings(:,:,i));
end


%% Generate feature list and select best assignment into groups

disp(['- Selecting best matched clustering.'])

% generate cell-wise feature lists
features.names = ["anteroposterior position (\mum)","mediolateral position (\mum)",...
    "number of pixels","compactness","footprint","aspect ratio",...
    "baseline fluorescence intensity","mean activity (\DeltaF/F)","s.d. of activity (\DeltaF/F)","skewness of activity",...
    "distance to cluster centre (\mum)","mean pw dist - FOV (\mum)","mean pw dist - group (\mum)","mean pw dist - cluster (\mum)",...
    "mean pw corr - FOV","mean pw corr - group","mean pw corr - cluster",...
    "cluster extent (\mum)"];

max_minpval = 0;
for i=1:p.clustAssignmentIters
    
    % collect cells
    temp1 = clusterings(:,groupings(:,1,i),i);
    this_selection1 = sort(temp1(:));
    this_selection1_unsorted = temp1(:);
    temp2 = clusterings(:,groupings(:,2,i),i);
    this_selection2 = sort(temp2(:));
    this_selection2_unsorted = temp2(:);
    this_selection12 = sort([this_selection1,this_selection2]);
    
    % calculate distances from cluster centres
    these_cluster_centres1 = zeros(2,p.numClustersPerGroup);
    for j=1:p.numClustersPerGroup
        these_cluster_centres1(1,j) = (min(f.med1(temp1(:,j))) + max(f.med1(temp1(:,j)))) / 2;
        these_cluster_centres1(2,j) = (min(f.med2(temp1(:,j))) + max(f.med2(temp1(:,j)))) / 2;
    end
    these_cluster_centres2 = zeros(2,p.numClustersPerGroup);
    for j=1:p.numClustersPerGroup
        these_cluster_centres2(1,j) = (min(f.med1(temp2(:,j))) + max(f.med1(temp2(:,j)))) / 2;
        these_cluster_centres2(2,j) = (min(f.med2(temp2(:,j))) + max(f.med2(temp2(:,j)))) / 2;
    end
    
    % generate cell-wise feature lists - 1
    these_features_cellw1 = zeros(p.numCellsPerGroup,17);
    for j=1:length(f_labels)
        these_features_cellw1(:,j) = f.(f_labels{j})(this_selection1_unsorted);
    end
    for j=1:p.numClustersPerGroup
        these_features_cellw1((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,11) = pdist2([these_features_cellw1((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,1),these_features_cellw1((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,2)],these_cluster_centres1(:,j)','euclidean');
    end 
    temp = pwd(this_selection1_unsorted,:);
    temp(find(temp==0))=NaN;
    these_features_cellw1(:,12) = nanmean(temp,2); % mean pwd to any cell in FOV
    temp = pwd(this_selection1_unsorted,this_selection1_unsorted);
    temp(1:length(this_selection1_unsorted)+1:end) = diag(NaN(length(this_selection1_unsorted)));
    these_features_cellw1(:,13) = nanmean(temp,2); % mean pwd to any cell in same stim group
    for j=1:p.numClustersPerGroup
        these_features_cellw1((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,14) = nanmean(temp((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,(j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster),2); % mean pwd to any cell in same cluster
    end
    temp = pwc(this_selection1_unsorted,:);
    temp(find(temp==0))=NaN;
    these_features_cellw1(:,15) = nanmean(temp,2); % mean pwc with any cell in FOV
    temp = pwc(this_selection1_unsorted,this_selection1_unsorted);
    temp(1:length(this_selection1_unsorted)+1:end) = diag(NaN(length(this_selection1_unsorted)));
    these_features_cellw1(:,16) = nanmean(temp,2); % mean pwc with any cell in same stim group
    for j=1:p.numClustersPerGroup
        these_features_cellw1((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,17) = nanmean(temp((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,(j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster),2); % mean pwc with any cell in same cluster
    end
    
    % generate cell-wise feature lists - 2
    these_features_cellw2 = zeros(p.numCellsPerGroup,17);
    for j=1:length(f_labels)
        these_features_cellw2(:,j) = f.(f_labels{j})(this_selection2_unsorted);
    end
    for j=1:p.numClustersPerGroup
        these_features_cellw2((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,11) = pdist2([these_features_cellw2((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,1),these_features_cellw2((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,2)],these_cluster_centres2(:,j)','euclidean');
    end 
    temp = pwd(this_selection2_unsorted,:);
    temp(find(temp==0))=NaN;
    these_features_cellw2(:,12) = nanmean(temp,2); % mean pwd to any cell in FOV
    temp = pwd(this_selection2_unsorted,this_selection2_unsorted);
    temp(1:length(this_selection2_unsorted)+1:end) = diag(NaN(length(this_selection2_unsorted)));
    these_features_cellw2(:,13) = nanmean(temp,2); % mean pwd to any cell in same stim group
    for j=1:p.numClustersPerGroup
        these_features_cellw2((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,14) = nanmean(temp((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,(j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster),2); % mean pwd to any cell in same cluster
    end
    temp = pwc(this_selection2_unsorted,:);
    temp(find(temp==0))=NaN;
    these_features_cellw2(:,15) = nanmean(temp,2); % mean pwc with any cell in FOV
    temp = pwc(this_selection2_unsorted,this_selection2_unsorted);
    temp(1:length(this_selection2_unsorted)+1:end) = diag(NaN(length(this_selection2_unsorted)));
    these_features_cellw2(:,16) = nanmean(temp,2); % mean pwc with any cell in same stim group
    for j=1:p.numClustersPerGroup
        these_features_cellw2((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,17) = nanmean(temp((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,(j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster),2); % mean pwc with any cell in same cluster
    end   
    
    % generate cluster-wise feature lists - 1
    these_features_clusterw1 = zeros(p.numClustersPerGroup,size(these_features_cellw1,2)+1);
    for j=1:p.numClustersPerGroup
        these_features_clusterw1(j,1:end-1) = nanmean(these_features_cellw1((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,:),1);
    end
    temp = pwd(this_selection1_unsorted,this_selection1_unsorted);
    temp(1:length(this_selection1_unsorted)+1:end) = diag(NaN(length(this_selection1_unsorted)));
    for j=1:p.numClustersPerGroup
        these_features_clusterw1(j,end) = nanmax(nanmax(temp((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,(j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster),[],2),[],1); % max pwd to any cell in same cluster
    end
    
    % generate cluster-wise feature lists - 2
    these_features_clusterw2 = zeros(p.numClustersPerGroup,size(these_features_cellw2,2)+1);
    for j=1:p.numClustersPerGroup
        these_features_clusterw2(j,1:end-1) = nanmean(these_features_cellw2((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,:),1);
    end
    temp = pwd(this_selection2_unsorted,this_selection2_unsorted);
    temp(1:length(this_selection2_unsorted)+1:end) = diag(NaN(length(this_selection2_unsorted)));
    for j=1:p.numClustersPerGroup
        these_features_clusterw2(j,end) = nanmax(nanmax(temp((j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster,(j-1)*p.numCellsPerCluster+1:j*p.numCellsPerCluster),[],2),[],1); % max pwd to any cell in same cluster
    end
       
    % calculate p values - group comparisons and cluster comparisons
    pval=NaN(18,4);
    for j=1:18
        if j<18
            pval(j,1) = ranksum(these_features_cellw1(:,j),these_features_cellw2(:,j));
            [~,pval(j,2)] = kstest2(these_features_cellw1(:,j),these_features_cellw2(:,j));
        end
        pval(j,3) = signrank(these_features_clusterw1(:,j),these_features_clusterw2(:,j));
        [~,pval(j,4)] = kstest2(these_features_clusterw1(:,j),these_features_clusterw2(:,j));
    end
    
    % update current best option
    if nanmin(pval(:)) > max_minpval
        disp(['--- current max. min. p-value: ',num2str(nanmin(pval(:)),2)])
        max_minpval = nanmin(pval(:));
        winner = i;
        clustering = clusterings(:,:,i);
        grouping = groupings(:,:,i);
        selection1 = this_selection1;
        selection2 = this_selection2;
        
        features.selection1_unsorted = this_selection1_unsorted;
        features.selection2_unsorted = this_selection2_unsorted;
        features.cellw1 = these_features_cellw1;
        features.cellw2 = these_features_cellw2;
        features.clusterw1 = these_features_clusterw1;
        features.clusterw2 = these_features_clusterw2;
        features.pval = pval;
    end
end
if max_minpval < 0.05
    warning('No assignment without significant differences between the groups was found.')
end


%% Finalise output (clustering, grouping, cluster_centres, galvo, sequences)

disp(['- Finalising output.'])

% crop additional clusters
temp = sort(grouping(:));
clustering = clustering(:,temp);
for i=1:p.numClustersPerGroup
    grouping(i,1) = find(temp==grouping(i,1));
    grouping(i,2) = find(temp==grouping(i,2));
end

% randomly assign final order of clusters for seq stim
temp = randperm(p.numClustersPerGroup);
grouping = grouping(temp,:);

% assign galvo positions
cluster_centres = zeros(2,p.numClustersPerGroup*2);
galvo = zeros(2,p.numClustersPerGroup*2);
for i=1:p.numClustersPerGroup*2
    cluster_centres(1,i) = (min(f.med1(clustering(:,i))) + max(f.med1(clustering(:,i)))) / 2;
    cluster_centres(2,i) = (min(f.med2(clustering(:,i))) + max(f.med2(clustering(:,i)))) / 2;
    
    temp2 = zeros(info.numPixels,info.numPixels);
    temp = zeros(numel(temp2),2);
    [temp(:,1),temp(:,2)] = find(temp2==0);
    
    temp2 = reshape(min(pdist2(temp,[f.med1(clustering(:,i)),f.med2(clustering(:,i))]),[],2) > (p.galvoDisplacement/info.pixelSize),info.numPixels,info.numPixels);
    this_map = reshape(pdist2(temp,cluster_centres(:,i)'),info.numPixels,info.numPixels);
    this_map(~temp2) = NaN;

    temp = nanmin(this_map(:));
    temp2=[];
    [temp2(1,:),temp2(2,:)] = find(this_map==temp);
    galvo(:,i) = temp2(:,randsample(1:size(temp2,2),1));
    
%     mat = this_map;
%     [r, c] = size(mat); 
%     figure; imagesc((1:c)+0.5, (1:r)+0.5, mat); hold on; plot(cluster_centres(2,i),cluster_centres(1,i),'ro'); plot(galvo(2,i),galvo(1,i),'r+'); hold off; colormap(parula); axis equal;
end


%%% new from 20210607 %%%

if strcmp(p.stimType,'seq')
    
    % define sequences clusters
    sequenceClusters = zeros(p.numClustersPerGroup,2);
    sequenceClusters(:,1) = grouping(:,1);    
    sequenceClusters(:,2) = grouping(randperm(p.numClustersPerGroup),2);
    
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
        temp(this_seq==1|this_seq==3) = temp2;
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


%% --- Plot results ---

disp(['- Plotting results.'])

default_figure();
close;


%% Plot - Target sequence

% this_selection1 = selection1;
% this_selection2 = selection2;
% temp1 = clustering(:,grouping(:,1));
% 
% F1 = figure;
% cmap = spring(p.numClustersPerGroup);
% hold on
% scatter(med2_all,med1_all,'.','MarkerEdgeColor',[0.8,0.8,0.8])
% l(1) = scatter(f.med2(this_selection1),f.med1(this_selection1),'d','MarkerEdgeColor','k');
% for n=1:p.numClustersPerGroup
%     scatter(f.med2(temp1(:,n)),f.med1(temp1(:,n)),'d','MarkerEdgeColor','k','MarkerFaceColor',cmap(n,:))
% end
% l(2) = scatter(f.med2(this_selection2),f.med1(this_selection2),'.','MarkerEdgeColor',p.groupCols(2,:));
% xline(ceil(p.discardFOVEdges/info.pixelSize),'LineStyle',':');
% xline(info.numPixels-ceil(p.discardFOVEdges/info.pixelSize),'LineStyle',':');
% xline(info.numPixels,'LineStyle','-');
% yline(ceil(p.discardFOVEdges/info.pixelSize),'LineStyle',':');
% yline(info.numPixels-ceil(p.discardFOVEdges/info.pixelSize),'LineStyle',':');
% yline(info.numPixels,'LineStyle','-');
% hold off
% pbaspect([1 1 1]);
% set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
% xlim([0,info.numPixels])
% xticks([0:info.numPixels/4:info.numPixels])
% xlabel('medial <-> lateral (pixels)')
% ylim([0,info.numPixels])
% yticks([0:info.numPixels/4:info.numPixels])
% ylabel('posterior <-> anterior (pixels)')
% legend(l,p.groupNames)
% colormap('spring'); cbr = colorbar; cbr.Ticks = 0:1/5:1; cbr.TickLabels = strsplit(num2str(0:1:5))';  cbr.Label.String = 'time (s)';
% title('Target sequence')


%% Plot - Target groups on FOV

this_selection1 = selection1;
this_selection2 = selection2;

F2 = figure;
hold on
scatter(med2_all,med1_all,'.','MarkerEdgeColor',[0.8,0.8,0.8])
l(1) = scatter(f.med2(this_selection1),f.med1(this_selection1),'d','MarkerEdgeColor','k','MarkerFaceColor',p.groupCols(1,:));
l(2) = scatter(f.med2(this_selection2),f.med1(this_selection2),'s','MarkerEdgeColor','k','MarkerFaceColor',p.groupCols(2,:));
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
legend(l,p.groupNames)
title('Target groups on FOV')


%% Plot - Counterbalancing of target groups

nrows = 4;
ncols = 5;
these_labels = {p.groupNames(1),p.groupNames(2)};
plot_positions = [1,6,2,7,12,17,3,8,13,18,11,4,9,14,5,10,15];
exc = [1,1,0,0,0,0,0,2,2,0,0,0,0,0,0,0,0];

F3 = figure;
for n=1:17
    subplot(nrows,ncols,plot_positions(n));
    if exc(n)==2
        v = violinplot([features.cellw1(:,n),features.cellw2(:,n)]*100,these_labels);
        ytickformat('percentage');
    else
        v = violinplot([features.cellw1(:,n),features.cellw2(:,n)],these_labels);
    end
    v(1).ViolinColor = p.groupCols(1,:); v(2).ViolinColor = p.groupCols(2,:); v(1).BoxColor = 'k'; v(2).BoxColor = 'k';
    ylabel(features.names(n));
    title(['ranksum: p=',num2str(features.pval(n,1),2),newline,'kstest: p=',num2str(features.pval(n,2),2)],'FontSize',10,'FontWeight','normal');
    
    if exc(n)==1
        ylim([0,floor(info.FOVsize)*info.numPixels/info.FOVsize]);
        yticks([0:(floor(info.FOVsize)*info.numPixels/info.FOVsize)/4:floor(info.FOVsize)*info.numPixels/info.FOVsize])
        yticklabels(strsplit(num2str([0:floor(info.FOVsize)/4:floor(info.FOVsize)])))
    end
    
end

suptitle('Counterbalancing of target groups');


%% Plot - Target clusters on FOV

nrows = 4;
ncols = 5;

temp1 = clustering(:,grouping(:,1));
temp2 = clustering(:,grouping(:,2));
temp1_2 = galvo(:,grouping(:,1));
temp2_2 = galvo(:,grouping(:,2));

F4 = figure;
for n=1:p.numClustersPerGroup
    subplot(nrows,ncols,n)
    hold on
    scatter(med2_all,med1_all,'.','MarkerEdgeColor',[0.8,0.8,0.8])
    scatter(f.med2(temp1(:,n)),f.med1(temp1(:,n)),'d','MarkerEdgeColor','k','MarkerFaceColor',p.groupCols(1,:))
    scatter(f.med2(temp2(:,n)),f.med1(temp2(:,n)),'s','MarkerEdgeColor','k','MarkerFaceColor',p.groupCols(2,:))
    scatter(temp1_2(2,n),temp1_2(1,n),'x','MarkerEdgeColor',p.groupCols(1,:),'LineWidth',1)
    scatter(temp2_2(2,n),temp2_2(1,n),'x','MarkerEdgeColor',p.groupCols(2,:),'LineWidth',1)
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
    %xticks([0:info.numPixels/4:info.numPixels])
    xticks([])
    xticklabels([])
    ylim([0,info.numPixels])
    %yticks([0:info.numPixels/4:info.numPixels])
    yticks([])
    yticklabels([])
    xlabel(['Cluster ',num2str(n)])
end
han=axes(F4,'visible','off'); han.Title.Visible='on'; han.XLabel.Visible='on'; han.YLabel.Visible='on';
xlabel(han,'medial <-> lateral');
ylabel(han,'posterior <-> anterior');
suptitle('Target clusters on FOV');


%% Plot - Counterbalancing of target clusters

nrows = 4;
ncols = 5;
these_labels = categorical({char(p.groupNames(1)),char(p.groupNames(2))});
these_labels = reordercats(these_labels,{char(p.groupNames(1)),char(p.groupNames(2))});
plot_positions = [1,6,2,7,12,17,3,8,13,18,11,4,9,14,5,10,15,16];
exc = [1,1,0,3,3,3,0,2,2,4,0,0,0,0,3,3,3,0];

F5 = figure;
for n=1:18
    subplot(nrows,ncols,plot_positions(n));
    hold on
    if exc(n)==2
        v = bar(these_labels,[mean(features.clusterw1(:,n)),mean(features.clusterw2(:,n))]*100);
        for i=1:p.numClustersPerGroup
            plot(these_labels,[features.clusterw1(i,n),features.clusterw2(i,n)]*100,'-k')
        end
        ytickformat('percentage');
    else
        v = bar(these_labels,[mean(features.clusterw1(:,n)),mean(features.clusterw2(:,n))]);
        for i=1:p.numClustersPerGroup
            plot(these_labels,[features.clusterw1(i,n),features.clusterw2(i,n)],'-k')
        end
    end
    hold off
    v.FaceColor = 'flat'; v.CData(1,:) = p.groupCols(1,:); v.CData(2,:) = p.groupCols(2,:);
    ylabel(features.names(n));
    title(['signrank: p=',num2str(features.pval(n,3),2),newline,'kstest: p=',num2str(features.pval(n,4),2)],'FontSize',10,'FontWeight','normal');
    
    if exc(n)==1
        ylim([0,floor(info.FOVsize)*info.numPixels/info.FOVsize]);
        yticks([0:(floor(info.FOVsize)*info.numPixels/info.FOVsize)/4:floor(info.FOVsize)*info.numPixels/info.FOVsize])
        yticklabels(strsplit(num2str([0:floor(info.FOVsize)/4:floor(info.FOVsize)])))
    elseif exc(n)==4
        ylim([floor(min(min([features.clusterw1(:,n),features.clusterw2(:,n)]))*10)/10,ceil(max(max([features.clusterw1(:,n),features.clusterw2(:,n)]))*10)/10]);
    elseif exc(n)==3
        ylim([floor(min(min([features.clusterw1(:,n),features.clusterw2(:,n)]))*100)/100,ceil(max(max([features.clusterw1(:,n),features.clusterw2(:,n)]))*100)/100]);
    elseif exc(n)==2
        ylim([floor(min(min([features.clusterw1(:,n),features.clusterw2(:,n)]*100))),ceil(max(max([features.clusterw1(:,n),features.clusterw2(:,n)]*100)))]);
    else
        ylim([floor(min(min([features.clusterw1(:,n),features.clusterw2(:,n)]))),ceil(max(max([features.clusterw1(:,n),features.clusterw2(:,n)])))]);
    end
    
end

suptitle('Counterbalancing of target clusters');


%% Save results

disp(['- Saving results.'])

trg.info = info;
trg.path = path;
trg.p = p;
trg.xml = xml;
trg.s2p = s2p;
trg.iscells = iscells;
trg.seq = seq;
trg.f = f;
trg.pwd = pwd;
trg.pwc = pwc;
trg.cmp_within = cmp_within;
trg.cmp_across = cmp_across;
trg.clustering = clustering;
trg.grouping = grouping;
trg.cluster_centres = cluster_centres;
trg.galvo = galvo;
trg.sequenceClusters = sequenceClusters;
trg.sequenceOrder = sequenceOrder;
trg.selection1 = selection1;
trg.selection2 = selection2;
trg.oris = oris;
trg.features = features;

if ~exist(path.savePath,'dir')
    mkdir(path.savePath);
end

temp = datestr(now,'YYYYmmdd_HHMMSS');
writeTargetXML(trg,[path.savePath,['Targeting_',info.animal,'_',temp,'.xml']]);
save([path.savePath,['Targeting_',info.animal,'_',temp,'.mat']],'trg','-v7.3');
save([path.savePath,['oris_',info.animal,'_',temp,'.mat']],'oris');
%savefig(F1,[path.savePath,['Targeting_',info.animal,'_',temp,'_F1.fig']]);
savefig(F2,[path.savePath,['Targeting_',info.animal,'_',temp,'_F2.fig']]);
savefig(F3,[path.savePath,['Targeting_',info.animal,'_',temp,'_F3.fig']]);
savefig(F4,[path.savePath,['Targeting_',info.animal,'_',temp,'_F4.fig']]);
savefig(F5,[path.savePath,['Targeting_',info.animal,'_',temp,'_F5.fig']]);

temp = toc;
disp(['- Done. Elapsed time: ',num2str(floor(temp/60)),' minutes and ',num2str(floor(rem(temp,60))),' seconds.']);
