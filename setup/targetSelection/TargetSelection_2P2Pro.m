%% Select data and parameters

% suite2p QuickCuration:
% 1) include all ROIs
% 2) compact: remove top n -> removes scanning artefacts (threshold: ~1.1)
% 3) skew: add top n back in, then remove bottom n -> selects for cells with sparse transients (threshold: ~1)
% 4) draw selection to remove big bad areas

clc;
clear;
tic;

info.animal             = 'E97'; %'Mufasa';
info.date               = '20210624'; %'20200520';
info.recording          = 'E97-20210624-spont'; %'Mufasa-20200520-5min';
info.rootPath           = 'B:\Data\MethodsPaper\'; %C:\Users\user\Documents\MATLAB\SniffinHippo\TargetSelection\TestDataForMethodsPaper\'; % 'G:\Data\SniffinHippo\';
% will save stuff here: e.g. '%ROOT%\2021\2021-04\2021-04-14\Cardano\Targeting\TargetSelection\'


%%% --- power parameters --- %%%

p.powerPerCell          = 3; % [mW]
p.powerPerCell_blast    = 3; % [mW]

info.calibration        = [...
    0,  1.9;...
    1,  1.9;...
    2,  2.1;...
    3,  3.1;...
    4,  5.4;...
    5,  9;...
    6,  13.7;...
    7,  19.7;...
    8,  26.8;...
    9,  35.2;...
    10, 44;...
    11, 55;...
    12, 66;...
    13, 79;...
    14, 92;...
    15, 106;...
    20, 183;...
    25, 275;...
    30, 359;...
    33, 405];
    % needs to be high enough for blast
    
%%% --------------------------- %%%


% general assignment
p.discardFOVEdges       = 50; % [um]
p.minDistanceAcross     = 8; % [um] (above 10 takes long, even above 8)

% cluster assignment
p.maxClusterExtent      = 450; % [um](below 400-450 can take very long)
p.minDistanceWithin     = 50; % [um]
p.clustAssignmentIters  = 3000;
p.maxSlippage_number    = 0.1; % fraction of cells to drop from max to reach good counterbalancing
p.maxSlippage_spec      = 0.1; % fraction of cells to drop from max to reach good counterbalancing
p.galvoDisplacement     = 50; % [um]

% groups
p.groupNames            = [ "PC-start", "PC-reward", "PCnum-high", "PCnum-middle", "PCnum-low", "PC-rand", "nonPC", "blast" ];
p.groupCols             = [ 192,202,50; 130,120,29; [linspace(136,241,3);linspace(19,143,3);linspace(80,177,3)]'; 90,188,173; 4,88,156; 150,150,150 ]/255;
p.groupOrder            = [1,2,5,6,7,3,4,8];
p.boundaries            = [ 20,60; 240,280; 100,200; 100,200; 100,200; 20,280; NaN,NaN; NaN,NaN ]; % groups can't be overlapping until incl. group 5
p.prcMiddle             = 2/3;
p.prcLow                = 0.5; % (fraction of prcMiddle -> 0.5 means 1/3 of number of cells in high)

% trial sequence
p.trialsPerBlock        = 35;
p.numBlocks             = 20;

% blast
p.do_blast              = 1;
p.numTrials_blast       = 50;
p.numTargets_blast      = 100;

% experiment conditions
info.FOVsize            = 1000; % [um]
info.numPixels          = 512;
p.smoothingSD           = 2; % [frames]
p.neuropilSubtraction   = true;

% xlm info
xml.name_root           = 'SLMPattern_';
xml.shape               = 'Ellipse';
xml.pxSpacing           = '1';
xml.durationMS          = '20'; % spiral duration [ms]
xml.iterations          = '1'; % spiral iterations
xml.prePatIdleMS        = '0';
xml.postPatIdleMS       = '0'; % time between patterns [ms]
xml.preIteIdleMS        = '0';
xml.postIteIdleMS       = '0';
xml.measurePowerMW      = '0';
xml.measurePowerMWPerUM2= '0';
xml.sequenceEpochCount  = '1';
xml.roiWidthUM          = 8; % spiral size [um]
xml.roiHeightUM         = 8; % spiral size [um]

%%% --- preparations ---

info.pixelSize          = info.FOVsize / info.numPixels;
xml.roiWidthPx          = num2str(xml.roiWidthUM / info.pixelSize);
xml.roiHeightPx         = num2str(xml.roiHeightUM / info.pixelSize);
if p.do_blast
    if p.numTargets_blast*p.powerPerCell_blast > info.calibration(end,2)
        error('Power calibration does not extend far enough.');
    end
    p.numClusters           = length(p.groupNames)-1;
    p.numClusters_withBlast = length(p.groupNames);
else
    p.numClusters           = length(p.groupNames)-1;
    p.numClusters_withBlast = length(p.groupNames)-1;
end

disp(['Target selection.'])


%% Load data

disp(['- Loading data.'])

info.year = info.date(1:4);
info.month = info.date(5:6);
info.day = info.date(7:8);
path.homeFolder = [info.rootPath,info.year,'\',info.year,'-',info.month,'\',info.year,'-',info.month,'-',info.day,'\',info.animal,'\','Targeting','\'];
path.s2pPath = [path.homeFolder,info.recording,'\suite2p\plane0\Fall.mat'];
path.iscellPath = [path.homeFolder,info.recording,'\suite2p\plane0\iscell.npy'];
path.pcaPath = [path.homeFolder,'Place_Cells.csv'];
% temp = dir([path.homeFolder,info.animal,'_',info.date,'_*_SEQ.txt']);
% path.seqPath = [path.homeFolder,temp.name];
path.savePath = [path.homeFolder,'TargetSelection\'];

% import pca file
if ispc
    pca = readtable(path.pcaPath);
else
    pca = readtable('/Users/Moritz/Documents/MATLAB/PhD/MethodsPaper/TestDataForMethodsPaper/2021/2021-04/2021-04-14/Cardano/Targeting/Place_Cells.csv');
end

% import s2p file
if ispc
    s2p = load(path.s2pPath);
else
    s2p = load('/Users/Moritz/Documents/MATLAB/PhD/MethodsPaper/TestDataForMethodsPaper/2021/2021-04/2021-04-14/Cardano/Targeting/Cardano-20210414-spont/suite2p/plane0/Fall.mat');
end

numROIs = length(s2p.stat);
if ispc
    temp = readNPY(path.iscellPath);
else
    temp = readNPY('/Users/Moritz/Documents/MATLAB/PhD/MethodsPaper/TestDataForMethodsPaper/2021/2021-04/2021-04-14/Cardano/Targeting/Cardano-20210414-spont/suite2p/plane0/iscell.npy');
end
iscells = find(temp(:,1)); %find(s2p.iscell(:,1));
if size(iscells,1)~=size(pca,1)
    error('Number of iscells differs between iscell.npy and Place_Cells.csv files')
end
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

% get place cell indices
placecells_all = iscells_all(find(~isnan(pca.Field_Peak_Location__cm_)));
placecells = intersect(placecells_all,iscells);
numPlaceCells = length(placecells);
placecells_all_iccords = find(~isnan(pca.Field_Peak_Location__cm_));
[~,temp] = setdiff(placecells_all,placecells);
placecells_iccords = setdiff(placecells_all_iccords,placecells_all_iccords(temp));


%% Process data

disp(['- Processing data.'])

if p.neuropilSubtraction
    F = s2p.F - 0.7*s2p.Fneu;
else
    F = s2p.F;
end
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
f.width = zeros(numCells,1);
f.amp = zeros(numCells,1);
f.var = zeros(numCells,1);
f.loc = zeros(numCells,1);
[~,temp] = intersect(iscells_all,iscells);
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
    f.mean(c,1) = mean(dFF(iscells(c),:)); %mean(s2p.F(iscells(c),:)-0.7*s2p.Fneu(iscells(c),:));
    f.std(c,1) = std(dFF(iscells(c),:)); %s2p.stat{1,iscells(c)}.std; % is calculated like this std(s2p.F(1,:)-0.7*s2p.Fneu(1,:))
    f.skew(c,1) = skewness(dFF(iscells(c),:)); %s2p.stat{1,iscells(c)}.skew; % is calculated like this skewness(s2p.F(1,:)-0.7*s2p.Fneu(1,:))
    
    % place cell features
    f.width(c,1) = pca.Field_Widths__cm_(temp(c));
    f.amp(c,1) = pca.Response_Amplittude(temp(c));
    try
        f.var(c,1) = pca.Response_Variability(temp(c));
    catch
        temp2 = str2num(cell2mat(pca.Response_Variability(temp(c))));
        if isempty(temp2)
            f.var(c,1) = NaN;
        else
            f.var(c,1) = temp2;
        end
    end
    f.loc(c,1) = pca.Field_Peak_Location__cm_(temp(c));
end
f_labels = fields(f);

if numPlaceCells ~= length(rmmissing(f.loc))
    error(['Some indexing problem with the place cell data.'])
end

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

pwd_z = squareform(zscore(squareform(pwd)));
temp = pwc;
temp(1:numCells+1:end) = 0;
pwc_z = squareform(zscore(squareform(temp)));

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


%% Identify cand.candidates for different place cell selections

% candidate cells
cand.candidates = zeros(10000,p.numClusters_withBlast);
for i=1:p.numClusters
    % relative to the iscells number
    if ~isnan(p.boundaries(i,1))
        temp = find((f_vec(:,end) >= p.boundaries(i,1)) & (f_vec(:,end) <= p.boundaries(i,2)));
    else
        temp = find(isnan(f_vec(:,end)));
    end
    cand.candidates(1:length(temp),i) = temp;
end
if p.do_blast
    temp = find(f.med1 < info.numPixels/2 + (p.maxClusterExtent/info.pixelSize)/2 & ...
        f.med1 > info.numPixels/2 - (p.maxClusterExtent/info.pixelSize)/2 & ...
        f.med2 < info.numPixels/2 + (p.maxClusterExtent/info.pixelSize)/2 & ...
        f.med2 > info.numPixels/2 - (p.maxClusterExtent/info.pixelSize)/2);
    cand.candidates(1:length(temp),8) = temp;
end
cand.candidates(~any(cand.candidates,2),:) = [];

% get candidate cliques for stim group 1 (PC_start) and stim group 2 (PC_reward)
cand.cliques_g1 = maximalCliques( cmp(nonzeros(cand.candidates(:,1)),nonzeros(cand.candidates(:,1))) );
temp = sum(cand.cliques_g1,1);
temp2 = floor(max(temp)*(1-p.maxSlippage_spec));
cand.candidate_cliques_g1 = find(temp >= temp2);
cand.cliques_g2 = maximalCliques( cmp(nonzeros(cand.candidates(:,2)),nonzeros(cand.candidates(:,2))) );
temp = sum(cand.cliques_g2,1);
temp2 = floor(max(temp)*(1-p.maxSlippage_spec));
cand.candidate_cliques_g2 = find(temp >= temp2);

% get combinations of compatible candidate cliques for stim group 1 (PC_start) and stim group 2 (PC_reward)
[temp1,temp2] = meshgrid(cand.candidate_cliques_g1,cand.candidate_cliques_g2);
temp = reshape(cat(2,temp1',temp2'),[],2);
cand.candidate_cliques_g1_g2 = [];
for i=1:10000 %size(temp,1)
    if all(squareform(cmp_across([cand.candidates(find(cand.cliques_g1(:,temp(i,1))),1);cand.candidates(find(cand.cliques_g2(:,temp(i,2))),2)],...
            [cand.candidates(find(cand.cliques_g1(:,temp(i,1))),1);cand.candidates(find(cand.cliques_g2(:,temp(i,2))),2)])))
        cand.candidate_cliques_g1_g2 = [cand.candidate_cliques_g1_g2; temp(i,:)];
    end
end
if isempty(cand.candidate_cliques_g1_g2)
    warning('Cliques for start and reward groups not compatible. Increase p.maxSlippage_spec.')
end

% get candidate cliques for stim group 3 (PCnum_high)
cand.cliques_g3 = maximalCliques( cmp(nonzeros(cand.candidates(:,3)),nonzeros(cand.candidates(:,3))) );
temp = sum(cand.cliques_g3,1);
temp2 = floor(max(temp)*(1-p.maxSlippage_number));
cand.candidate_cliques_g3 = find(temp >= temp2);

% get candidate cliques for stim group 8 (blast)
if p.do_blast
    cand.cliques_blast = maximalCliques( cmp_across(nonzeros(cand.candidates(:,8)),nonzeros(cand.candidates(:,8))) );
    cand.cliques_blast = cand.cliques_blast(:,1);
end

if p.do_blast
    maxCells = max([max(sum(cand.cliques_g3,1)),max(sum(cand.cliques_g1,1)),max(sum(cand.cliques_g2,1)),p.numTargets_blast]);
else
    maxCells = max([max(sum(cand.cliques_g3,1)),max(sum(cand.cliques_g1,1)),max(sum(cand.cliques_g2,1))]);
end


%% Assign cells into tentative clusters

disp(['- Identifying allowed clusterings.'])

n=0;
this_cmp = 0;
clusterings = zeros(maxCells,p.numClusters_withBlast,p.clustAssignmentIters);
while n<p.clustAssignmentIters
    try

        % initialise
        this_pool = 1:numCells;
        this_clustering = zeros(maxCells,p.numClusters_withBlast);

        % randomly draw clusters from cmp matrices
        for j=1:p.numClusters
            
            % --- Draw 7 groups with different place cell selections ---
            
            % PC_start and PC_reward clusters
            if j==1
                temp = randsample(1:size(cand.candidate_cliques_g1_g2,1),1);
                temp1 = cand.candidates(find(cand.cliques_g1(:,cand.candidate_cliques_g1_g2(temp,1))),1);
                temp2 = cand.candidates(find(cand.cliques_g2(:,cand.candidate_cliques_g1_g2(temp,2))),2);
                temp = min([length(temp1),length(temp2)]);
                this_clustering(1:temp,1) = sort(randsample(temp1,temp));
                this_clustering(1:temp,2) = sort(randsample(temp2,temp));
                this_pool = setdiff(this_pool,this_clustering(:));
            elseif j==2
                % already done in j==1
            
            % PCnum_high, PCnum_middle and PCnum_low clusters
            elseif j==3
                while this_cmp==0
                    temp = cand.candidates(find(cand.cliques_g3(:,randsample(cand.candidate_cliques_g3,1))),j);
                    if isempty(setdiff(temp,this_pool)) && all(squareform(cmp_across([nonzeros(this_clustering(:));temp],[nonzeros(this_clustering(:));temp])))
                        this_clustering(1:length(temp),j) = sort(temp);
                        this_pool = setdiff(this_pool,this_clustering(:));
                        this_cmp=1;
                    end
                end
                this_cmp=0;
            elseif j==4   
                temp = randsample(nonzeros(this_clustering(:,3)),round(length(nonzeros(this_clustering(:,3)))*p.prcMiddle));
                this_clustering(1:length(temp),j) = sort(temp);
            elseif j==5
                temp = randsample(nonzeros(this_clustering(:,4)),round(length(nonzeros(this_clustering(:,4)))*p.prcLow));
                this_clustering(1:length(temp),j) = sort(temp);
            
            elseif j==8 && p.do_blast
                % follows below
                
            % PC_all and nonPC clusters
            else
                for k=1:length(nonzeros(this_clustering(:,1)))
                        this_specific_pool = intersect(this_pool,nonzeros(cand.candidates(:,j)));

                        % draw seeds
                        if j==1 && k==1
                            this_clustering(k,j) = randsample(this_specific_pool,1);
                        elseif k==1
                            [~,~,temp] = find(this_clustering(:));
                            these_compatible_across = find( all(cmp_across(temp,:),1) );
                            this_clustering(k,j) = randsample( intersect(these_compatible_across,this_specific_pool) ,1);

                        % draw compatible cells
                        else
                            [~,~,temp] = find(this_clustering(:));
                            these_compatible_across = find( all(cmp_across(temp,:),1) );     
                            [~,~,temp2] = find(this_clustering(:,j));
                            these_compatible_within = find( all(cmp_within(temp2,:),1) );                  
                            this_clustering(k,j) = randsample( intersect(intersect(these_compatible_across,these_compatible_within),this_specific_pool) ,1);
                            if ~all(squareform(cmp(nonzeros(this_clustering(:,j)),nonzeros(this_clustering(:,j)))))
                                error(['Just drew a cluster that doesnt fulfill criteria. j, k: ',num2str(j),',',num2str(k)])
                            end
                        end
                        % remove selected cells from pool
                        this_pool = setdiff(this_pool,this_clustering(:));
                end
                % sort
                this_clustering(1:length(nonzeros(this_clustering(:,1))),j) = sort(this_clustering(1:length(nonzeros(this_clustering(:,1))),j));
            end
            
            % sanity checks
            if ~all(squareform(cmp(nonzeros(this_clustering(:,j)),nonzeros(this_clustering(:,j)))))
                error(['Just drew a cluster that isnt compatible within cluster. j: ',num2str(j),'.'])
            end
        end
        if ~all(squareform(cmp_across(unique(nonzeros(this_clustering(:))),unique(nonzeros(this_clustering(:))))))
            error(['Just drew a clustering that isnt compatible across clusters.'])
        end
        
        n=n+1;
        if mod(n,100)==0
            disp(['--- ',num2str(n),' allowed clusterings identified.'])
        end
        clusterings(:,:,n) = this_clustering;
    catch
        disp('Jumped into catch')
        if ~all(squareform(cmp(nonzeros(this_clustering(:,j)),nonzeros(this_clustering(:,j)))))
            error(['Catch: Just drew a cluster that doesnt fulfill criteria. j, k: ',num2str(j),',',num2str(k)])
        end
    end
end
if p.do_blast
    for i=1:p.clustAssignmentIters
        temp = randsample(cand.candidates(find(cand.cliques_blast),8),p.numTargets_blast);
        clusterings(1:length(temp),8,i) = sort(temp);
    end
end


%% Generate feature list and select best matched clustering

disp(['- Selecting best matched clustering.'])

% generate cell-wise feature lists
features.names = ["anteroposterior position (\mum)","mediolateral position (\mum)",...
    "number of pixels","compactness","footprint","aspect ratio",...
    "baseline fluorescence intensity","mean activity (\DeltaF/F)","s.d. of activity (\DeltaF/F)","skewness of activity",...
    "PF width (cm)","PF amplitude","PF variability","PF peak location (cm)",...
    "distance to cluster centre (\mum)","mean pw dist - FOV (\mum)","mean pw dist - exp (\mum)","mean pw dist - cluster (\mum)",...
    "mean pw corr - FOV","mean pw corr - exp","mean pw corr - cluster",...
    "cluster extent (\mum)"];

max_minpval = 0;
for i=1:p.clustAssignmentIters
    
    temp = clusterings(:,:,i);
    this_selection_unsorted = nonzeros(temp(:));
    this_selection = sort(nonzeros(temp(:)));
    these_numCellsPerCluster = zeros(1,p.numClusters);
    for j=1:p.numClusters
        these_numCellsPerCluster(j) = length(nonzeros(temp(:,j)));
    end
    these_cumNumCellsPerCluster = cumsum(these_numCellsPerCluster);
    
    % get cluster centres
    these_cluster_centres = zeros(2,p.numClusters);
    for j=1:p.numClusters
        these_cluster_centres(1,j) = (min(f.med1(nonzeros(clusterings(:,j,i)))) + max(f.med1(nonzeros(clusterings(:,j,i))))) / 2;
        these_cluster_centres(2,j) = (min(f.med2(nonzeros(clusterings(:,j,i)))) + max(f.med2(nonzeros(clusterings(:,j,i))))) / 2;
    end
    
    % generate cell-wise feature lists
    these_features_cellw = zeros(length(this_selection_unsorted),21);
    for j=1:length(f_labels)
        these_features_cellw(:,j) = f.(f_labels{j})(this_selection_unsorted);
    end
    for j=1:p.numClusters
        if j==1
            these_features_cellw(1:these_cumNumCellsPerCluster(j),15) = pdist2([these_features_cellw(1:these_cumNumCellsPerCluster(j),1),these_features_cellw(1:these_cumNumCellsPerCluster(j),2)],these_cluster_centres(:,j)','euclidean');
        else
            these_features_cellw(these_cumNumCellsPerCluster(j-1)+1:these_cumNumCellsPerCluster(j),15) = pdist2([these_features_cellw(these_cumNumCellsPerCluster(j-1)+1:these_cumNumCellsPerCluster(j),1),these_features_cellw(these_cumNumCellsPerCluster(j-1)+1:these_cumNumCellsPerCluster(j),2)],these_cluster_centres(:,j)','euclidean');
        end
    end 
    temp = pwd(this_selection_unsorted,:);
    temp(find(temp==0))=NaN;
    these_features_cellw(:,16) = nanmean(temp,2); % mean pwd to any cell in FOV
    temp = pwd(this_selection_unsorted,this_selection_unsorted);
    temp(1:length(this_selection_unsorted)+1:end) = diag(NaN(length(this_selection_unsorted)));
    these_features_cellw(:,17) = nanmean(temp,2); % mean pwd to any other stimmed cell
    for j=1:p.numClusters
        if j==1
            these_features_cellw(1:these_cumNumCellsPerCluster(j),18) = nanmean(temp(1:these_cumNumCellsPerCluster(j),1:these_cumNumCellsPerCluster(j)),2); % mean pwd to any cell in same cluster            
        else
            these_features_cellw(these_cumNumCellsPerCluster(j-1)+1:these_cumNumCellsPerCluster(j),18) = nanmean(temp(these_cumNumCellsPerCluster(j-1)+1:these_cumNumCellsPerCluster(j),these_cumNumCellsPerCluster(j-1)+1:these_cumNumCellsPerCluster(j)),2); % mean pwd to any cell in same cluster            
        end
    end
    temp = pwc(this_selection_unsorted,:);
    temp(find(temp==0))=NaN;
    these_features_cellw(:,19) = nanmean(temp,2); % mean pwc with any cell in FOV
    temp = pwc(this_selection_unsorted,this_selection_unsorted);
    temp(1:length(this_selection_unsorted)+1:end) = diag(NaN(length(this_selection_unsorted)));
    these_features_cellw(:,20) = nanmean(temp,2); % mean pwc with any other stimmed cell
    for j=1:p.numClusters
        if j==1
            these_features_cellw(1:these_cumNumCellsPerCluster(j),21) = nanmean(temp(1:these_cumNumCellsPerCluster(j),1:these_cumNumCellsPerCluster(j)),2); % mean pwc with any cell in same cluster
        else
            these_features_cellw(these_cumNumCellsPerCluster(j-1)+1:these_cumNumCellsPerCluster(j),21) = nanmean(temp(these_cumNumCellsPerCluster(j-1)+1:these_cumNumCellsPerCluster(j),these_cumNumCellsPerCluster(j-1)+1:these_cumNumCellsPerCluster(j)),2); % mean pwc with any cell in same cluster
        end
    end
    
    % generate cluster-wise feature lists
    these_features_clusterw = zeros(p.numClusters,size(these_features_cellw,2)+1);
    for j=1:p.numClusters
        if j==1
            these_features_clusterw(j,1:end-1) = nanmean(these_features_cellw(1:these_cumNumCellsPerCluster(j),:),1);
        else
            these_features_clusterw(j,1:end-1) = nanmean(these_features_cellw(these_cumNumCellsPerCluster(j-1)+1:these_cumNumCellsPerCluster(j),:),1);
        end
    end
    temp = pwd(this_selection_unsorted,this_selection_unsorted);
    temp(1:length(this_selection_unsorted)+1:end) = diag(NaN(length(this_selection_unsorted)));
    for j=1:p.numClusters
        if j==1
            these_features_clusterw(j,end) = nanmax(nanmax(temp(1:these_cumNumCellsPerCluster(j),1:these_cumNumCellsPerCluster(j)),[],2),[],1); % max pwd to any cell in same cluster
        else
            these_features_clusterw(j,end) = nanmax(nanmax(temp(these_cumNumCellsPerCluster(j-1)+1:these_cumNumCellsPerCluster(j),these_cumNumCellsPerCluster(j-1)+1:these_cumNumCellsPerCluster(j)),[],2),[],1); % max pwd to any cell in same cluster
        end
    end
       
    % calculate p values - group comparisons and cluster comparisons
    pval=NaN(21,1);
    these_features = nan([size(this_clustering),21]);
    for k=1:21
        this_feature = nan(size(this_clustering));
        
        for j=1:p.numClusters
            if j==1
                temp = these_features_cellw(1:these_cumNumCellsPerCluster(j),k);
            else
                temp = these_features_cellw(these_cumNumCellsPerCluster(j-1)+1:these_cumNumCellsPerCluster(j),k);
            end
            
            these_features(1:length(temp),j,k) = temp;

            % skipping features that shouldn't be used for counterbalancing
            if  ismember(k,[1,2]) ||...
                (ismember(k,[14,19:21]) && ismember(j,[1,2,6,7])) ||...
                (ismember(k,16:18) && ismember(j,[3,4,5])) 
                temp = nan(length(temp),1);
            end

            this_feature(1:length(temp),j) = temp;
        end
        
        % calculate p value
        pval(k) = kruskalwallis(this_feature,categorical(p.groupNames),'off');
        features.names(k);
    end
    
    % update current best option
    if nanmin(pval(:)) > max_minpval
        disp(['--- current max. min. p-value: ',num2str(nanmin(pval(:)),2)])
        max_minpval = nanmin(pval(:));
        winner = i;
        clustering = clusterings(:,:,i);
        selection = this_selection;
        selection_unsorted = this_selection_unsorted;
        
        features.features = these_features;
        features.cluster_centres = these_cluster_centres;
        features.selection_unsorted = this_selection_unsorted;
        features.cellw = these_features_cellw;
        features.clusterw = these_features_clusterw;
        features.pval = pval;
    end
end
if max_minpval < 0.05
    warning('No assignment without significant differences between the groups was found.')
end


%% Finalise output (clustering, grouping, cluster_centres, galvo, sequences)

disp(['- Finalising output.'])

% various
clustering(~any(clustering,2),:) = [];

% assign galvo positions
cluster_centres = zeros(2,p.numClusters_withBlast);
galvo = zeros(2,p.numClusters_withBlast);
for i=1:p.numClusters_withBlast
    cluster_centres(1,i) = (min(f.med1(nonzeros(clustering(:,i)))) + max(f.med1(nonzeros(clustering(:,i))))) / 2;
    cluster_centres(2,i) = (min(f.med2(nonzeros(clustering(:,i)))) + max(f.med2(nonzeros(clustering(:,i))))) / 2;
    
    temp2 = zeros(info.numPixels,info.numPixels);
    temp = zeros(numel(temp2),2);
    [temp(:,1),temp(:,2)] = find(temp2==0);
    
    temp2 = reshape(min(pdist2(temp,[f.med1(nonzeros(clustering(:,i))),f.med2(nonzeros(clustering(:,i)))]),[],2) > (p.galvoDisplacement/info.pixelSize),info.numPixels,info.numPixels);
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

% generate trial sequence
seq = [];
for i=1:p.numBlocks
    this_stim = [];
    for j=1:p.trialsPerBlock/p.numClusters
        this_stim = cat(2,this_stim,1:p.numClusters);
    end
    this_stim = this_stim(randperm(p.trialsPerBlock));
    seq = [seq,this_stim];
end
if p.do_blast
    seq = [seq,ones(1,p.numTrials_blast)*8];
end

% prepare stim order for STA movie maker
oris = seq;

% assign stim powers
temp1 = 0:0.1:max(info.calibration(:,1));
temp2 = interp1(info.calibration(:,1),info.calibration(:,2),temp1);
xml.power = {}; % power per cluster [% of 1V output]
for i=1:p.numClusters_withBlast
    this_num = length(nonzeros(clustering(:,i)));
    if p.do_blast && i==8
        this_mW = this_num*p.powerPerCell_blast;
    else
        this_mW = this_num*p.powerPerCell;
    end
    [a,temp] = min(abs(temp2-this_mW));
    temp = temp1(temp);
    xml.power{i} = num2str(temp,2);
end


%% --- Plot results ---

disp(['- Plotting results.'])

default_figure();
close;


%% Plot - Target groups on FOV

this_selection = selection;

F1 = figure;
hold on
for i=1:p.numClusters
    scatter(med2_all,med1_all,'.','MarkerEdgeColor',[0.8,0.8,0.8])
    l(i) = scatter(f.med2(nonzeros(clustering(:,i))),f.med1(nonzeros(clustering(:,i))),'d','MarkerEdgeColor','k','MarkerFaceColor',p.groupCols(i,:));
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
legend(l,p.groupNames(1:p.numClusters))
title('Target clusters on FOV')


%% Plot - Counterbalancing of target groups

nrows = 4;
ncols = 6;
these_labels = {'1','2','3','4','5','6','7'}; %{p.groupNames(1),p.groupNames(2),p.groupNames(3),p.groupNames(4),p.groupNames(5),p.groupNames(6),p.groupNames(7)};
plot_positions = [1,6+1,2,7+1,12+2,17+3,3,8+1,13+2,18+3, 6,12,18,24, 11+2,4,9+1,14+2,5,10+1,15+2];
exc = [1,1,0,0,0,0,0,2,2,0, 0,0,0,0, 0,0,0,0,0,0,0];

F2 = figure;
for n=1:21
    subplot(nrows,ncols,plot_positions(n));
    if exc(n)==2
        v = violinplot([features.features(:,find(p.groupOrder==1),n),features.features(:,find(p.groupOrder==2),n),features.features(:,find(p.groupOrder==3),n),features.features(:,find(p.groupOrder==4),n),features.features(:,find(p.groupOrder==5),n),features.features(:,find(p.groupOrder==6),n),features.features(:,find(p.groupOrder==7),n)]*100,these_labels);
        ytickformat('percentage');
    else
        v = violinplot([features.features(:,find(p.groupOrder==1),n),features.features(:,find(p.groupOrder==2),n),features.features(:,find(p.groupOrder==3),n),features.features(:,find(p.groupOrder==4),n),features.features(:,find(p.groupOrder==5),n),features.features(:,find(p.groupOrder==6),n),features.features(:,find(p.groupOrder==7),n)],these_labels);
    end
    for i=1:p.numClusters
        temp = find(p.groupOrder==i);
        if ismember(n,[1,2]) ||...
                (ismember(n,[14,19:21]) && ismember(temp,[1,2,6,7])) ||...
                (ismember(n,16:18) && ismember(temp,[3,4,5])) 
            v(i).ViolinColor = [0.5,0.5,0.5];
        else
            v(i).ViolinColor = p.groupCols(find(p.groupOrder==i),:);
        end
        v(i).BoxColor = 'k';
    end
    ylabel(features.names(n));
    title(['KW-test: p=',num2str(features.pval(n),2)],'FontSize',10,'FontWeight','normal');
    
    if exc(n)==1
        ylim([0,floor(info.FOVsize)*info.numPixels/info.FOVsize]);
        yticks([0:(floor(info.FOVsize)*info.numPixels/info.FOVsize)/4:floor(info.FOVsize)*info.numPixels/info.FOVsize])
        yticklabels(strsplit(num2str([0:floor(info.FOVsize)/4:floor(info.FOVsize)])))
    end
    
end

suptitle('Counterbalancing of target groups');


%% Plot - Target clusters on FOV

nrows = 2;
ncols = 4;

F3 = figure;
for n=1:p.numClusters_withBlast
    subplot(nrows,ncols,p.groupOrder(n))
    hold on
    scatter(med2_all,med1_all,'.','MarkerEdgeColor',[0.8,0.8,0.8])
    scatter(f.med2(nonzeros(clustering(:,n))),f.med1(nonzeros(clustering(:,n))),'d','MarkerEdgeColor','k','MarkerFaceColor',p.groupCols(n,:))
    scatter(galvo(2,n),galvo(1,n),'x','MarkerEdgeColor',p.groupCols(n,:),'LineWidth',1)
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
    xlabel([num2str(n),' - ',char(p.groupNames(n))])
end
han=axes(F3,'visible','off'); han.Title.Visible='on'; han.XLabel.Visible='on'; han.YLabel.Visible='on';
xlabel(han,'medial <-> lateral');
ylabel(han,'posterior <-> anterior');
suptitle('Target clusters on FOV');


%% Plot - Target clusters on FOV

F4 = figure;
hold on
for i=1:p.numClusters
    temp = rmmissing(features.features(:,i,14));
    if ~isempty(temp)
        xline(p.boundaries(i,1),'LineStyle',':');
        xline(p.boundaries(i,2),'LineStyle',':');
        plot(temp,ones(length(temp),1)*p.groupOrder(i),'d','MarkerEdgeColor','k','MarkerFaceColor',p.groupCols(i,:))
    end 
end
hold off
xlim([0,400])
ylim([0,p.numClusters+1])
xlabel('Place field peak position (cm)');
yticks([1,2,3,4,5,6,7]);
yticklabels(p.groupNames([1,2,6,7,3,4,5]));
set(gca, 'YDir','reverse');

    
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

trg.cluster_centres = cluster_centres;
trg.galvo = galvo;
trg.selection = selection;
trg.selection_unsorted = selection_unsorted;
trg.oris = oris;
trg.features = features;
trg.pca = pca;
trg.cand = cand;

if ~exist(path.savePath,'dir')
    mkdir(path.savePath);
end

temp = datestr(now,'YYYYmmdd_HHMMSS');
writeTargetXML_2P2Pro(trg,[path.savePath,['Targeting_',info.animal,'_',temp,'.xml']]);
save([path.savePath,['Targeting_',info.animal,'_',temp,'.mat']],'trg','-v7.3');
save([path.savePath,['oris_',info.animal,'_',temp,'.mat']],'oris');
savefig(F1,[path.savePath,['Targeting_',info.animal,'_',temp,'_F1.fig']]);
savefig(F2,[path.savePath,['Targeting_',info.animal,'_',temp,'_F2.fig']]);
savefig(F3,[path.savePath,['Targeting_',info.animal,'_',temp,'_F3.fig']]);
savefig(F4,[path.savePath,['Targeting_',info.animal,'_',temp,'_F4.fig']]);

temp = toc;
disp(['- Done. Elapsed time: ',num2str(floor(temp/60)),' minutes and ',num2str(floor(rem(temp,60))),' seconds.']);

