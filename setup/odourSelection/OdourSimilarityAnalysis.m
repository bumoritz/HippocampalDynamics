%% Load data

temp = readmatrix('C:\Users\Moritz\Desktop\Odourant Similarity Analysis\Descriptors.csv');
data_raw = temp(2:end,:);

% temp2 = readmatrix('C:\Users\Moritz\Desktop\OdourantDescriptors.xlsx');
% data_raw = temp2(:,2:end);

labels = {'Eucalyptol','IsoamylAcetate','Pinene','EthylButyrate','Heptanone','Octanal','Terpinene','IsobutylPropionate','BenzylAcetate','Benzaldehyde',...
    'Nonanone','Dimethoxybenzene','MethylBenzoate','Allylanisole','Acetophenone','EthylValerate','Limonene','Heptanol','Pentanol','Fenchone',...
    'Valeraldehyde','PropylAcetate','Cineole','Octenol','HexylAcetate','Decenal','Citronellol','Eugenol','Guaiacol','IsobutylAcetate',...
    'EthylAcetate','MethylAcetate','Geraniol','Dimethyloctadienol','Hexanone','Pentadione','EthylPropionate','Pentanone','ButylFormate','PropylFormate',...
    'PhenethylAlcohol','HexanoicAcid','AceticAcid','Prenol','Methylthiopropanol','MethylButyrate'};


%% Preparations

[numOdours, numDescriptors] = size(data_raw);

data_zscored = zscore(data_raw,0,1);


%% PCA

% [pca_coeff,pca_score,pca_latent,pca_tsquared,pca_explained,pca_mu] = pca(data_zscored);
% 
% figure;
% scatter(pca_score(:,1),pca_score(:,2))

%% Calculate pair-wise correlations

data = data_zscored; % data_zscored, pca_score(:,1:6)

correlations = zeros(numOdours);
for i=1:numOdours
    for j=1:numOdours
        correlations(i,j) = corr(data(i,:)',data(j,:)','Type','Pearson','Rows','Complete');
    end
end
distances = 1-correlations;


%% Heatmap of all odourants

figure;
heatmap(distances);
colormap('pink')
title('chemical distance')


%% Heatmap of all odourants with hierarchical clustering (old)

cm = struct('GroupNumber',{27,31,25,21,24,33,29,18,32},...
    'Annotation',{'woody','almond','waxy','glue','vomit','fruity','fruity','Christmas','piss'},...
    'Color',{rgb('Brown'),rgb('Orange'),rgb('Yellow'),rgb('Purple'),rgb('Blue'),rgb('Red'),rgb('Red'),rgb('DarkOrange'),rgb('DarkBlue')});

cgo = clustergram(distances,...
    'Colormap','pink','Symmetric',false,...
    'ColumnLabels',labels,'RowLabels',labels);
set(cgo,'ColumnGroupMarker',cm,'RowGroupMarker',cm);


%% Heatmap of all odourants with hierarchical clustering (new)

% cm = struct('GroupNumber',{27,31,25,21,24,33,29,18,32},...
%     'Annotation',{'woody','almond','waxy','glue','vomit','fruity','fruity','Christmas','piss'},...
%     'Color',{rgb('Brown'),rgb('Orange'),rgb('Yellow'),rgb('Purple'),rgb('Blue'),rgb('Red'),rgb('Red'),rgb('DarkOrange'),rgb('DarkBlue')});

cm = struct('GroupNumber',{34,17,26,31,32,33,34,35,36},...%numOdours-2*k:numOdours-k-1,...
    'Annotation',{'9','8','7','6','5','4','3','2','1'},...
    'Color',{rgb('Brown'),rgb('Orange'),rgb('Yellow'),rgb('Purple'),rgb('Blue'),rgb('Red'),rgb('DarkOrange'),rgb('DarkBlue'),rgb('Green')});

cgo = clustergram(distances,...
    'Colormap','pink','Symmetric',false,...
    'ColumnLabels',labels,'RowLabels',labels);
set(cgo,'ColumnGroupMarker',cm,'RowGroupMarker',cm);


%% Identify smallest distances

neighbours = distances<0.25;
[temp,temp2] = find(neighbours);
neighbours = [temp,temp2];

% distance < 0.2
% 1 (Eucalyptol)        - 23 (Cineole)              
% 7 (Terpinene)         - 17 (Limonene)
% 19 (Pentanol)         - 21 (Valeraldehyde)
% 20 (Fenchone)         - 23 (Cineole)
% 21 (Valeraldehyde)    - 40 (PropylFormate)
% 22 (PropylAcetate)    - 31 (EthylAcetate)         !
% 22 (PropylAcetate)    - 37 (EthylPropionate)
% 27 (Citronellol)      - 33 (Geraniol)
% 31 (EthylAcetate)     - 32 (MethylAcetate)        !
% 31 (EthylAcetate)     - 38 (Pentanone)
% 33 (Geraniol)         - 34 (Dimethyloctadienol)   
% 39 (ButylFormate)     - 40 (PropylFormate)        !

% distance < 0.2
% 1 (Eucalyptol)        - 20 (Fenchone)  
% 14 (Allylanisole)     - 28 (Eugenol)  
% 21 (Valeraldehyde)    - 38 (Pentanone)
% 22 (PropylAcetate)    - 35 (Hexanone)
% 27 (Citronellol)      - 34 (Dimethyloctadienol)
% 31 (EthylAcetate)     - 37 (EthylPropionate)      !
% 31 (EthylAcetate)     - 38 (Pentanone)
% 32 (MethylAcetate)    - 38 (Pentanone)
% 27 (Citronellol)      - 34 (Dimethyloctadienol)
% 35 (Hexanone)         - 38 (Pentanone)            !


%% Heatmap of first odourant candidates

selection_firstOdCand = [31,32, 39,40]; %[22,31, 31,32, 31,37, 35,38, 39,40];

correlations_firstOdCand = zeros(length(selection_firstOdCand));
for i=selection_firstOdCand
    for j=selection_firstOdCand
        correlations_firstOdCand(i,j) = corr(data(i,:)',data(j,:)','Type','Pearson','Rows','Complete');
    end
end
distances_firstOdCand = 1-correlations_firstOdCand;

figure;
heatmap(distances_firstOdCand);
colormap('pink')
title('chemical distance')


%% k-means clustering

k = 9;
[kmeans_idx,kmeans_C,kmeans_sumd,kmeans_D] = kmeans(data,k,'Replicates',1000,'Options',statset('Display','final'));

%% UMPA embedding

[reduction,umap,clusterIdentifiers,extras]=run_umap(data_zscored,'n_components',2,'metric','correlation');

%% Work in progress

A1_label = 'MethylAcetate';
B1_label = 'BenzylAcetate';
X1_label = 'EthylAcetate';
Y1_label = 'Pinene';

A2_label = 'PropylFormate';
B2_label = 'Guaiacol'
X2_label = 'ButylFormate';
Y2_label = 'Acetophenone';

% A3_label = 
% B3_label =
% X3_label = 
% Y3_label =

A1 = 32;
B1 = 9;
X1 = 31;
Y1 = 4;

A2 = 40;
B2 = 29;
X2 = 39;
Y2 = 15;

% A3 = 
% B3 =
% X3 = 
% Y3 =

data = data_zscored([A1,X1,B1,Y1,A2,X2,B2,Y2],:); % data_zscored, pca_score(:,1:6)

correlations = zeros(8);
for i=1:8
    for j=1:8
        correlations(i,j) = corr(data(i,:)',data(j,:)','Type','Pearson','Rows','Complete');
    end
end
distances = 1-correlations;

figure;
heatmap(distances);
colormap('pink')
title('chemical distance')


%% Maximise differences within and across set

firstOdours = [40,39,35,38,46,37];
excludedOdours = [5,6,7,11,12,14,18,20,23,24,25,26,27,28,33,34,41,42,45];
numSelectedOdours = 12;

% zscore after selecting 27 compounds
data_raw_woExluded = data_raw;
data_raw_woExluded(excludedOdours,:) = NaN;
data_woExluded_zscored = nanzscore(data_raw_woExluded,0,1);

data = data_woExluded_zscored; %data_zscored;

% list all possible combinations
numAdditionalOdours = numSelectedOdours - length(firstOdours);
possibleAdditionalOdours = setdiff(1:numOdours,[firstOdours,excludedOdours]);
combs.combinationList = nchoosek(possibleAdditionalOdours,numAdditionalOdours);
combs.numCombinations = size(combs.combinationList,1);

% go through every possible combination
combs.combinationMeanDistance_all = zeros(combs.numCombinations,1);
combs.combinationMeanDistance_woDiagonal = zeros(combs.numCombinations,1);
combs.combinationMeanDistance_woFirstOdours = zeros(combs.numCombinations,1);
combs.combinationMinDistance_all = zeros(combs.numCombinations,1);
combs.combinationMinDistance_woDiagonal = zeros(combs.numCombinations,1);
combs.combinationMinDistance_woFirstOdours = zeros(combs.numCombinations,1);
for c=1:combs.numCombinations

    this_selection = combs.combinationList(c,:);
    these_odours = [firstOdours,this_selection];
    this_data = data(these_odours,:);

    % calculate correlation matrix
    these_correlations = zeros(numSelectedOdours);
    for i=1:numSelectedOdours
        for j=1:numSelectedOdours
            these_correlations(i,j) = corr(this_data(i,:)',this_data(j,:)','Type','Pearson','Rows','Complete');
        end
    end
    these_distances = 1-these_correlations;
    
%     figure;
%     heatmap(these_distances);
%     colormap('pink');
%     title('chemical distance');

    % remove irrelevant entries
    these_distances_woDiagonal = these_distances;
    these_distances_woDiagonal(1:numSelectedOdours+1:end) = diag(NaN(numSelectedOdours));
    these_distances_woFirstOdours = these_distances_woDiagonal;
    for k=0:length(firstOdours)/2-1
        these_distances_woFirstOdours(2*k+1,2*k+2) = NaN;
        these_distances_woFirstOdours(2*k+2,2*k+1) = NaN;
    end

    combs.combinationMeanDistance_all(c) = nanmean(these_distances(:));
    combs.combinationMeanDistance_woDiagonal(c) = nanmean(these_distances_woDiagonal(:));
    combs.combinationMeanDistance_woFirstOdours(c) = nanmean(these_distances_woFirstOdours(:));
    
    combs.combinationMinDistance_all(c) = nanmin(these_distances(:));
    combs.combinationMinDistance_woDiagonal(c) = nanmin(these_distances_woDiagonal(:));
    combs.combinationMinDistance_woFirstOdours(c) = nanmin(these_distances_woFirstOdours(:));
end


%% Test winning combinations
% [a,b]=max(combs.combinationMeanDistance_woFirstOdours(1:45196)); -> 1.1648 at 3636; 12 odours -> 1.0942 at 264
% [a,b]=max(combs.combinationMinDistance_woFirstOdours(1:45196)); -> 0.2986 at 40; 12 odours -> 0.5626 at 7575

selectedCombination = 7575;

firstOdours = [40,39,35,38,46,37];

this_selection = combs.combinationList(selectedCombination,:);
these_odours = [firstOdours(1:2),this_selection(1:2),firstOdours(3:4),this_selection(3:4),firstOdours(5:6),this_selection(5:6),]; %temp
this_data = data(these_odours,:);

% calculate correlation matrix
these_correlations = zeros(numSelectedOdours);
for i=1:numSelectedOdours
    for j=1:numSelectedOdours
        these_correlations(i,j) = corr(this_data(i,:)',this_data(j,:)','Type','Pearson','Rows','Complete');
    end
end
these_distances = 1-these_correlations;
    
% plot heatmap
figure;
heatmap(these_distances);
colormap('pink');
caxis([0,0.8])
%xticklabels({'MethylAcetate','EthylAcetate','Eucalyptol','Terpinene','PropylFormate','ButylFormate','Nonanone','Eugenol'})

set(gca,'XData',{'Ai','Xi','Bi','Yi','Aii','Xii','Bii','Yii','Aiii','Xiii','Biii','Yiii'})

% [1,2,3,9,10,29]
%set(gca,'YData',{'PropylFormate','ButylFormate','Eucalyptol','IsoamylAcetate','Hexanone','Pentanone','Pinene','BenzylAcetate','MethylButyrate','EthylPropionate','Benzaldehyde','Guaiacol'})

% [1,4,9,10,30,36]
set(gca,'YData',{'PropylFormate','ButylFormate','Eucalyptol','EthylButyrate','Hexanone','Pentanone','BenzylAcetate','Benzaldehyde','MethylButyrate','EthylPropionate','IsobutylAcetate','Pentanedione'})


title(['chemical distance, combination ',num2str(selectedCombination)]);


%% Plot arbitrary combinations

% load('C:\Users\Moritz\Desktop\Odourant Similarity Analysis\new to replace 2nd ones\idcs_local.mat')
% load('C:\Users\Moritz\Desktop\Odourant Similarity Analysis\new to replace 2nd ones\data_zscored.mat')
% data = data_zscored(idcs_local,:);

% classic
these_odours = [46,37,3,10,...
    38,35,1,17,... 
    39,40,29,8];

% option 1
these_odours = [46,37,3,10,...
    38,35,2,17,... 
    39,40,16,8];
% option 2
these_odours = [46,37,3,10,...
    38,35,2,17,... 
    39,40,13,8];
% option 3
these_odours = [46,37,3,10,...
    38,35,16,17,... 
    39,40,13,8];

this_data = data(these_odours,:);

% calculate correlation matrix
these_correlations = zeros(length(these_odours));
for i=1:length(these_odours)
    for j=1:length(these_odours)
        these_correlations(i,j) = corr(this_data(i,:)',this_data(j,:)','Type','Pearson','Rows','Complete');
    end
end
these_distances = 1-these_correlations;
    
% plot heatmap
figure;
heatmap(these_distances);
colormap('pink');

set(gca,'ColorScaling','log')
temp = these_distances(:);
caxis([log(min(temp(temp>0.001))),log(max(these_distances(:)))])

%caxis([0,1.3])

% set(gca,'XData',{'MethylButyrate','EthylPropionate','Pinene','Benzaldehyde','Pentanone','Hexanone','Eucalyptol','Limonene','ButylFormate','PropylFormate','Guaiacol','IsobutylPropionate'})
% set(gca,'XData',{'MethylButyrate','EthylPropionate','Pinene','Benzaldehyde','Pentanone','Hexanone','Eucalyptol','Limonene','ButylFormate','PropylFormate','Guaiacol','IsobutylPropionate'})
% title('Original - Eucalyptol, Guaiacol');
%title('chemical distance');

% set(gca,'XData',{'MethylButyrate','EthylPropionate','Pinene','Benzaldehyde','Pentanone','Hexanone','IsoamylAcetate','Limonene','ButylFormate','PropylFormate','EthylValerate','IsobutylPropionate'})
% set(gca,'XData',{'MethylButyrate','EthylPropionate','Pinene','Benzaldehyde','Pentanone','Hexanone','IsoamylAcetate','Limonene','ButylFormate','PropylFormate','EthylValerate','IsobutylPropionate'})
% title('Option 1 - IsoamylAcetate, EthylValerate');

% set(gca,'XData',{'MethylButyrate','EthylPropionate','Pinene','Benzaldehyde','Pentanone','Hexanone','IsoamylAcetate','Limonene','ButylFormate','PropylFormate','MethylBenzoate','IsobutylPropionate'})
% set(gca,'XData',{'MethylButyrate','EthylPropionate','Pinene','Benzaldehyde','Pentanone','Hexanone','IsoamylAcetate','Limonene','ButylFormate','PropylFormate','MethylBenzoate','IsobutylPropionate'})
% title('Option 2 - IsoamylAcetate, MethylBenzoate');

set(gca,'XData',{'MethylButyrate','EthylPropionate','Pinene','Benzaldehyde','Pentanone','Hexanone','EthylValerate','Limonene','ButylFormate','PropylFormate','MethylBenzoate','IsobutylPropionate'})
set(gca,'XData',{'MethylButyrate','EthylPropionate','Pinene','Benzaldehyde','Pentanone','Hexanone','EthylValerate','Limonene','ButylFormate','PropylFormate','MethylBenzoate','IsobutylPropionate'})
title('Option 3 - EthylValerate, MethylBenzoate');



%% Plot arbitrary combinations

these_odours = [1,6];

this_data = data(these_odours,:);

% calculate correlation matrix
these_correlations = zeros(length(these_odours));
for i=1:length(these_odours)
    for j=1:length(these_odours)
        these_correlations(i,j) = corr(this_data(i,:)',this_data(j,:)','Type','Pearson','Rows','Complete');
    end
end
these_distances = 1-these_correlations;
    
% plot heatmap
figure;
heatmap(these_distances);
colormap('pink');
%caxis([0,0.6])
title('chemical distance');